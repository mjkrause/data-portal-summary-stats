#!/usr/bin/env python3

import sys
import os
import time
import shutil
import gzip
import pandas as pd
import scanpy as sc
import requests
import logging
from zipfile import ZipFile
from typing import List, Dict, Any
from tempfile import TemporaryDirectory
from more_itertools import first
from src.settings import endpoints


logger = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)


class MatrixSummaryStats:

    def __init__(self, bucket_name, key, client):
        self.client = client
        self.s3_bucket_name = bucket_name
        self.s3_key = key
        self.project_UUIDs = self.get_project_UUIDs()
        self.project_UUID = None
        self.projdir = os.getcwd()
        self.project_title = None
        self.tmpdir = TemporaryDirectory()
        self.matrix_zipfile_name = None
        self.matrix_response = None
        self.matrix_path = None
        self.genes_path = None
        self.barcodes_path = None

    def get_expression_matrix(self) -> None:
        status_response = self._request_matrix()
        assert status_response.status_code == 200
        s3_download_url = status_response.json()['matrix_url']
        logger.info(f'Download URL for matrix is {s3_download_url}')
        self.matrix_response = requests.get(s3_download_url, stream=True)
        self.matrix_zipfile_name = os.path.basename(s3_download_url)

    def unzip_files(self) -> None:
        logger.info('Decompressing files...')
        os.chdir(self.tmpdir.name)
        logger.info(f'Writing to temporary directory {self.tmpdir.name}')
        with open(self.matrix_zipfile_name, 'wb') as matrix_zip_file:
            shutil.copyfileobj(self.matrix_response.raw, matrix_zip_file)
        ZipFile(self.matrix_zipfile_name).extractall()
        self.matrix_path = first(os.path.splitext(self.matrix_zipfile_name))  # remove ".zip"
        os.chdir(self.matrix_path)
        files = os.listdir('.')

        assert len(files) == 3
        assert 'cells.tsv.gz' in files
        assert 'genes.tsv.gz' in files
        assert 'matrix.mtx.gz' in files

        logger.info('Now preprocessing matrix...')
        for file in files:
            logger.info(f'Processing file {file}.')
            with gzip.open(file, 'rb') as f_in:
                if f_in.name == 'genes.tsv.gz' or f_in.name == 'cells.tsv.gz':
                    self.preprocessing(f_in)  # writes files to disk
                elif f_in.name == 'matrix.mtx.gz':
                    outfile = first(os.path.splitext(f_in.name))
                    with open(outfile, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            os.remove(file)
        files = os.listdir('.')
        assert 'barcodes.tsv' in files
        assert 'genes.tsv' in files
        logger.info(f'Processed files {files}')
        self.genes_path = os.path.abspath('genes.tsv')
        self.barcodes_path = os.path.abspath('barcodes.tsv')
        os.chdir(self.tmpdir.name)

    def create_images(self) -> None:
        # ToDo: parameterize..., and add violin plot
        logger.info('Creating figures...')
        logger.info(f'Path to matrix files is {self.matrix_path}')

        # 1. Figure: highest-expressing genes.
        adata = sc.read_10x_mtx(self.matrix_path, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        sc.pl.highest_expr_genes(adata, n_top=20, save='.png', show=False)  # write to disk

        # 2. Figure: highest-variable genes:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e3)
        sc.pp.log1p(adata)  # logarithmize
        adata.raw = adata   # save raw data
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata, save='.png', show=False)  # write to disk

    def upload_figs_to_s3(self):
        os.chdir(self.tmpdir.name)
        figures = os.listdir('figures/')
        if figures is None:
            return
        else:
            for figure in figures:
                key = self.s3_key + self.project_UUID + '/' + figure
                logger.info(f'Uploading {figure} to S3 bucket {self.s3_bucket_name} into {key}')
                self.client.upload_file(Filename=f'figures/{figure}',
                                        Bucket=self.s3_bucket_name,
                                        Key=key)
        logger.info('...done uploading figures')
        self.tmpdir.cleanup()
        logger.info('Temporary directory destroyed')

    def _request_matrix(self) -> requests.models.Response:
        # ToDo: how to find project UUIDs?
        # Parameters.
        feature = 'gene'
        format_ = 'mtx'
        project_title = 'Single cell transcriptome analysis of human pancreas'  # get single project
        project_document_id = 'cddab57b-6868-4be4-806f-395ed9dd635a'
        project_field_name = 'project.provenance.document_id'
        #project_field_name = 'project.project_core.project_short_name'
        min_cell_count = 300
        min_cell_count_field = 'genes_detected'

        hca_matrix_service_url = endpoints['hca_matrix_service_url']

        list_projects_url = hca_matrix_service_url + '/filters/' + project_field_name
        assert project_document_id in requests.get(list_projects_url).json()['cell_counts']
        self.project_title = project_title
        self.project_UUID = project_document_id

        payload = {
            'feature': feature,
            'format': format_,
            'filter': {
                'op': 'and',
                'value': [
                    {
                        'op': '=',
                        'value': project_document_id,
                        'field': project_field_name
                    },
                    {
                        'op': '>=',
                        'value': min_cell_count,
                        'field': min_cell_count_field
                    }
                ]
            }
        }
        logger.info(f'Requesting expression matrix for project {project_title}...')
        response = requests.post(hca_matrix_service_url + '/matrix', json=payload)
        self.check_response(response)

        while True:
            status_response = requests.get(hca_matrix_service_url + '/matrix/' +
                                           response.json()['request_id'])
            if status_response.json()['status'] == 'Complete':
                break
            elif status_response.json()['status'] == 'In Progress':
                logger.info(f'Matrix request status: {status_response.json()["status"]}...')
                time.sleep(30)
            else:
                sys.exit(f'Matrix Service request status is: {status_response.json()["status"]}')
        logger.info('Successfully requested matrix.')

        return status_response

    @staticmethod
    def preprocessing(fileobj: gzip.GzipFile) -> None:
        """Remove headers and columns to prepare TSV files in mtx directory to make them
        compatible for use with ScanPy methods."""
        f = pd.read_table(fileobj, sep='\t')  # Pandas dataframe
        col_to_keep = []
        if fileobj.name == 'genes.tsv.gz':
            col_to_keep = ['featurekey', 'featurename']
            assert col_to_keep[0] in f.columns
            assert col_to_keep[1] in f.columns
        elif fileobj.name == 'cells.tsv.gz':
            fileobj.name = 'barcodes.tsv.gz'
            col_to_keep = 'cellkey'
            assert col_to_keep in f.columns
        else:
            raise ValueError('Expected genes.tsv.gz and cells.tsv.gz in directory.')
        f_new = f[col_to_keep]
        # Write to file without column or row headers.
        f_new.to_csv(first(os.path.splitext(fileobj.name)), index=False, header=False, sep='\t')

    def get_project_UUIDs(self) -> List[Dict[str,Any]]:
        response = requests.get(endpoints['azul_projects_url'])
        self.check_response(response)
        hits = response.json()['hits']
        logger.info(f"Found {len(hits)} projects from Azul's /repository/projects endpoint.")
        return [{'project_title': first(hits[i]['projects'])['projectTitle'],
                 'project_UUID': hits[i]['entryId']} for i in range(len(hits))]

    @staticmethod
    def check_response(response: requests.Response) -> None:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            logger.info(f'{str(err)}')
