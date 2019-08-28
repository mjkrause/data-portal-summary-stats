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
import urllib
import matplotlib
from zipfile import ZipFile
from typing import List, Dict, Any
from tempfile import TemporaryDirectory
from botocore.exceptions import ClientError
from more_itertools import first
from src.settings import endpoints
from src.utils import convert_size


logger = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)

# See https://stackoverflow.com/questions/27147300/
# matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
# for why the following line is required.
matplotlib.use('Agg')


class MatrixSummaryStats:

    def __init__(self, deployment: str, bucket_name: str, key: str, client: str,
                 source_matrix: str, project_field_name, min_cell_count: int):
        self.deployment = deployment
        self.s3_bucket_name = bucket_name
        self.s3_key = key
        self.client = client
        self.source_matrix = source_matrix
        self.project_field_name = project_field_name  # needed if matrices requests from service
        self.min_cell_count = min_cell_count  # for field genes_detected in matrix filter
        self.projdir = os.getcwd()
        self.project_uuid = None
        self.tmpdir = None
        self.matrix_zipfile_name = None
        self.matrix_response = None
        self.matrix_path = None  # rel. path to unzipped directory
        self.genes_path = None
        self.barcodes_path = None

    def get_project_uuids_from_matrix_service(self) -> list:
        """Return list matrix directory names (with prefix keys) from matrix service"""
        response = requests.get(endpoints[self.deployment]['hca_matrix_service_url'] +
                                'filters' + '/' + self.project_field_name)
        self.check_response(response)

        return list(response.json()['cell_counts'].keys())

    def get_canned_matrix_filenames_from_S3(self) -> list:
        """Return list of canned matrix directory names (with prefix keys)
        contained in S3 bucket."""
        prefix = 'project-assets/project-matrices'
        postfix = '.mtx.zip'
        matrix_files = []
        try:
            response = self.client.list_objects_v2(Bucket=self.s3_bucket_name)
            for obj in response['Contents']:
                if obj['Key'].startswith(prefix) and obj['Key'].endswith(postfix):
                    matrix_files.append(obj['Key'])
        except ClientError as cerr:
            logger.info(cerr)

        return matrix_files

    def get_project_uuids_from_azul(self) -> List[Dict[str, Any]]:
        """Get all project UUIDs from Azul by using the service APIs search_after query parameter."""
        url = endpoints[self.deployment]['azul_projects_url']
        project_uuids = []
        search_after_string = ''
        first_iter = True
        while True:
            if first_iter:
                response = requests.get(url)
                first_iter = False
            else:
                response = requests.get(url + '?' + search_after_string)
            try:
                hits = response.json()['hits']
            except KeyError:

                return project_uuids

            search_after_string = self.make_string_url_compliant_for_search_after(
                project_title = first(hits[-1]['projects'])['projectTitle'],
                document_id = hits[-1]['entryId'])
            project_uuids = project_uuids + [
                {'project_title': first(hits[i]['projects'])['projectTitle'],
                 'project_UUID': hits[i]['entryId']} for i in range(len(hits))
            ]

    def get_expression_matrix_from_service(self, projectID=None) -> None:
        self.tmpdir = TemporaryDirectory()
        status_response = self._request_matrix(projectID)
        assert status_response.status_code == 200
        s3_download_url = status_response.json()['matrix_url']
        logger.info(f'Download URL for matrix is {s3_download_url}')
        self.matrix_response = requests.get(s3_download_url, stream=True)
        self.matrix_zipfile_name = os.path.basename(s3_download_url)

    def download_canned_expression_matrix_from_S3(self, mtx_file) -> list:
        """Download matrix directory into local temporary directory and return as list."""
        self.tmpdir = TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        self.matrix_zipfile_name = os.path.basename(mtx_file)
        self.client.download_file(Bucket=self.s3_bucket_name,
                                  Key=mtx_file,
                                  Filename=self.tmpdir.name + '/' + self.matrix_zipfile_name)
        assert self.matrix_zipfile_name in os.listdir('.')  # confirm successful download
        size = os.path.getsize(self.matrix_zipfile_name)    # in bytes
        logger.info(f'Size of {self.matrix_zipfile_name}: {convert_size(size)}')

        return os.listdir(self.tmpdir.name)

    def write_response_content_to_file_and_unzip(self) -> list:
        """Return list of canned matrix directory names (with prefix keys)
        contained in S3 bucket."""
        logger.info('Decompressing files...')
        os.chdir(self.tmpdir.name)
        logger.info(f'Writing to temporary directory {self.tmpdir.name}')
        if self.matrix_response:  # only if matrix is from matrix service
            with open(self.matrix_zipfile_name, 'wb') as matrix_zip_file:
                shutil.copyfileobj(self.matrix_response.raw, matrix_zip_file)
        size = os.path.getsize(self.matrix_zipfile_name)  # in bytes
        logger.info(f'Size of {self.matrix_zipfile_name}: {convert_size(size)}')
        self.get_matrix_path_and_project_uuid()
        self._unzip_files()
        os.chdir(self.matrix_path)
        mtx_files = os.listdir('.')
        if self.source_matrix == 'fresh':
            assert len(mtx_files) == 3
        elif self.source_matrix == 'canned':
            assert len(mtx_files) >= 3 and len(mtx_files) <= 4  # may contain extra readme file
        assert 'cells.tsv.gz' in mtx_files
        assert 'genes.tsv.gz' in mtx_files
        assert 'matrix.mtx.gz' in mtx_files
        os.chdir(self.tmpdir.name)
        os.remove(self.matrix_zipfile_name)

        return mtx_files

    def _unzip_files(self):
        """Unzip each file in zip_file individually and remove source immediately."""
        os.mkdir(os.path.join(self.tmpdir.name, self.matrix_path))
        with ZipFile(self.matrix_zipfile_name) as zip_file:
            for member in zip_file.namelist():
                filename = os.path.basename(member)
                if not filename:
                    continue  # skip directories
                source = zip_file.open(member)
                target = open(os.path.join(self.tmpdir.name, self.matrix_path, filename), 'wb')
                with source, target:
                    shutil.copyfileobj(source, target)

    def create_images(self) -> None:
        # ToDo: parameterize..., and add violin plot
        logger.info('Creating figures...')
        figure_format = '.png'
        logger.info(f'Figures saved in {figure_format} format.')
        logger.info(f'Path to matrix files is {self.matrix_path}')

        # 1. Figure: highest-expressing genes.
        adata = sc.read_10x_mtx(self.matrix_path, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        sc.pl.highest_expr_genes(adata, n_top=20, save=figure_format, show=False)  # write to disk

        # 2. Figure: highest-variable genes:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e3)
        sc.pp.log1p(adata)  # logarithmize
        adata.raw = adata   # save raw data
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata, save=figure_format, show=False)  # write to disk

    def upload_figs_to_s3(self) -> None:
        os.chdir(self.tmpdir.name)
        figures = os.listdir('figures/')
        if figures is None:
            return
        else:
            for figure in figures:
                key = self.s3_key + self.project_uuid + '/' + figure
                logger.info(f'Uploading {figure} to S3 bucket {self.s3_bucket_name} into {key}')
                self.client.upload_file(Filename=f'figures/{figure}',
                                        Bucket=self.s3_bucket_name,
                                        Key=key)
        logger.info('...done uploading figures')
        self.cleanup_tmpdir()
        logger.info('Temporary directory removed')

    def _request_matrix(self, project_document_id: str = None) -> requests.models.Response:
        # ToDo: how to find project UUIDs?
        # Parameters.
        feature = 'gene'
        format_ = 'mtx'
        project_field_name = self.project_field_name
        min_cell_count = self.min_cell_count
        min_cell_count_field = 'genes_detected'

        hca_matrix_service_url = endpoints[self.deployment]['hca_matrix_service_url']

        list_projects_url = hca_matrix_service_url + 'filters' + '/' + project_field_name
        assert project_document_id in requests.get(list_projects_url).json()['cell_counts']
        self.project_uuid = project_document_id

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
        logger.info(f'Requesting expression matrix for project document ID {project_document_id}')
        logger.info(f'Request payload and filter settings: {payload}')
        response = requests.post(hca_matrix_service_url + 'matrix', json=payload)
        self.check_response(response)
        minute_counter = 0
        while True:
            status_response = requests.get(hca_matrix_service_url + 'matrix' + '/' +
                                           response.json()['request_id'])
            if status_response.json()['status'] == 'Complete':
                break
            elif status_response.json()['status'] == 'In Progress':
                logger.info(f'Matrix request status: {status_response.json()["status"]}...')
                time.sleep(30)
                minute_counter += 0.5
            else:
                sys.exit(f'Matrix Service request status is: {status_response.json()["status"]}')
        logger.info(f'Successfully requested matrix in {minute_counter} min.')

        return status_response

    def process_mtx_files(self):
        logger.info('Preprocessing matrix files...')
        os.chdir(self.matrix_path)
        mtx_files = os.listdir('.')
        for mtx_file in mtx_files:
            logger.info(f'Processing file {mtx_file}.')
            if mtx_file.endswith('gz'):
                with gzip.open(mtx_file, 'rb') as f_in:
                    if f_in.name == 'genes.tsv.gz' or f_in.name == 'cells.tsv.gz':
                        self.preprocessing(f_in)  # writes files to disk
                    elif f_in.name == 'matrix.mtx.gz':
                        outfile = first(os.path.splitext(f_in.name))
                        with open(outfile, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
            os.remove(mtx_file)
        files = os.listdir('.')
        assert 'barcodes.tsv' in files
        assert 'genes.tsv' in files
        logger.info(f'Processed files {files}')
        self.genes_path = os.path.abspath('genes.tsv')
        self.barcodes_path = os.path.abspath('barcodes.tsv')
        os.chdir(self.tmpdir.name)

    @staticmethod
    def preprocessing(fileobj: gzip.GzipFile) -> None:
        """Remove headers and columns to prepare TSV files in mtx directory to make them
        compatible for use with ScanPy methods."""
        f = pd.read_table(fileobj, sep='\t')  # Pandas dataframe
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

    @staticmethod
    def check_response(response: requests.Response) -> None:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            logger.info(f'{str(err)}')

    @staticmethod
    def make_string_url_compliant_for_search_after(project_title: str, document_id: str) -> str:
        """Return input string to be URL compliant, i.e., replacing special characters by their
         corresponding hexadecimal representation."""
        return 'search_after=' + \
               urllib.parse.quote(project_title) + \
               '&search_after_uid=doc' + \
               urllib.parse.quote('#' + document_id)

    def cleanup_tmpdir(self):
        self.tmpdir.cleanup()

    def get_matrix_path_and_project_uuid(self):
        """Create name of project ID with .mtx extension."""
        try:
            # For some reason the mtx directory is named by the project UUID for canned files.
            if self.source_matrix == 'canned':
                self.project_uuid = first(self.matrix_zipfile_name.split('.'))
            self.matrix_path = self.project_uuid + '.mtx'
        except ValueError:
            self.matrix_path = []

