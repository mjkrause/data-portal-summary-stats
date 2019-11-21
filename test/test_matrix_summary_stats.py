#!/usr/bin/env python3

import unittest
import boto3
import botocore
import os
import shutil
import csv
import responses
import filecmp
from moto import mock_s3
from zipfile import ZipFile
from more_itertools import one
from tempfile import TemporaryDirectory
from src.matrix_summary_stats import MatrixSummaryStats
from src.settings import endpoints, s3_bucket_info


TEST_DIR = os.path.dirname(__file__)
DEPLOYMENT = 'dev'
_BUCKET = f'{DEPLOYMENT}.project-assets.data.humancellatlas.org'  # is source and sink bucket
_SRC_PREFIX = 'project-assets/project-matrices'
AWS_REGION = 'us-west-2'
AWS_ACCESS_KEY_ID = 'fake_access_key'
AWS_SECRET_ACCESS_KEY = 'fake_secret_key'

# Tests were designed with first matrix in mind; second done is larger and can
# help rule out issues caused by how tiny the first one is, but will always fail
# certain tests.
MOCK_MATRIX_ID = 'e7d811e2-832a-4452-85a5-989e2f8267bf' if True else 'dfdd4a14-4ca0-4108-a6f7-52fa2540c61f'
MOCK_MATRIX_FILE = f'{MOCK_MATRIX_ID}.mtx.zip'


class TestMatrixSummaryStatsMatrixServiceResponses(unittest.TestCase):

    def setUp(self) -> None:
        os.chdir(TEST_DIR)
        with ZipFile(MOCK_MATRIX_FILE, 'r') as zipObj:
            zipObj.extractall()

        self.mss = MatrixSummaryStats(
            deployment='prod',
            bucket_name=_BUCKET,
            key=_SRC_PREFIX,
            client='fake_client',
            source_matrix='canned',
            project_field_name='project.provenance.document_id',
            min_gene_count=1200
        )

    def tearDown(self) -> None:
        os.chdir(TEST_DIR)
        shutil.rmtree(f'{MOCK_MATRIX_ID}.mtx/')
        if os.path.isfile('mtx_readme.md'):
            os.remove('mtx_readme.md')
        if self.mss.tmpdir is not None:
            self.mss.cleanup_tmpdir()
        del self.mss

    @responses.activate
    def test_get_project_uuids_from_matrix_service(self):

        # Happy path.
        responses.add(
            responses.GET,
            'https://matrix.data.humancellatlas.org/v1/filters/project.provenance.document_id',
            json={
                'cell_counts': {
                    '08e7b6ba-5825-47e9-be2d-7978533c5f8c': 2,
                    '425efc0a-d3fe-4fab-9f74-3bd829ebdf01': 1,
                    '84e01f35-4f77-4e07-ac7b-058f545d782a': 1
                },
                'field_description': 'Unique identifier for overall project.',
                'field_name': 'project.provenance.document_id',
                'field_type': 'categorical'
            }
        )
        expected_list = ['08e7b6ba-5825-47e9-be2d-7978533c5f8c',
                         '425efc0a-d3fe-4fab-9f74-3bd829ebdf01',
                         '84e01f35-4f77-4e07-ac7b-058f545d782a']
        observed_list = self.mss.get_project_uuids_from_matrix_service()
        self.assertEqual(expected_list, observed_list)

    @responses.activate
    def test_get_expression_matrix_from_service(self):
        # Happy path.

        # URL fixures:
        matrix_service_url = 'https://matrix.data.humancellatlas.org/v1/'
        url_get1 = os.path.join(matrix_service_url, 'filters/project.provenance.document_id')
        url_post = os.path.join(matrix_service_url, 'matrix')
        url_get2 = os.path.join(matrix_service_url + 'matrix/3c44455c-d751-4fc9-a119-3a55c05f8990')
        matrix_url = 'https://s3.amazonaws.com/fake-matrix-service-results/0/424242/13-13.mtx.zip'

        # First response: GET (to satisfy the assertion):
        responses.add(
            responses.GET, url_get1, json={
                'cell_counts': {
                    MOCK_MATRIX_ID: 2,
                },
            }
        )
        # Second response: POST:
        responses.add(
            responses.POST, url_post,
            json={'message': 'Job started.',
                  'non_human_request_ids': {},
                  'request_id': '3c44455c-d751-4fc9-a119-3a55c05f8990',
                  'status': 'In Progress'},
            status=202,
        )

        # Third response: GET returns in progress:
        responses.add(
            responses.GET, url_get2,
            json={'status': 'In Progress'}
        )

        # Fourth response: GET returns "complete":
        responses.add(
            responses.GET, url_get2,
            json={'eta': '',
                  'matrix_url': matrix_url,
                  'message': 'Request 3c44455c-d751-4fc9-a119-3a55c05f8990 has successfully ...',
                  'request_id': '3c44455c-d751-4fc9-a119-3a55c05f8990',
                  'status': 'Complete'},
            status=200
        )

        # Fifth response: GET in get_expression_matrix_from_service
        responses.add(responses.GET, matrix_url, status=200, stream=True)

        self.mss.get_expression_matrix_from_service(projectID=MOCK_MATRIX_ID)

        self.assertEqual(self.mss.matrix_response.status_code, 200)
        self.assertEqual(self.mss.matrix_zipfile_name, '13-13.mtx.zip')

@mock_s3
class TestMatrixSummaryStatsS3(unittest.TestCase):
    # adapted from
    # https://medium.com/@l.peppoloni/how-to-mock-s3-services-in-python-tests-dd5851842946
    def setUp(self) -> None:
        bucket_name = s3_bucket_info[DEPLOYMENT]['bucket_name']
        key = s3_bucket_info[DEPLOYMENT]['key_prefix']

        try:
            s3 = boto3.resource(
                's3',
                region_name=AWS_REGION,
                aws_access_key_id=AWS_SECRET_ACCESS_KEY,
                aws_secret_access_key=AWS_SECRET_ACCESS_KEY
            )
            s3.meta.client.head_bucket(Bucket=_BUCKET)  # source bucket
        except botocore.exceptions.ClientError:
            pass
        else:
            err = f'{_BUCKET} should not exist'
            raise EnvironmentError(err)

        client = boto3.client(
            's3',
            region_name=AWS_REGION,
            aws_access_key_id=AWS_ACCESS_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCESS_KEY
        )
        # Create mock S3 bucket and copy test matrix files to it.
        client.create_bucket(Bucket=_BUCKET)
        self.current_dir = os.path.dirname(__file__)
        _upload_fixtures(bucket=_BUCKET, fixtures_dir=self.current_dir)

        self.mss = MatrixSummaryStats(
            deployment=DEPLOYMENT,
            bucket_name=bucket_name,
            key=key,
            client=client,
            source_matrix='canned',
            project_field_name='project.provenance.document_id',
            min_gene_count=1200
        )

    def tearDown(self) -> None:
        if self.mss.tmpdir is not None:
            self.mss.cleanup_tmpdir()
        del self.mss
        s3 = boto3.resource(
            's3',
            region_name=AWS_REGION,
            aws_access_key_id=AWS_SECRET_ACCESS_KEY,
            aws_secret_access_key=AWS_SECRET_ACCESS_KEY
            )
        bucket = s3.Bucket(_BUCKET)
        for key in bucket.objects.all():
            key.delete()
        bucket.delete()

    def test_get_canned_matrix_filenames(self):
        matrix_files = self.mss.get_canned_matrix_filenames_from_s3()
        expected = [f'project-assets/project-matrices/{MOCK_MATRIX_FILE}']
        self.assertEqual(expected, matrix_files)

    def test_download_canned_expression_matrix_from_S3(self):
        self.assertEqual(_BUCKET, self.mss.s3_bucket_name)
        mtx_file = os.path.join(_SRC_PREFIX, MOCK_MATRIX_FILE)
        downloaded_mtx_files = self.mss.download_canned_expression_matrix_from_s3(mtx_file)
        self.assertTrue(os.path.isdir(self.mss.tmpdir.name))

        expected_mtx_files = [MOCK_MATRIX_FILE]
        self.assertEqual(sorted(expected_mtx_files), sorted(downloaded_mtx_files))

    def test_upload_figs_to_S3(self):
        self.mss.tmpdir = TemporaryDirectory()
        shutil.copyfile(os.path.join(TEST_DIR, MOCK_MATRIX_FILE),
                        os.path.join(self.mss.tmpdir.name, MOCK_MATRIX_FILE))
        self.mss.matrix_response = None
        self.mss.source_matrix = 'canned'
        self.mss.matrix_zipfile_name = MOCK_MATRIX_FILE

        observed = self.mss.write_response_content_to_file_and_unzip()
        expected = {'genes.tsv.gz', 'cells.tsv.gz', 'matrix.mtx.gz'}
        self.assertTrue(set(observed).issuperset(expected))

        self.mss.process_mtx_files()
        self.mss.create_images()
        fig_path = os.path.join(self.mss.tmpdir.name, 'figures')
        figures_dir = os.listdir(fig_path)
        self.assertTrue('filter_genes_dispersion.png' in figures_dir)
        self.assertTrue('highest_expr_genes.png' in figures_dir)

        expected = os.path.join(self.current_dir, 'figures')
        self.assertTrue(filecmp.cmpfiles(expected, fig_path, common='filter_genes_dispersion.png'))
        self.assertTrue(filecmp.cmpfiles(expected, fig_path, common='highest_expr_genes.png'))

        self.mss.upload_figs_to_s3()

        keys = []
        s3 = boto3.resource(
            's3',
            region_name=AWS_REGION,
            aws_access_key_id=AWS_SECRET_ACCESS_KEY,
            aws_secret_access_key=AWS_SECRET_ACCESS_KEY
            )
        bucket = s3.Bucket(_BUCKET)
        for key in bucket.objects.all():
            keys.append(key)

        filter_genes_dispersion_fig = s3.ObjectSummary(
            bucket_name='dev.project-assets.data.humancellatlas.org',
            key=f'project-assets/project-stats/{MOCK_MATRIX_ID}/filter_genes_dispersion.png')
        self.assertTrue(filter_genes_dispersion_fig in keys)

        highest_exp_genes_fig = s3.ObjectSummary(
            bucket_name='dev.project-assets.data.humancellatlas.org',
            key=f'project-assets/project-stats/{MOCK_MATRIX_ID}/highest_expr_genes.png')
        self.assertTrue(highest_exp_genes_fig in keys)


class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        os.chdir(TEST_DIR)
        with ZipFile(MOCK_MATRIX_FILE, 'r') as zipObj:
            zipObj.extractall()

        self.mss = MatrixSummaryStats(
            deployment=DEPLOYMENT,
            bucket_name=_BUCKET,
            key=_SRC_PREFIX,
            client='fake_client',
            source_matrix='canned',
            project_field_name='project.provenance.document_id',
            min_gene_count=1200
        )

    def tearDown(self) -> None:
        os.chdir(TEST_DIR)
        shutil.rmtree(f'{MOCK_MATRIX_ID}.mtx/')
        if os.path.isfile('mtx_readme.md'):
            os.remove('mtx_readme.md')
        if self.mss.tmpdir is not None:
            self.mss.cleanup_tmpdir()
        del self.mss

    def test_deployment_selection(self):
        self.mss.deployment = 'dev'
        expected_url = 'https://matrix.dev.data.humancellatlas.org/v1/'
        self.assertEqual(endpoints[self.mss.deployment]['hca_matrix_service_url'], expected_url)
        expected_url = 'dev.project-assets.data.humancellatlas.org'
        self.assertEqual(s3_bucket_info[self.mss.deployment]['bucket_name'], expected_url)

        self.mss.deployment = 'integration'
        expected_url = 'https://matrix.integration.data.humancellatlas.org/v1/'
        self.assertEqual(endpoints[self.mss.deployment]['hca_matrix_service_url'], expected_url)
        expected_url = 'integration.project-assets.data.humancellatlas.org'
        self.assertEqual(s3_bucket_info[self.mss.deployment]['bucket_name'], expected_url)

        self.mss.deployment = 'staging'
        expected_url = 'https://matrix.staging.data.humancellatlas.org/v1/'
        self.assertEqual(endpoints[self.mss.deployment]['hca_matrix_service_url'], expected_url)
        expected_url = 'staging.project-assets.data.humancellatlas.org'
        self.assertEqual(s3_bucket_info[self.mss.deployment]['bucket_name'], expected_url)

        self.mss.deployment = 'prod'
        expected_url = 'https://matrix.data.humancellatlas.org/v1/'
        self.assertEqual(endpoints[self.mss.deployment]['hca_matrix_service_url'], expected_url)

    def test_make_string_url_compliant_for_search_after(self):
        project_title = 'Single Cell Transcriptomics of a Human Kidney Allograft Biopsy ' \
                        'Defines a Diverse Inflammatory Response'
        document_id = '027c51c6-0719-469f-a7f5-640fe57cbece'
        observed = self.mss.make_string_url_compliant_for_search_after(project_title, document_id)
        expected = 'search_after=Single%20Cell%20Transcriptomics%20of%20a%20Human%20Kidney%' \
                   '20Allograft%20Biopsy%20Defines%20a%20Diverse%20Inflammatory%20Response' \
                   '&search_after_uid=doc%23027c51c6-0719-469f-a7f5-640fe57cbece'
        self.assertEqual(observed, expected)

    def test_process_mtx_files(self):
        # Implicitly tests method "preprocessing".
        self.mss.tmpdir = TemporaryDirectory()
        self.mss.matrix_path = f'{MOCK_MATRIX_ID}.mtx'
        self.mss.process_mtx_files()
        self.assertIsNotNone(self.mss.tmpdir.name)

        with open(self.mss.barcodes_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            first_row = reader.__next__()
        expected_first_row = '00ca0d37-b787-41a4-be59-2aff5b13b0bd'
        self.assertEqual(expected_first_row, one(first_row))  # should only contain one column

        with open(self.mss.genes_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            first_row = reader.__next__()
        expected_first_row = ['ENSG00000000003', 'TSPAN6']
        self.assertEqual(expected_first_row, first_row)

def _upload_fixtures(bucket: str, fixtures_dir: str) -> None:
    client = boto3.client('s3')
    client.upload_file(
        Filename=os.path.join(fixtures_dir, MOCK_MATRIX_FILE),
        Bucket=bucket,
        Key=os.path.join(_SRC_PREFIX, MOCK_MATRIX_FILE)
    )


if __name__ == "__main__":
    unittest.main()
