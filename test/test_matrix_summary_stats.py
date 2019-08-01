#!/usr/bin/env python3

import unittest
import boto3
import os
from src.matrix_summary_stats import MatrixSummaryStats
from src.settings import endpoints, s3_bucket_info



class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        bucket_name = 'foo_bucket'
        key = 'foo_key'
        client = 'foo_client'
        deployment = 'dev'
        source_matrix = 'baz'
        self.mss = MatrixSummaryStats(
            deployment=deployment,
            bucket_name=bucket_name,
            key=key,
            client=client,
            source_matrix=source_matrix
        )

    def tearDown(self) -> None:
        if self.mss.tmpdir is not None:
            self.mss.cleanup_tmpdir()

    def test_get_expression_matrix(self):
        pass

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

    # def test_integration(self):
    #     deployment = 'dev'
    #     bucket_name = s3_bucket_info[deployment]['bucket_name']
    #     key = s3_bucket_info[deployment]['key']
    #     client = boto3.client('s3')
    #     project_field_name = 'project.provenance.document_id'
    #
    #
    #     # FRESH
    #     # source_matrix = 'fresh'
    #     # mss = MatrixSummaryStats(deployment, bucket_name, key, client, source_matrix, project_field_name)
    #     # mss.get_expression_matrix_from_service('2043c65a-1cf8-4828-a656-9e247d4e64f1')
    #     # L = mss.get_project_uuids_from_azul()
    #
    #     # CANNED
    #     source_matrix = 'canned'
    #     mss = MatrixSummaryStats(deployment, bucket_name, key, client, source_matrix, project_field_name)
    #     mtx_files = mss.get_canned_matrix_filenames_from_S3()
    #     print(mtx_files[0])
    #     file = '005d611a-14d5-4fbf-846e-571a1f874f70.mtx.zip'
    #     L = mss.download_canned_expression_matrix_from_S3(file)
    #     f = mss.write_response_content_to_file_and_unzip()
    #
    #     # BOTH
    #     f = mss.write_response_content_to_file_and_unzip()
    #     mss.process_mtx_files()
    #     mss.create_images()
    #     mss.upload_figs_to_s3()


