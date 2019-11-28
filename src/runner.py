#!/usr/bin/env python3

import argparse
import logging
import time

from src import settings
from src.matrix_preparer import MatrixPreparer
from src.matrix_provider import (
    FreshMatrixProvider,
    CannedMatrixProvider,
)
from src.matrix_summary_stats import MatrixSummaryStats
from src.s3_service import S3Service
from src.utils import TemporaryDirectoryChange

log = logging.getLogger(__name__)


def run_data_portal_summary_stats(args: argparse.Namespace):
    log.info(f'\nGenerating per-project summary statistics of matrix data from '
             f'{args.environ} deployment environment.\n')

    s3_bucket_info = settings.s3_bucket_info[args.envrion]
    endpoints = settings.endpoints[args.environ]

    s3 = S3Service(s3_bucket_info['bucket_name'], s3_bucket_info['key_prefix'])

    # Temporary work-around for matrices that can't be processed for various reasons.
    do_not_process = s3.get_blacklist() if args.blacklist else []

    if args.source == 'fresh':
        provider = FreshMatrixProvider(blacklist=do_not_process,
                                       min_gene_count=args.min_gene_count,
                                       project_field_name='project.provenance.document_id',
                                       hca_endpoint=endpoints['hca_matrix_service_url'],
                                       azul_endpoint=endpoints['azul_projects_url'])
    elif args.source == 'canned':
        provider = CannedMatrixProvider(blacklist=do_not_process,
                                        s3_service=s3)
    else:
        assert False

    iter_matrices = iter(provider)
    while True:
        with TemporaryDirectoryChange() as tempdir:
            try:
                mtx_info = next(iter_matrices)
            except StopIteration:
                break
            log.info(f'Writing to temporary directory {tempdir}')
            log.info(f'Processing matrix for project {mtx_info.project_uuid} ({mtx_info.source})')

            preparer = MatrixPreparer(mtx_info)
            preparer.unzip()
            preparer.preprocess()

            for pp_mtx_info in preparer.separate():
                log.info(f'Generating stats for {pp_mtx_info.extract_path}')
                mss = MatrixSummaryStats(pp_mtx_info)
                mss.create_images()
                s3.upload_figures(mtx_info)
                log.info('...done uploading figures')

                # This logic was in Krause's code, no idea why
                if mtx_info.source == 'fresh':
                    time.sleep(15)
