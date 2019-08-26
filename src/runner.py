#!/usr/bin/env python3

import os
import sys
import boto3
import time
import argparse
import logging
from more_itertools import first
from src.settings import s3_bucket_info
from src.matrix_summary_stats import MatrixSummaryStats
from src.utils import get_blacklist_from_s3


logger = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)

_MSS = None  # mss is matrix summary stats, create global instance
_MSS_CLIENT = None  # boto3 client


def run_data_portal_summary_stats(args: argparse.Namespace):
    input_args = {
        'deployment': args.environ,
        'project_field_name': 'project.provenance.document_id',
        'source_of_matrix': args.source
    }
    # Temporary work-around for matrices that can't be processed for various reasons.
    if args.blacklist == 'false':
        do_not_process = []
    else:
        assert args.blacklist == 'true'
        client = mss_client()
        do_not_process = get_blacklist_from_s3(
            client=client,
            bucket=s3_bucket_info[args.environ]['bucket_name'],
            key='blacklist'
        )

    logger.info(f'\nGenerating per-project summary statistics of matrix data from '
                f'{get_mss(**input_args).deployment} deployment environment.\n')
    if args.source == 'fresh':
        project_IDs = get_mss(**input_args).get_project_uuids_from_matrix_service()
        project_info_azul = get_mss(**input_args).get_project_uuids_from_azul()
        for project_ID in project_IDs:
            logger.info(f'\nProcessing matrix from service with project ID {project_ID}')
            if project_ID in do_not_process:
                logger.info(f'   ...skipping project ID {project_ID}')
                continue
            project_title = get_project_title(project_info_azul, project_ID)
            if project_title:
                logger.info(f'Project title: {project_title}')
            else:
                logger.info(f'No project title found for project ID {project_ID} in Azul')
            get_mss(**input_args).get_expression_matrix_from_service(project_ID)
            get_mss(**input_args).write_response_content_to_file_and_unzip()
            get_mss(**input_args).process_mtx_files()
            get_mss(**input_args).create_images()
            get_mss(**input_args).upload_figs_to_s3()
            time.sleep(15)
    elif args.source == 'canned':
        mtx_files = get_mss(**input_args).get_canned_matrix_filenames_from_S3()
        for mtx_file in mtx_files:
            zipfile_ID = first(os.path.basename(mtx_file).split('.'))
            if zipfile_ID in do_not_process:
                logger.info(f'\nSkipping processing of project ID {zipfile_ID}')
                continue
            logger.info(f'\nDownloading canned matrix file with project ID {zipfile_ID} from S3.')
            canned_file = get_mss(**input_args).download_canned_expression_matrix_from_S3(mtx_file)
            if canned_file == []:
                continue
            logger.info(f'Processing matrix with project ID {zipfile_ID} from S3')
            get_mss(**input_args).write_response_content_to_file_and_unzip()
            get_mss(**input_args).process_mtx_files()
            get_mss(**input_args).create_images()
            get_mss(**input_args).upload_figs_to_s3()


def mss_client():
    global _MSS_CLIENT
    if _MSS_CLIENT is None:
        _MSS_CLIENT = boto3.client('s3')
    return _MSS_CLIENT


def get_mss(**input_args):
    global _MSS
    if _MSS is None:
        _MSS = MatrixSummaryStats(
            deployment=input_args['deployment'],
            bucket_name=s3_bucket_info[input_args['deployment']]['bucket_name'],
            key=s3_bucket_info[input_args['deployment']]['key'],
            client=mss_client(),
            source_matrix=input_args['source_of_matrix'],
            project_field_name=input_args['project_field_name']
        )
    return _MSS


def get_project_title(proj_infos: list, project_id: str) -> str:
    for proj_info in proj_infos:
        for key, val in proj_info.items():
            if val == project_id:
                return proj_info['project_title']