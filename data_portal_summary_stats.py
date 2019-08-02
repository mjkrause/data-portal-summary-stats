#/usr/env/python3

import sys
import os
import logging
import boto3
from src.settings import s3_bucket_info
from src.matrix_summary_stats import MatrixSummaryStats


# Set up logging
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=f'{ROOT_DIR}/{__file__}.log',
                    filemode='w')
# These libraries make a lot of debug-level log messages which make the log file hard to read
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)

_MSS = None  # mss is matrix summary stats, create global instance
_MSS_CLIENT = None  # boto3 client


def main():
    logger.info('Downloading matrix, ETF, uploading images to S3...')
    get_mss().get_expression_matrix()
    get_mss().unzip_files()
    get_mss().create_images()
    get_mss().upload_figs_to_s3()


def mss_client():
    global _MSS_CLIENT
    if _MSS_CLIENT is None:
        _MSS_CLIENT = boto3.client('s3')
    return _MSS_CLIENT


def get_mss():
    global _MSS
    if _MSS is None:
        _MSS = MatrixSummaryStats(s3_bucket_info['bucket_name'], s3_bucket_info['key'], mss_client())
    return _MSS


if __name__ == "__main__":
    main()