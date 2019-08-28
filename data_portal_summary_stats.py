#/usr/env/python3

import sys
import os
import argparse
import logging
from src.runner import run_data_portal_summary_stats


# Set up logging
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=f'{ROOT_DIR}/{__file__.strip(".py")}.log',
                    filemode='w')
# These libraries make a lot of debug-level log messages which make the log file hard to read
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def main():
    parser = argparse.ArgumentParser(
        description='Per-project summary statistics and figures.'
    )
    args_group = parser.add_argument_group(title='arguments')
    args_group.add_argument(
        '--environ',
        default='dev',
        choices=('dev', 'integration', 'staging', 'prod'),
        type=str,
        help='Deployment environment (default: "dev") from which '
             'matrix data are requested to create summary '
             'statistics.'
    )
    args_group.add_argument(
        '--source',
        default='fresh',
        choices=('fresh', 'canned'),
        type=str,
        help='Source of matrix files. "fresh" (default) denotes requesting '
             'matrix files from the matrix service. "canned" denotes '
             'downloading already created matrix files from AWS S3.'
    )
    args_group.add_argument(
        '--blacklist',
        default='false',
        choices=('false', 'true'),
        type=str,
        help='Skip files with project IDs listed in a file named '
             '"blacklist" during processing. Default: false.'
    )
    args_group.add_argument(
        '--min_cell_count',
        default=1000,
        choices=range(300, 2000),
        metavar="[300-2000]",
        type=int,
        help='The minimal cell count in the "genes_detected" field '
             '(range from 300 to 2000). Default is 1000.'
    )
    args_group.set_defaults(environ='dev', source='fresh', blacklist='False', min_cell_count=1000)
    args = parser.parse_args()

    run_data_portal_summary_stats(args)


if __name__ == "__main__":
    main()
