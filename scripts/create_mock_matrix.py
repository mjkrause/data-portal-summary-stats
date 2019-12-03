#!./.venv/bin/python
import logging
import os

from src import Config
from src.matrix_preparer import MatrixPreparer
from src.matrix_provider import CannedMatrixProvider
from src.s3_service import S3Service
from src.utils import DirectoryChange
from test.tempdir_test_case import MockMatrixTestCase

log = logging.getLogger(__name__)


def main():
    os.environ['AWS_DEFAULT_PROFILE'] = 'hca-id'
    s3 = S3Service(Config('dev'))
    provider = CannedMatrixProvider(s3_service=s3)
    with DirectoryChange('../test/'):
        log.info('Entering test directory')
        mtxinfo = next(iter(provider))
        log.info(f'Obtained matrix: {mtxinfo.zip_path}')
        prep = MatrixPreparer(mtxinfo)
        log.info('Processing...')
        prep.unzip()
        prep.prune(0.05)
        prep.rezip(remove_dir=True, zip_path=MockMatrixTestCase.mock_matrix)
        log.info('Finished.')


if __name__ == '__main__':
    main()
