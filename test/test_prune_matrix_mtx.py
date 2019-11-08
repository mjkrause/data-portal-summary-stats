#/usr/bin/env python3

import unittest
import os
import shutil
import gzip
from tempfile import TemporaryDirectory
from zipfile import ZipFile, is_zipfile
from src.prune_matrix_mtx import prune_matrix_mtx
from src.utils import remove_extension

@unittest.skip
class TestPruneMatrixMtx(unittest.TestCase):
    """
    This test requires a zipped matrix file as returned by the HCA matrix service
    in the test directory. Once a pruned matrix file has been created using "prune_matrix_mtx.py",
    this test should not be executed anymore. Use the unittest.skip decorator to prevent this
    test from being executed.
    """
    def setUp(self) -> None:
        self.src_file = 'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx.zip'
        self.bak_file = 'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx.zip.bak'
        shutil.copyfile(self.src_file, self.bak_file)

    def tearDown(self) -> None:
        os.remove(self.src_file)
        shutil.copyfile(self.bak_file, remove_extension(self.bak_file, 'bak'))
        os.rename(self.bak_file, remove_extension(self.bak_file, 'bak'))

    def test_prune_matrix_mtx(self):

        zipinfo = ZipFile(self.src_file)
        statinfo = os.stat(self.src_file)

        expected_namelist = ['e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/genes.tsv.gz',
                             'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/matrix.mtx.gz',
                             'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/cells.tsv.gz']

        original_namelist = zipinfo.namelist()
        self.assertEqual(sorted(expected_namelist), sorted(original_namelist))
        original_file_size = statinfo.st_size

        percent_prune = 0.97
        os.chdir('..')
        prune_matrix_mtx(self.src_file, percent_prune)

        self.assertTrue(is_zipfile(self.src_file))
        zipinfo = ZipFile(self.src_file)
        statinfo = os.stat(self.src_file)

        expected_namelist = ['e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/genes.tsv.gz',
                             'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/matrix.mtx.gz',
                             'e7d811e2-832a-4452-85a5-989e2f8267bf.mtx/cells.tsv.gz']
        pruned_namelist = zipinfo.namelist()
        self.assertEqual(sorted(expected_namelist), sorted(pruned_namelist))
        pruned_file_size = statinfo.st_size

        # Test whether pruned file has been reduced in size by at least 75% (0.25 * original size).
        self.assertGreaterEqual(0.25, pruned_file_size/original_file_size)

        currdir = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            src_file = os.path.join(currdir, self.src_file)
            dst_file = os.path.join(tmpdir, self.src_file)
            shutil.copyfile(src_file, dst_file)
            with ZipFile(dst_file, 'r') as zipObj:
                zipObj.extractall(tmpdir)
            dst_file = remove_extension(dst_file, 'zip')
            matrix_file = os.path.join(dst_file, 'matrix.mtx.gz')
            # Test whether the number of lines reported in header of matrix.mtx is correct.
            second_header_line = ''
            observed_row_count = 0
            for line in gzip.open(matrix_file).readlines():
                if observed_row_count == 1:
                    second_header_line = line.decode('utf-8')
                observed_row_count += 1
            observed_second_header_row = second_header_line.strip('\r\n').split(' ')
            expected_second_header_row = ['63925', '2544', str(observed_row_count - 2)]
            self.assertEqual(expected_second_header_row, observed_second_header_row)

if __name__ == '__main__':
    unittest.main()
