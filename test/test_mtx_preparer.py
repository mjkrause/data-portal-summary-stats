import logging

import os
import unittest

import pandas as pd
import scanpy as sc

from src.matrix_preparer import MatrixPreparer
from test.tempdir_test_case import (
    MockMatrixTestCase,
)

log = logging.getLogger(__name__)


class TestMatrixPreparer(MockMatrixTestCase):

    def setUp(self):
        super().setUp()
        self.preparer = MatrixPreparer(self.info)

    def test_strip_version(self):
        test_cases = [
            # real cases
            ('optimus_v1.3.1', 'optimus'),
            ('smartseq2_v2.3.0', 'smartseq2')
        ]

        for versioned, unversioned in test_cases:
            self.assertEqual(self.preparer.strip_version_suffix(versioned), unversioned)

    def test_unzip(self):
        self.preparer.unzip()

        self.assertFalse(os.path.exists(self.info.zip_path))
        self.assertEqual(set(os.listdir(self.info.extract_path)), set(MatrixPreparer.unproc_files.values()))

    def test_prepprocess(self):

        self.preparer.unzip()
        self.preparer.preprocess()
        self.assertEqual(set(os.listdir(self.info.extract_path)), set(MatrixPreparer.proc_files.values()))

        for filekey in ['barcodes', 'genes']:
            filename = os.path.join(self.info.extract_path, MatrixPreparer.proc_files[filekey])
            df = pd.read_csv(filename, sep='\t', index_col=None, header=None)
            #  Unfortuantely Scanpy requires us to remove the headers so we lose
            #  a lot of info that would be useful in verification
            self.assertEqual(len(df.columns), len(MatrixPreparer.file_columns[filekey]))

    def test_separate_homogeneous(self):
        self.preparer.unzip()
        self.preparer.preprocess()
        # mock matrix has only two versions of hte same protocol so it will not
        # be separated here.
        sep_infos = self.preparer.separate(strip_version=True)
        self.assertEqual(len(sep_infos), 1)
        self.assertEqual(sep_infos[0], self.info)
        self.assertEqual(self.info.lib_con_method, 'smartseq2')

        sc.read_10x_mtx(self.info.extract_path)

    def test_separate_heterogeneous(self):

        expected_libcon_methods = {'smartseq2_v2.3.0', 'smartseq2_v2.4.0'}

        self.preparer.unzip()
        self.preparer.preprocess()
        # likewise, this will keep protocol versions separate from each other
        sep_infos = self.preparer.separate(strip_version=False)

        observed_libcon_methods = {i.lib_con_method for i in sep_infos}

        self.assertEqual(observed_libcon_methods, expected_libcon_methods)

        for sep_info in sep_infos:
            self.assertTrue(os.path.isdir(sep_info.extract_path))
            self.assertEqual(sep_info.extract_path, os.path.join(self.info.extract_path, sep_info.lib_con_method))
            for filekey, filename in MatrixPreparer.proc_files.items():
                filename = os.path.join(sep_info.extract_path, filename)
                if filekey == 'genes':
                    self.assertTrue(os.path.islink(filename))
                else:
                    self.assertTrue(os.path.isfile(filename))
            sc.read_10x_mtx(sep_info.extract_path)


if __name__ == '__main__':
    unittest.main()
