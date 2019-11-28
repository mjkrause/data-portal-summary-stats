#!/usr/bin/env python3

import unittest
import os
import filecmp

from more_itertools import first

from src.matrix_preparer import MatrixPreparer
from src.matrix_summary_stats import MatrixSummaryStats
from test.tempdir_test_case import MockMatrixTestCase


class TestMatrixSummaryStats(MockMatrixTestCase):

    def setUp(self) -> None:
        super().setUp()
        preparer = MatrixPreparer(self.info)
        preparer.unzip()
        preparer.preprocess()
        new_info = first(preparer.separate())
        self.mss = MatrixSummaryStats(new_info)

    def test_create_images(self):
        self.mss.create_images()
        fig_path = 'figures'
        figures_dir = os.listdir(fig_path)
        self.assertTrue('filter_genes_dispersion.png' in figures_dir)
        self.assertTrue('highest_expr_genes.png' in figures_dir)

        expected = 'figures'
        self.assertTrue(filecmp.cmpfiles(expected, fig_path, common='filter_genes_dispersion.png'))
        self.assertTrue(filecmp.cmpfiles(expected, fig_path, common='highest_expr_genes.png'))

        # TODO test more figures


if __name__ == "__main__":
    unittest.main()
