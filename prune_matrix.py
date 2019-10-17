#/usr/bin/env python3

import sys
import os
from src.prune_matrix_mtx import prune_matrix_mtx

def main():
    """
    Utility to reduce the size of an mtx matrix market file by some percentage by randomly
    excluding rows in "matrix.mtx" from the original file. Run this method if you need an HCA
    matrix file in matrix market format for testing.
    The method removes randomly rows from "matrix.mtx". The magnitude of that process is controlled
    by parameter percent_prune. The utility expects a zipped matrix market file downloaded
    from the HCA matrix service in the test directory.

    WARNING: this is an in-place operation and replaces the input file.

    :parameter mtx_dir_name:   original, zipped "mtx" file that contains "genes.tsv.gz",
                               "cells.tsv.gz", and "matrix.mtx.gz" in a directory of the same
                               name as the zip-file.
    :parameter percent_prune:  percentage of pruning, float in range [0 .. 1); the larger
                               this parameter the smaller the output file
                               (e.g., 0.95 leaves ~25% of the original file).

    Run like so (use correct Zipfile name and percent-accepted value):
        $ python src/prune_matrix_mtx.py test/e7d811e2-832a-4452-85a5-989e2f8267bf.mtx.zip 0.97
    from project root.
    """

    mtx_dir_name = sys.argv[1]
    if len(sys.argv) > 2:
        percent_prune = float(sys.argv[2])
    else:
        percent_prune = None
    mtx_dir_name = os.path.basename(mtx_dir_name)
    prune_matrix_mtx(mtx_dir_name=mtx_dir_name, percent_prune=percent_prune)


if __name__ == '__main__':
    main()