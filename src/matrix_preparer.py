from collections import defaultdict
from functools import lru_cache
import gzip
import logging
import os

from more_itertools import (
    first,
    one,
)
import numpy as np
import pandas as pd
import shutil
import re
from typing import (
    List,
    Optional,
)
from zipfile import ZipFile

from src.matrix_info import MatrixInfo
from src.utils import DirectoryChange

log = logging.getLogger(__name__)


class MatrixPreparer:
    library_method_column = 'analysis_protocol.protocol_core.protocol_id'

    unproc_files = {
        'genes': 'genes.tsv.gz',
        'barcodes': 'cells.tsv.gz',
        'matrix': 'matrix.mtx.gz'
    }

    proc_files = {
        'genes': 'genes.tsv',
        'barcodes': 'barcodes.tsv',
        'matrix': 'matrix.mtx'
    }

    file_columns = {
        'matrix': None,
        'barcodes': ['cellkey', library_method_column],
        'genes': ['featurekey', 'featurename']
    }

    def __init__(self, mtx_info: MatrixInfo):
        self.info = mtx_info

    def unzip(self) -> None:
        """
        Extract files from top-level zip archive and remove archive.
        """
        log.info(f'Unzipping {self.info.zip_path}')
        os.mkdir(self.info.extract_path)
        with ZipFile(self.info.zip_path) as zipfile:
            for member in zipfile.namelist():
                filename = os.path.basename(member)
                source = zipfile.open(member)
                target = open(os.path.join(self.info.extract_path, filename), 'wb')
                with source, target:
                    shutil.copyfileobj(source, target)
        os.remove(self.info.zip_path)

        extracted_files = os.listdir(self.info.extract_path)
        assert 3 <= len(extracted_files) <= 4  # optional readme
        assert all(filename in extracted_files for filename in self.unproc_files.values())

    def preprocess(self):
        """
        Extract gzip files and transform for ScanPy compatibility.
        """
        # def process_mtx_files(self):

        with DirectoryChange(self.info.extract_path):
            for key, gzfilename in self.unproc_files.items():
                log.info(f'Processing file {gzfilename}')
                with gzip.open(gzfilename, 'rb') as gzfile:
                    self._preprocess(gzfile,
                                     keep_cols=self.file_columns[key],
                                     rename=self.proc_files[key])
                os.remove(gzfilename)

    def separate(self, strip_version: Optional[bool] = True) -> List[MatrixInfo]:
        """
        Split matrix into independent entities based on library construction method.'

        If the LCM is homogeneous, no change is made to the directory structure.
        Otherwise, a new directory is created within the extraction directoy for
        every observed library prep method and populated with the subset of
        barcodes.tsv corresponding to that method.
        Links are created for genes.tsv and matrix.mtx, which remain in the
        top-level extraction dir.
        :param strip_version: if true, best effort is made to strip version
        identifiers from the end of the LCM strings to conglomerate different
        versions of the same protocol.
        :return: list of paths to matrix directories.
        """
        barcodes_file = self.proc_files['barcodes']
        # Column headers have been removed so it takes a bit of footwork to
        # access the lib prep method field in pandas
        library_method_column_index = 1
        with DirectoryChange(self.info.extract_path):
            barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
            lib_con_data = barcodes.iloc[:, library_method_column_index]
            barcodes.drop(columns=barcodes.columns[library_method_column_index], inplace=True)
            if strip_version:
                lib_con_data = lib_con_data.apply(self.strip_version_suffix)
            distinct_method_count = lib_con_data.nunique()
            assert distinct_method_count > 0
            if distinct_method_count == 1:
                self.info.lib_con_method = first(lib_con_data)
                result = [self.info]
                self._write_tsv(barcodes_file, barcodes)
            else:
                result = []
                for lib_con_method, group in barcodes.groupby(lib_con_data.iloc.__getitem__):
                    os.mkdir(lib_con_method)
                    result.append(
                        MatrixInfo(source=self.info.source,
                                   project_uuid=self.info.project_uuid,
                                   zip_path=None,
                                   extract_path=os.path.join(self.info.extract_path, lib_con_method),
                                   lib_con_method=lib_con_method))

                    with DirectoryChange(lib_con_method):
                        self._write_tsv(barcodes_file, group)

                        genes_file = self.proc_files['genes']
                        os.symlink(f'../{genes_file}', genes_file)

                        # edit matrix file to match barcode subset
                        dst_matrix_file = self.proc_files['matrix']
                        src_matrix_file = f'../{dst_matrix_file}'
                        # scanpy needs MatrixMarket header even though pandas can't read it
                        with open(src_matrix_file) as f:
                            header = f.readline()
                        mtx = pd.read_csv(src_matrix_file, header=None, comment='%', sep=' ')
                        barcodes_column = mtx.iloc[1:, 1]
                        # align 0-based index to 1-based barcode line numbers
                        group.index += 1
                        keeprows = barcodes_column.map(group.index.__contains__)
                        # drop barcodes from other groups. +1 because of size header
                        mtx.drop(labels=1+np.flatnonzero(~keeprows), inplace=True)
                        barcodes_column = mtx.iloc[1:, 1]
                        # realign barcode references in matrix file to 1-based
                        # line numbers.
                        collapser = lru_cache(maxsize=1)(lambda i: one(np.flatnonzero(group.index == i))+1)
                        mtx.iloc[1:, 1] = barcodes_column.map(collapser)
                        # update size info
                        mtx.iloc[0, 1] = group.index.size
                        mtx.iloc[0, 2] = mtx.index.size - 1
                        with open(dst_matrix_file, 'w') as f:
                            f.write(header)
                        self._write_tsv(dst_matrix_file, mtx, mode='a')
        return result

    @classmethod
    def strip_version_suffix(cls, s: str):
        """
        >>> MatrixPreparer.strip_version_suffix('smartseq2_v2.4.0')
        'smartseq2'
        """
        version_suffix = r'_v\d+(?:\.\d+)*$'
        return re.sub(version_suffix, '', s)

    @classmethod
    def _preprocess(cls,
                    gzfile: gzip.GzipFile,
                    keep_cols: Optional[List[str]] = None,
                    rename: Optional[str] = None) -> None:
        """
        :param gzfile: gzipped file
        :param keep_cols: columns to keep. None to skip column selection.
        :param rename: name of the unzipped, processed file to write to
        """
        outfile_name = gzfile.name.rstrip('.gz') if rename is None else rename
        if keep_cols is None:
            with open(outfile_name, 'wb') as outfile:
                shutil.copyfileobj(gzfile, outfile)
        else:
            df = pd.read_table(gzfile, sep='\t')
            assert all(col in df.columns for col in keep_cols)
            df = df[keep_cols]
            cls._write_tsv(outfile_name, df)

    @classmethod
    def _write_tsv(cls, path: str, df: pd.DataFrame, **kwargs):
        df.to_csv(path, index=False, header=False, sep='\t', **kwargs)
