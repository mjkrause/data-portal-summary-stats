from functools import lru_cache
import gzip
import logging
import os
import warnings

from more_itertools import (
    first,
    one,
)
import numpy as np
import pandas as pd
import shutil
from typing import (
    List,
    Optional,
    Tuple,
    Callable,
    Union,
)
from zipfile import ZipFile

from src.matrix_info import MatrixInfo
from src.utils import (
    DirectoryChange,
    remove_ext,
)

log = logging.getLogger(__name__)


class MatrixPreparer:
    library_method_column = 'analysis_protocol.protocol_core.protocol_id'

    zipped_files = {
        'genes': 'genes.tsv.gz',
        'barcodes': 'cells.tsv.gz',
        'matrix': 'matrix.mtx.gz'
    }

    unproc_files = {
        'genes': 'genes.tsv',
        'barcodes': 'cells.tsv',
        'matrix': 'matrix.mtx'
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
        Extract files from top-level zip archive, uncompress .gz files, and
        remove archive.
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
        assert all(filename in extracted_files for filename in self.zipped_files.values())

        with DirectoryChange(self.info.extract_path):
            for gzfilename in self.zipped_files.values():
                with gzip.open(gzfilename, 'rb') as gzfile:
                    filename = remove_ext(gzfilename, '.gz')
                    with open(filename, 'wb') as outfile:
                        shutil.copyfileobj(gzfile, outfile)
                        os.remove(gzfilename)

    def preprocess(self):
        """
        Extract gzip files and transform for ScanPy compatibility.
        """
        with DirectoryChange(self.info.extract_path):
            for filekey, filename in self.unproc_files.items():
                log.info(f'Processing file {filename}')
                self._preprocess(remove_ext(filename, '.gz'),
                                 keep_cols=self.file_columns[filekey],
                                 rename=self.proc_files[filekey])

    def prune(self, keep_frac: float) -> None:
        """
        Shrink matrix by removing entries at random to reach a desired fraction
        of the original size.
        :param keep_frac: the fraction of entries to keep, e.g. 0.05 to remove
        95% of the entries.
        :return: nothing, operations occurs on disk.
        """
        if not (0 < keep_frac <= 1):
            raise ValueError(f'Invalid prune fraction: {keep_frac}')
        mtx_path = os.path.join(self.info.extract_path, self.proc_files['matrix'])
        header, mtx = self._read_matrixmarket(mtx_path)
        entries = mtx.index[1:]
        keep_rows = np.random.choice(entries, round(keep_frac * entries.size), replace=False)
        mtx = self._filter_matrix_entries(mtx, pd.Series(keep_rows))
        self._write_matrixmarket(mtx_path, header, mtx)

    def separate(self, strip_version: bool = True) -> List[MatrixInfo]:
        """
        Split matrix into independent entities based on library construction method.

        If the LCM is homogeneous, no change is made to the directory structure.
        Otherwise, a new directory is created within the extraction directory
        for every observed library prep method and populated with the subset of
        matrix.mtx corresponding to that method.
        Links are created for the row and column tsv files, which remain in the
        top-level extraction dir.
        :param strip_version: if true, best effort is made to strip version
        identifiers from the end of the LCM strings to conglomerate different
        versions of the same protocol.
        :return: list of MatrixInfo objects describing the results of the
        separation.
        """
        # Column headers have been removed so it takes a bit of footwork to
        # access the lib prep method field in pandas
        library_method_column_index = 1
        barcode_column_index = 1

        log.info('Separating by library construction method...')

        with DirectoryChange(self.info.extract_path):
            barcodes_file = self.proc_files['barcodes']
            barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
            lib_con_data = barcodes.pop(barcodes.columns[library_method_column_index])
            self._write_tsv(barcodes_file, barcodes)

            if strip_version:
                lib_con_data = lib_con_data.map(self.strip_version_suffix)

            lib_con_methods = lib_con_data.unique()
            assert len(lib_con_methods) > 0

            if len(lib_con_methods) == 1:
                log.info('Homogeneous LCM.')
                self.info.lib_con_method = first(lib_con_data)
                return [self.info]
            else:
                for lib_con_method in lib_con_methods:
                    log.info(f'Consolidating {lib_con_method} cells')
                    os.mkdir(lib_con_method)
                    with DirectoryChange(lib_con_method):
                        for filekey, filename in self.proc_files.items():
                            if filekey == 'matrix':
                                def matrix_filter(entry):
                                    barcode_lineno = entry[barcode_column_index]
                                    entry_method = lib_con_data.iloc[barcode_lineno - 1]
                                    return entry_method == lib_con_method

                                header, mtx = self._read_matrixmarket(f'../{filename}')
                                mtx = self._filter_matrix_entries(mtx, matrix_filter)
                                self._write_matrixmarket(filename, header, mtx)
                            else:
                                os.symlink(f'../{filename}', filename)
                return [MatrixInfo(source=self.info.source,
                                   project_uuid=self.info.project_uuid,
                                   zip_path=None,
                                   extract_path=os.path.join(self.info.extract_path, lib_con_method),
                                   lib_con_method=lib_con_method)
                        for lib_con_method in lib_con_methods]

    def rezip(self, zip_path: str = None, remove_dir: bool = False) -> None:
        """
        Compress unprocessed files to zip archive.
        This is the inverse operation of `unzip`.
        :param remove_dir: whether the extraction directory should be removed
        after zipping.
        :param zip_path: path of resulting zip archive. Updates matrix info if
        provided.
        """
        if zip_path is not None:
            self.info.zip_path = zip_path

        with ZipFile(self.info.zip_path, 'w') as zipfile:
            with DirectoryChange(self.info.extract_path):
                for gzfilename, filename in zip(self.zipped_files.values(),
                                                self.unproc_files.values()):
                    with open(filename, 'rb') as infile:
                        with gzip.open(gzfilename, 'wb') as gzfile:
                            shutil.copyfileobj(infile, gzfile)
                            os.remove(filename)
                    zipfile.write(gzfilename)

        if remove_dir:
            shutil.rmtree(self.info.extract_path)

    @classmethod
    def strip_version_suffix(cls, protocol: str):
        """
        >>> MatrixPreparer.strip_version_suffix('smartseq2_v2.4.0')
        'smartseq2'
        """
        parts = protocol.split('_')
        if len(parts) != 2:
            raise RuntimeError(f'Library protocol {protocol} does not properly define a version')
        return parts[0]

    @classmethod
    def _filter_matrix_axis(cls,
                            mtx: pd.DataFrame,
                            axis_idx_subset: pd.Index,
                            axis: int):
        """
        Remove entries in a matrix market matrix file to reflect a reduced set
        of values is a matrix market row or column file.
        This is not necessary to remove entries from the matrix. It is only
        necessary when entries have been removed from the row or column file and
        the matrix file must be updated to avoid incorrect line number
        references.
        :param mtx: matrix market entry dataframe.
        :param axis_idx_subset: index containing a subset of the line numbers
        referenced in either the rows or columns of the matrix file.
        :param axis: 1 for rows, 2 for columns.
        :return: matrix with entries referencing absent row or column numbers
        removed and size info updated.
        """
        # align index to 1-based line numbers in axis (row or column) file
        axis_idx_subset = axis_idx_subset + 1

        # remove entries referencing missing row/columns
        mtx = cls._filter_matrix_entries(mtx, lambda entry: entry[axis] in axis_idx_subset)

        # realign references in matrix file to 1-based line numbers in axis file.
        collapser = lru_cache(maxsize=1)(lambda i: one(np.flatnonzero(axis_idx_subset == i)) + 1)
        mtx.iloc[1:, 1] = mtx.iloc[1:, 1].map(collapser)

        # update size of axis file
        mtx.iloc[0, axis] = axis_idx_subset.size

        return mtx

    @classmethod
    def _filter_matrix_entries(cls,
                               mtx: pd.DataFrame,
                               which_rows: Union[pd.Series, Callable[[pd.Series], bool]]) -> pd.DataFrame:
        """
        Remove rows from a matrix-market format matrix file, leaving the row and
        column files unchanged.
        :param mtx: dataframe of matrix entries in matrix market format.
        :param which_rows: predicate function or Series of integer indices.
        :return: the modified matrix.
        """
        if callable(which_rows):
            keep_labels = mtx.iloc[1:, :].apply(which_rows, axis=1)
            mtx = mtx.drop(labels=1 + np.flatnonzero(~keep_labels))
        else:
            mtx = mtx.iloc[pd.concat([pd.Series(0), which_rows]), :]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Pandas complains about this modification because it at this point
            # `mtx` is a slice of of the original parameter and the following
            # assignment will not propagate past the slice. Which is exactly
            # what is intended here.
            mtx.iloc[0, 2] = mtx.index.size - 1
        return mtx

    @classmethod
    def _preprocess(cls,
                    infile_name: str,
                    keep_cols: Optional[List[str]] = None,
                    rename: Optional[str] = None) -> None:
        """
        :param infile_name: name of un-gzipped, unprocessed file
        :param keep_cols: columns to keep. None to skip column selection.
        :param rename: name of the unzipped, processed file to write to
        """
        if infile_name == cls.unproc_files['matrix']:
            # ScanPy requires gene expression to be int for some reason
            # TODO Normalize?
            header, df = cls._read_matrixmarket(infile_name)
            df.iloc[:, 2] = np.rint(df.iloc[:, 2]).astype(int)
            cls._write_matrixmarket(infile_name, header, df)
        elif keep_cols is not None:
            df = pd.read_table(infile_name, sep='\t')
            assert all(col in df.columns for col in keep_cols)
            df = df[keep_cols]
            cls._write_tsv(infile_name, df)

        if rename is not None:
            os.rename(infile_name, rename)

    @classmethod
    def _write_tsv(cls, path: str, df: pd.DataFrame, **pd_kwargs) -> None:
        df.to_csv(path, index=False, header=False, sep='\t', **pd_kwargs)

    @classmethod
    def _read_matrixmarket(cls, path: str, **pd_kwargs) -> Tuple[str, pd.DataFrame]:
        with open(path) as f:
            header = f.readline()
        assert header.startswith('%')
        df = pd.read_csv(path, header=None, comment='%', sep=' ', **pd_kwargs)
        return header, df

    @classmethod
    def _write_matrixmarket(cls, path: str, header: str, df: pd.DataFrame, **pd_kwargs) -> None:
        # ScanPy needs MatrixMarket header even though pandas can't read it
        with open(path, 'w') as f:
            f.write(header)
        df.to_csv(path, index=False, header=False, sep=' ', mode='a', **pd_kwargs)
