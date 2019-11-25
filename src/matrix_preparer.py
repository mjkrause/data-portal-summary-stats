import gzip
import logging
import os

import pandas as pd
import shutil
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
        assert all(filename in extracted_files for filename in self.unproc_files)

    def preprocess(self):
        """
        Extract gzip files and transform for ScanPy compatibility.
        """
        # def process_mtx_files(self):

        file_columns = {
            self.unproc_files['matrix']: None,
            self.unproc_files['barcodes']: ['cellkey', self.library_method_column],
            self.unproc_files['genes']: ['featurekey, featurename']
        }

        with DirectoryChange(self.info.extract_path):
            for key, gzfilename in self.unproc_files.items():
                log.info(f'Processing file {gzfilename}')
                with gzip.open(gzfilename, 'rb') as gzfile:
                    self._preprocess(gzfile,
                                     keep_cols=file_columns[key],
                                     rename=self.proc_files[key])
                os.remove(gzfilename)

    def separate(self) -> List[MatrixInfo]:
        """
        Split matrix into independent entities based on library construction method.'

        If the LCM is homogeneous, no change is made to the directory structure.
        Otherwise, a new directory is created within the extraction directoy for
        every observed library prep method and populated with the subset of
        barcodes.tsv corresponding to that method.
        Links are created for genes.tsv and matrix.mtx, which remain in the
        top-level extraction dir.
        :return: list of paths to matrix directories.
        """
        barcodes_file = self.proc_files['barcodes']
        with DirectoryChange(self.info.extract_path):
            df = pd.read_csv(barcodes_file, sep='\t')
            distinct_method_count = df[self.library_method_column].nunique()
            assert distinct_method_count > 0
            if distinct_method_count == 1:
                self.info.lib_prep_method = df[self.library_method_column].iloc[0]
                result = [self.info]
            else:
                result = []
                for method, group in df.groupby(self.library_method_column):
                    os.mkdir(method)
                    with DirectoryChange(method):
                        group.to_csv(barcodes_file, index=False, header=False, sep='\t')
                        for filename in (self.proc_files['genes'],
                                         self.proc_files['matrix']):
                            os.symlink(f'../{filename}', filename)

                    result.append(
                        MatrixInfo(source=self.info.source,
                                   project_uuid=self.info.project_uuid,
                                   zip_path=None,
                                   extract_path=os.path.join(self.info.extract_path, method),
                                   lib_prep_method=method))

        return result

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
            df.to_csv(outfile_name, index=False, header=False, sep='\t')
