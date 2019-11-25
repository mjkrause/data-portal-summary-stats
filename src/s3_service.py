import logging
import os
from typing import List

import boto3

from src.matrix_info import MatrixInfo

log = logging.getLogger(__name__)


class S3Service:

    def __init__(self, bucket_name: str, key_prefix: str):
        self.client = boto3.client('s3')
        self.bucket_name = bucket_name
        self.key_prefix = key_prefix

    def list_bucket(self):
        return self.client.list_objects_v2(Bucket=self.bucket_name,
                                           Prefix=self.key_prefix)

    def download(self, filename):
        self.client.download_file(Bucket=self.bucket_name,
                                  Key=filename,
                                  Filename=filename)

    def get_blacklist(self) -> List[str]:
        response = self.client.get_object(Bucket=self.bucket_name, Key='blacklist')
        bytes_string = response['Body'].read()
        return bytes_string.decode().strip('\n').split('\n')

    def upload_figures(self, mtx_info: MatrixInfo) -> None:
        figures = os.listdir('figures/')
        for figure in figures:
            key = self.key_prefix + mtx_info.project_uuid + '/' + figure
            log.info(f'Uploading {figure} to S3 bucket {self.bucket_name} as {key}')
            self.client.upload_file(Filename=f'figures/{figure}',
                                    Bucket=self.bucket_name,
                                    Key=key)
