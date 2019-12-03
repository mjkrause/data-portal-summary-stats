import logging
import os
from typing import (
    List,
    Dict,
)

import boto3

from src.matrix_info import MatrixInfo

log = logging.getLogger(__name__)


class S3Service:

    def __init__(self, bucket_name: str, key_prefixes: Dict[str, str]):
        log.info('Initializing S3 service...')
        self.client = boto3.client('s3')
        self.bucket_name = bucket_name
        self.key_prefixes = key_prefixes

    def list_bucket(self, target: str) -> List[str]:
        response = self.client.list_objects_v2(Bucket=self.bucket_name,
                                               Prefix=self.key_prefixes[target])
        return [obj['Key'] for obj in response.get('Contents', [])]

    def download(self, target: str, filename: str) -> None:
        self.client.download_file(Bucket=self.bucket_name,
                                  Key=self.key_prefixes[target] + filename,
                                  Filename=filename)

    def get_blacklist(self) -> List[str]:
        response = self.client.get_object(Bucket=self.bucket_name, Key='blacklist')
        bytes_string = response['Body'].read()
        return bytes_string.decode().strip('\n').split('\n')

    def upload_figures(self, mtx_info: MatrixInfo) -> None:
        figures = os.listdir('figures/')
        for figure in figures:
            key = self.key_prefixes['figures'] + mtx_info.project_uuid + '/' + figure
            log.info(f'Uploading {figure} to S3 bucket {self.bucket_name} as {key}')
            self.client.upload_file(Filename=f'figures/{figure}',
                                    Bucket=self.bucket_name,
                                    Key=key)
