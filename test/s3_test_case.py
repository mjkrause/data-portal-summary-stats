import unittest

import boto3
from moto import (
    mock_sts,
    mock_s3,
)

from dpss.config import config


@mock_s3
@mock_sts
class S3TestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.managed_buckets = (config.s3_matrix_bucket_name, config.s3_figure_bucket_name)
        self.client = boto3.client('s3')
        for bucket in self.managed_buckets:
            self.client.create_bucket(Bucket=bucket)

    def tearDown(self) -> None:
        for bucket in self.managed_buckets:
            response = self.client.list_objects_v2(Bucket=bucket)
            keys = [obj['Key'] for obj in response.get('Contents', [])]
            for key in keys:
                self.client.delete_object(Bucket=bucket,
                                          Key=key)
            self.client.delete_bucket(Bucket=bucket)
