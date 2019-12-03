import unittest

import boto3
from moto import (
    mock_sts,
    mock_s3,
)

from src import Config


@mock_s3
@mock_sts
class S3TestCase(unittest.TestCase):
    config = Config('dev')

    @property
    def bucket_name(self):
        return self.config.s3_bucket_name

    def setUp(self) -> None:
        self.client = boto3.client('s3')
        self.client.create_bucket(Bucket=self.bucket_name)

    def tearDown(self) -> None:
        response = self.client.list_objects_v2(Bucket=self.bucket_name)
        keys = [obj['Key'] for obj in response['Contents']]
        for key in keys:
            self.client.delete_object(Bucket=self.bucket_name,
                                      Key=key)
        self.client.delete_bucket(Bucket=self.bucket_name)
