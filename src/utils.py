#!/usr/bin/env python3

import sys
import math
import boto3
import logging
from more_itertools import first

logger = logging.getLogger(__name__)


# The following is adapted from
# https://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python/14822210#14822210
def convert_size(size_bytes: float) -> str:
    if size_bytes == 0:
        return "0 B"
    order_of_magnitude = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)

    return f'{s} {order_of_magnitude[i]}'


def get_blacklist() -> list:
    do_not_process = []
    try:
        with open('blacklist', 'r') as fp:
            for line in fp:
                do_not_process.append(line.strip('\n'))
    except FileNotFoundError:
        sys.exit('\n     File "blacklist" not found. Please create and populate it.')

    return do_not_process


def get_blacklist_from_s3(client: boto3.client, bucket: str, key: str) -> list:
    response = client.get_object(Bucket=bucket, Key=key)
    bytes_string = response['Body'].read()

    return bytes_string.decode().strip('\n').split('\n')


def remove_extension(filename: str, ext: str) -> str:
    """Removes one extension in a filename following character "."."""
    assert '.' in filename

    return first(filename.split(f'.{ext}'))

