#!/usr/bin/env python3
import os
import sys
import math
from typing import List

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


def get_blacklist() -> List[str]:
    try:
        with open('blacklist', 'r') as fp:
            return [line.rstrip('\n') for line in fp]
    except FileNotFoundError:
        sys.exit('\n     File "blacklist" not found. Please create and populate it.')


def remove_extension(filename: str, ext: str) -> str:
    """Removes one extension in a filename following character "."."""
    assert '.' in filename
    return first(filename.split(f'.{ext}'))


def file_id(path):
    return first(os.path.basename(path).split('.'))


class DirectoryChange:
    """
    Context manager facilitating an undoable temporary switch to another working
    directory.
    """
    def __init__(self, new_dir):
        self.new_dir = new_dir

    def __enter__(self):
        self.old_dir = os.getcwd()
        os.chdir(self.new_dir)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.old_dir)
