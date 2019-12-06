import os
import math
from tempfile import TemporaryDirectory
from typing import List

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
    with open('blacklist', 'r') as fp:
        return [line.rstrip('\n') for line in fp]


def file_id(path):
    """
    Filename without preceding directories or extensions.
    """
    return first(os.path.basename(path).split('.'))


def remove_ext(path, ext):
    """
    Remove a file extension. No effect if provided extension is missing.
    """
    parts = path.rsplit(ext, 1)
    if len(parts) == 2 and parts[1] == '':
        return first(parts)
    else:
        return path


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
        return self.new_dir

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        os.chdir(self.old_dir)


class TemporaryDirectoryChange(DirectoryChange):
    """
    Directory change context manager that creates and cleans up a temporary
    directory.
    """

    def __init__(self):
        self.tmp = TemporaryDirectory()
        super().__init__(self.tmp.name)

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        super().__exit__(exc_type, exc_val, exc_tb)
        self.tmp.cleanup()
