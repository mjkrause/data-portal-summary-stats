#!/usr/bin/env python3

import sys
import math
import logging

logger = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)  # console handler to output to stdout
ch.setLevel(logging.INFO)
logger.addHandler(ch)

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