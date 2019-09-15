#!/usr/bin/env python3

import os
import unittest
from src.utils import convert_size, get_blacklist


class TestUtils(unittest.TestCase):

    def setUp(self) -> None:

        # Create data to test blacklist and write file to project root.
        self.blacklisted = ['979417ef-59d4-420d-b60f-a1c4c349d8a4',
                            '6719e02c-dae2-4e12-8978-3f28f01c9cb6',
                            '41e4d281-1505-4b31-8df4-1fed05e1ef79']
        file = 'blacklist'
        with open(file, 'w+') as fp:
            for uid in self.blacklisted:
                fp.writelines(uid + '\n')

    def tearDown(self) -> None:
        os.remove('blacklist')

    def test_convert_size(self):
        self.assertEqual(convert_size(0), '0 B')
        self.assertEqual(convert_size(1024), '1.0 KB')
        self.assertEqual(convert_size(2124569754), '1.98 GB')

    def test_get_blacklisted_files(self):
        do_not_process = get_blacklist()  # happy path
        self.assertEqual(self.blacklisted, do_not_process)


if __name__ == '__main__':
    unittest.main()
