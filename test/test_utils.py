#!/usr/bin/env python3

import os
import unittest
from src.utils import convert_size, get_blacklist, remove_extension


class TestUtils(unittest.TestCase):

    def setUp(self) -> None:

        self.proj_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Create data to test blacklist and write file to project root.
        self.blacklisted = ['979417ef-59d4-420d-b60f-a1c4c349d8a4',
                            '6719e02c-dae2-4e12-8978-3f28f01c9cb6',
                            '41e4d281-1505-4b31-8df4-1fed05e1ef79']
        file = 'blacklist'
        with open(os.path.join(self.proj_root, file), 'w+') as fp:
            for uid in self.blacklisted:
                fp.writelines(uid + '\n')

    def tearDown(self) -> None:
        blacklist = os.path.join(self.proj_root, 'blacklist')
        if os.path.isfile(blacklist):
            os.remove(blacklist)

    def test_convert_size(self):
        self.assertEqual(convert_size(0), '0 B')
        self.assertEqual(convert_size(1024), '1.0 KB')
        self.assertEqual(convert_size(2124569754), '1.98 GB')

    def test_get_blacklisted_files(self):
        os.chdir(self.proj_root)
        do_not_process = get_blacklist()  # happy path
        self.assertEqual(self.blacklisted, do_not_process)

    def test_get_blacklisted_files_no_blacklist(self):
        os.chdir(self.proj_root)
        assert 'blacklist' in os.listdir()
        os.remove('blacklist')
        with self.assertRaises(SystemExit):
            get_blacklist()

    def test_remove_extension(self):
        filename = 'file.ext1.ext2'
        expected = 'file.ext1'
        observed = remove_extension(filename, 'ext2')
        self.assertEqual(expected, observed)

        with self.assertRaises(AssertionError):
            filename = 'file'
            remove_extension(filename, 'ext')

if __name__ == '__main__':
    unittest.main()
