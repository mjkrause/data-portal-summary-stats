#!/usr/bin/env python3

import os
from tempfile import TemporaryDirectory
import unittest
from dpss.utils import (
    convert_size,
    file_id,
    DirectoryChange,
    remove_ext,
)


class TestUtils(unittest.TestCase):

    def test_convert_size(self):
        self.assertEqual(convert_size(0), '0 B')
        self.assertEqual(convert_size(1024), '1.0 KB')
        self.assertEqual(convert_size(2124569754), '1.98 GB')

    def test_file_id(self):
        self.assertEqual(file_id('/fish/dish'), 'dish')
        self.assertEqual(file_id('dish.exe'), 'dish')
        self.assertEqual(file_id('foo/bar/baz.egg'), 'baz')
        self.assertEqual(file_id('dotted.dir/file'), 'file')
        self.assertEqual(file_id('eggs.ham.spam'), 'eggs')

    def test_remove_ext(self):
        self.assertEqual(remove_ext('x.zip', '.zip'), 'x')
        self.assertEqual(remove_ext('x.zip', '.whoops'), 'x.zip')
        self.assertEqual(remove_ext('x/y/z.ship.zip', '.zip'), 'x/y/z.ship')
        self.assertEqual(remove_ext('x/y/z.ship.zip', '.ship'), 'x/y/z.ship.zip')

    def test_directory_change(self):
        old_dir = os.getcwd()
        with TemporaryDirectory() as test_dir_1, TemporaryDirectory() as test_dir_2:
            with DirectoryChange(test_dir_1) as new_dir:
                self.assertEqual(new_dir, test_dir_1)
                os.chdir(test_dir_2)  # go somewhere else
        self.assertEqual(os.getcwd(), old_dir)


if __name__ == '__main__':
    unittest.main()
