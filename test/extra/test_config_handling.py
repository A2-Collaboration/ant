#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_config_handling
----------------------------------

Tests for `parse_config` module.
"""

import unittest
from os.path import abspath, isfile, join as pjoin
import sys, shutil, tempfile
sys.path.append(abspath('extra/SimulationBlaster'))

from parse_config import Settings, read_config


class TestConfigHandling(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_dir = tempfile.mkdtemp()
        cls.settings = pjoin(cls.test_dir, 'sim_settings')
        cls.settings2 = pjoin(cls.test_dir, 'test_settings')

    def test_export_config(self):
        settings = Settings()
        self.assertTrue(settings.export(self.settings))
        self.assertTrue(isfile(self.settings))
        # test returning error on accidental overwrite
        self.assertFalse(settings.export(self.settings))

    def test_read_config(self):
        config = None
        with open(self.settings, 'r') as conf:
            config, _ = read_config(conf)
        self.assertIsNotNone(config)

    def test_changing_config_values(self):
        settings = Settings()
        test_val = 'OUTPUT_PATH'
        path = settings.get(test_val)
        new_path = pjoin(path, 'subdir')
        settings.set(test_val, new_path)
        self.assertEqual(settings.get(test_val), new_path)
        # test exporting value and reading in the change value, too
        self.assertTrue(settings.export(self.settings2))
        config = None
        with open(self.settings2, 'r') as conf:
            config, _ = read_config(conf)
        self.assertIsNotNone(config)
        self.assertEqual(config.get(test_val), new_path)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.test_dir)


if __name__ == '__main__':
    unittest.main()
