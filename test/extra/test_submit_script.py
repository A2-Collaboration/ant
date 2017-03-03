#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_submit_script
----------------------------------

Tests for some `Ant_simulation_blaster` routines.
"""

import unittest
from os.path import abspath, isfile, join as pjoin
import sys, shutil, tempfile
sys.path.append(abspath('extra/SimulationBlaster'))

from Ant_simulation_blaster import *


class TestConfigHandling(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_dir = tempfile.mkdtemp()
        cls.channels = [["p eta' [g g]", 10, 1000],
                ['"Cocktail"', 10, 1000],
                ['"Gun: g g p"', 10, 1000]]

    def test_decay_strings(self):
        pluto = get_decay_string(self.channels[0][0])
        cocktail = get_decay_string(self.channels[1][0])
        mcgun = get_decay_string(self.channels[2][0])
        self.assertEqual(pluto, 'etap_gg')
        self.assertEqual(cocktail, 'cocktail')
        self.assertEqual(mcgun, 'ggp-gun')

    def test_channels_sanity_check(self):
        self.assertTrue(sanity_check_channels('Ant-pluto', [self.channels[0]]))
        self.assertTrue(sanity_check_channels('Ant-cocktail', [self.channels[1]]))
        self.assertTrue(sanity_check_channels('Ant-mcgun', [self.channels[2]]))
        self.assertTrue(sanity_check_channels('fallback_own_generator', self.channels))

    def test_mc_generator_command_preparation(self):
        settings = Settings()
        settings.set('OUTPUT_PATH', self.test_dir)
        settings.set('A2_GEANT_PATH', self.test_dir)
        check_directories(settings, force=True)
        settings.set('Emin', 1000)
        settings.set('Emax', 1600)
        settings.set('COCKTAIL_SETUP', 'Setup_2014_07_EPT_Prod')
        settings.set('GUN_THETA', '20 160')
        settings.set('GUN_OPENING', 50)

        # test Ant-cocktail command creation
        cocktail = [self.channels[1][0], 'cocktail_out', self.channels[1][2]]
        cocktail_cmd = create_mcgen_cmd(settings, 'Ant-cocktail', *cocktail)
        cocktail_expected = 'Ant-cocktail -o cocktail_out -n %d -s Setup_2014_07_EPT_Prod' \
                % self.channels[1][2]
        cocktail_expected += ' ' + settings.get('AddFlags')
        self.assertEqual(cocktail_cmd, cocktail_expected)

        # test Ant-mcgun command creation
        mcgun = [self.channels[2][0], 'mcgun_out', self.channels[2][2]]
        mcgun_cmd = create_mcgen_cmd(settings, 'Ant-mcgun', *mcgun)
        mcgun_expected = 'Ant-mcgun -p g -p g -p p -o mcgun_out -n %d --EMin %f --EMax %f' \
                % (self.channels[2][2], 1000, 1600)
        mcgun_expected += ' --thetaMin {} --thetaMax {} --openingAngle {}' \
                .format(*settings.get('GUN_THETA').split(), settings.get('GUN_OPENING'))
        mcgun_expected += ' ' + settings.get('AddFlags')
        self.assertEqual(mcgun_cmd, mcgun_expected)

        # test Ant-pluto command creation
        mcgen_data = get_path(settings.get('OUTPUT_PATH'), settings.get('MCGEN_DATA'))
        pluto_decay = get_decay_string(self.channels[0][0])
        pluto_file = get_path(mcgen_data, get_file_name(MCGEN_PREFIX, pluto_decay, 1))
        expected_file = get_path(mcgen_data, MCGEN_PREFIX + '_etap_gg_0001.root')
        pluto = [self.channels[0][0], pluto_file, self.channels[0][2]]
        pluto_cmd = create_mcgen_cmd(settings, 'Ant-pluto', *pluto)
        pluto_expected = 'Ant-pluto --reaction %s -o %s -n %d --Emin %f --Emax %f --no-bulk' \
                % (self.channels[0][0], expected_file, self.channels[0][2], 1000, 1600)
        pluto_expected += ' ' + settings.get('AddFlags')
        self.assertEqual(pluto_cmd, pluto_expected)

        # test geant command creation
        settings.set('GeantFlags', "'s/regex/looking/g' options 's/another/regex/' --more flags")
        mcgen_file = '/some/mc/input'
        geant_file = '/geant/simulation/output'
        geant_cmd = create_geant_cmd(settings, '/path/to/runGeant.sh', mcgen_file, geant_file)
        geant_expected = '/path/to/runGeant.sh %s %s' % (mcgen_file, geant_file)
        geant_expected += " 's/regex/looking/g'  's/another/regex/' options --more flags"
        self.assertEqual(geant_cmd, geant_expected)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.test_dir)


if __name__ == '__main__':
    unittest.main()
