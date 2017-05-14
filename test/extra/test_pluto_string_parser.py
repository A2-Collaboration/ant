#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pluto_string_parser
----------------------------------

Tests for `parse_pluto_string` module.
"""

import unittest
from os.path import abspath
import sys
sys.path.append(abspath('extra/SimulationBlaster'))

import parse_pluto_string


class TestPlutoStringParser(unittest.TestCase):

    def test_pluto_string(self):
        channel = "eta' [pi0 [g g] pi0 [dilepton [e+ e-] g] eta [pi0 [g g] pi0 [dilepton [e+ e-] g] pi0 [g g]]]"
        initial_state = parse_pluto_string.get_initial_state(channel)
        self.assertEqual(len(initial_state), 1)
        self.assertEqual(*initial_state, "eta'")
        decay_string = parse_pluto_string.get_decay_string(channel, 1)
        self.assertEqual(decay_string, 'etap_etapi0pi0_e+e+e-e-8g')
        decay_string_more = parse_pluto_string.get_decay_string(channel, 2)
        self.assertEqual(decay_string_more, 'etap_etapi0pi0_dilepton3g3pi0_e+e+e-e-8g')
        channel_regex_needed = "eta' [g rho0 [g pi0 [g g]]]"
        decay_string_regex = parse_pluto_string.get_decay_string(channel_regex_needed, 1)
        self.assertEqual(decay_string_regex, 'etap_grho0_4g')


if __name__ == '__main__':
    unittest.main()
