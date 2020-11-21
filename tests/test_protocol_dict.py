#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for protocol dictionaries"""  # pylint: disable=invalid-name
import glob
from aiida_lsmo.workchains.cp2k_multistage_protocols import PROTOCOL_DIR, load_isotherm_protocol


def test_protocol_dict():
    """Check that the protocol dicts are valid"""
    for protocol_yaml in glob.glob(str(PROTOCOL_DIR / '*.yaml')):
        _ = load_isotherm_protocol(path=protocol_yaml)
