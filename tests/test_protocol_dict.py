#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for protocol dictionaries"""  # pylint: disable=invalid-name
from aiida_lsmo.workchains.cp2k_multistage_protocols import PROTOCOL_DIR, load_isotherm_protocol


def test_protocol_dict():
    """Check that the protocol dicts are valid"""
    for protocol_yaml in PROTOCOL_DIR.glob('*.yaml'):
        tag = protocol_yaml.stem
        _ = load_isotherm_protocol(tag=tag)
