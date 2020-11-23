#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for molecule dictionaries"""  # pylint: disable=invalid-name

from aiida import orm
from aiida_lsmo.workchains.isotherm import get_molecule_dict


def test_molecule_dict():
    """Check that the molecule dict is valid"""
    molecule_dict = get_molecule_dict(orm.Str('xenon')).get_dict()

    assert molecule_dict['singlebead']
