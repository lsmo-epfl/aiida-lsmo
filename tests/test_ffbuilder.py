#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for ff_builder"""  # pylint: disable=invalid-name
from aiida.orm import Dict, SinglefileData
from aiida.plugins import CalculationFactory
from aiida import engine

# Calculation objects
FFBuilder = CalculationFactory('lsmo.ff_builder')


def test_ff_builder():
    """Test force-field builder"""
    ff_parameters = Dict(
        {
            'ff_framework': 'UFF',
            'ff_molecules': {
                'CO2': 'TraPPE',
                'N2': 'TraPPE',
            },
            'shifted': False,
            'tail_corrections': True,
            'mixing_rule': 'Lorentz-Berthelot',
            'separate_interactions': True
        })

    results = engine.run(FFBuilder, ff_parameters)

    for pfile in results.values():
        assert isinstance(pfile, SinglefileData)

    assert results['ff_def'].get_attribute('filename') == 'force_field.def'

    with results['ff_def'].open('force_field.def') as handle:
        ff_def = handle.read()

    assert 'n2' in ff_def
    assert 'co2' in ff_def
    assert 'lennard-jones' in ff_def
