#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for ff_builder"""
from __future__ import absolute_import
from __future__ import print_function
from aiida.orm import Dict
from aiida.plugins import CalculationFactory
from aiida.engine import run_get_node

# Calculation objects
FFBuilder = CalculationFactory("lsmo.ff_builder")  # pylint: disable=invalid-name

ff_parameters = Dict( # pylint: disable=invalid-name
    dict={
        'ff_framework': 'Zeo-test',
        'ff_molecules': {
            'Na': 'Zeo-test',
            'CO2': 'TraPPE',
        },
        'shifted': False,
        'tail_corrections': True,
        'mixing_rule': 'Lorentz-Berthelot',
        'separate_interactions': False
    })

results, node = run_get_node(FFBuilder, ff_parameters)  # pylint: disable=invalid-name
print("Terminated ff_builder calcfunction, pk:", node.pk)
for key, val in results.items():
    #filepath = os.path.join(val._repository._get_base_folder().abspath, val.filename)
    print("Output:", val.pk, key)
