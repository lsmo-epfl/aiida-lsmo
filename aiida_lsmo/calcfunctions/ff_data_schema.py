# -*- coding: utf-8 -*-
"""Voluptuous schema for ff_data.yml"""
from voluptuous import Schema, Optional, Any

__all__ = ('FF_DATA_SCHEMA',)

NUMBER = Any(int, float)
FRAMEWORK_FF_SCHEMA = Schema(
    {
        'description': str,
        'atom_types': dict,
    },
    required=True,
)

MOLECULE_FF_SCHEMA = Schema(
    {
        'description': str,
        'atom_types': {
            str: {
                'force_field':
                    Any(['lennard-jones', NUMBER, NUMBER], ['feynman-hibbs-lennard-jones', NUMBER, NUMBER, NUMBER],
                        ['none'], 'dummy_separate'),
                'pseudo_atom':
                    list,
                Optional('force_field_mix'):
                    list,
            }
        },
        'atomic_positions': Any([[str, NUMBER, NUMBER, NUMBER]], 'flexible')
    },
    required=True)

MOLECULE_SCHEMA = Schema({
    'critical_constants': {
        'tc': NUMBER,
        'pc': NUMBER,
        'af': NUMBER,
    },
    str: MOLECULE_FF_SCHEMA,
},
                         required=True)

FF_DATA_SCHEMA = Schema(
    {
        'framework': {
            str: FRAMEWORK_FF_SCHEMA,
        },
        str: MOLECULE_SCHEMA,
    },
    required=True,
)
