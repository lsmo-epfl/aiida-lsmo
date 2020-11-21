# -*- coding: utf-8 -*-
"""Protocols for multi-stage CP2K calculations"""
from pathlib import Path
from voluptuous import Schema, Optional, Any
import ruamel.yaml as yaml  # does not convert OFF to False

__all__ = ('PROTOCOL_DIR', 'ISOTHERM_PROTOCOL_SCHEMA')

PROTOCOL_DIR = Path(__file__).resolve().parent

NUMBER = Any(int, float)
ELEMENT = Any('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
              'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
              'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
              'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
              'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn')
SETTINGS_SCHEMA = Schema({
    str: dict,
})
ISOTHERM_PROTOCOL_SCHEMA = Schema(
    {
        'protocol_description': str,
        'initial_magnetization': {
            ELEMENT: int,
        },
        'basis_set': {
            ELEMENT: str,
        },
        'pseudopotential': {
            ELEMENT: str,
        },
        'bandgap_thr_ev': NUMBER,
        'settings_0': SETTINGS_SCHEMA,
        Optional('settings_1'): SETTINGS_SCHEMA,
        'stage_0': SETTINGS_SCHEMA,
        Optional('stage_1'): SETTINGS_SCHEMA,
        Optional('stage_2'): SETTINGS_SCHEMA,
    },
    required=True)


def load_isotherm_protocol(*, path=None, tag=None):
    """Load isotherm protocol from yaml file (with validation)."""

    if path is not None:
        yaml_file = path
    elif tag is not None:
        yaml_file = PROTOCOL_DIR / (tag + '.yaml')
    else:
        raise ValueError('Provide either path or tag.')

    with open(yaml_file, 'r') as stream:
        protocol_dict = yaml.safe_load(stream)
        ISOTHERM_PROTOCOL_SCHEMA(protocol_dict)

    return protocol_dict
