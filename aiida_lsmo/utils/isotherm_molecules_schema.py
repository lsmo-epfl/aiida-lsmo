# -*- coding: utf-8 -*-
"""Voluptuous schema for isotherm_molecules.yaml"""
from voluptuous import Schema, Optional, Any

__all__ = ('ISOTHERM_MOLECULES_SCHEMA')

NUMBER = Any(int, float)
MOLECULE_SCHEMA = Schema(
    {
        'name': str,
        'forcefield': str,
        'molsatdens': NUMBER,
        'proberad': NUMBER,
        'singlebead': bool,
        'charged': bool,
        Optional('rosenbluth'): NUMBER,
        Optional('pressure_zero'): NUMBER,
    },
    required=True)

ISOTHERM_MOLECULES_SCHEMA = Schema({
    str: MOLECULE_SCHEMA,
})
