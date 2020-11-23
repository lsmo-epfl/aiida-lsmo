# -*- coding: utf-8 -*-
"""Schemas for validating input parameters of workchains.

Defines a couple of building blocks that are reused by many workchains.
"""
from voluptuous import Schema, Required, Any

NUMBER = Any(int, float)

FF_PARAMETERS_VALIDATOR = Schema({
    Required('ff_framework', default='UFF'):
        str,  # Forcefield of the structure  (used also as a definition of ff.rad for zeopp)
    Required('ff_separate_interactions', default=False):
        bool,  # if true use only ff_framework for framework-molecule interactions in the FFBuilder
    Required('ff_mixing_rule', default='Lorentz-Berthelot'): Any('Lorentz-Berthelot', 'Jorgensen'),  # Mixing rule
    Required('ff_tail_corrections', default=True): bool,  # Apply tail corrections.
    Required('ff_shifted', default=False): bool,  # Shift or truncate the potential at cutoff.
    Required('ff_cutoff', default=12.0): NUMBER,  # CutOff truncation for the VdW interactions (Angstrom).
})
