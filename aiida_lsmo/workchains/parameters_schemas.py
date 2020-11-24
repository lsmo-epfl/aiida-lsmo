# -*- coding: utf-8 -*-
"""Schemas for validating input parameters of workchains.

Defines a couple of building blocks that are reused by many workchains.
"""
__all__ = ('Required', 'Optional', 'NUMBER', 'FF_PARAMETERS_VALIDATOR')
import voluptuous


def show_description(cls):
    """Adds description to representation of voluptuous Marker.

    This makes the description show up in the sphinx autodoc.
    """

    def __repr__(self):
        if self.description:
            return f'{self.__class__.__name__}({repr(self.schema)}, description={repr(self.description)})'
        return repr(self.schema)

    setattr(cls, '__repr__', __repr__)
    return cls


Required = show_description(voluptuous.Required)
Optional = show_description(voluptuous.Optional)
Any = voluptuous.Any
Schema = voluptuous.Schema

NUMBER = voluptuous.Any(int, float)

FF_PARAMETERS_VALIDATOR = voluptuous.Schema({
    Required('ff_framework',
             default='UFF',
             description='Forcefield of the structure  (used also as a definition of ff.rad for zeopp)'):
        str,
    Required('ff_separate_interactions',
             default=False,
             description='if true use only ff_framework for framework-molecule interactions in the FFBuilder'):
        bool,
    Required('ff_mixing_rule', default='Lorentz-Berthelot', description='Mixing rule'):
        Any('Lorentz-Berthelot', 'Jorgensen'),
    Required('ff_tail_corrections', default=True, description='Apply tail corrections.'):
        bool,
    Required('ff_shifted', default=False, description='Shift or truncate the potential at cutoff.'):
        bool,
    Required('ff_cutoff', default=12.0, description='CutOff truncation for the VdW interactions (Angstrom).'):
        NUMBER,
})
