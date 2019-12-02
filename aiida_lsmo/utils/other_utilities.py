# -*- coding: utf-8 -*-
"""Other utilities"""

from __future__ import absolute_import
from aiida.orm import Dict, CifData
from aiida.engine import calcfunction


def dict_merge(dct, merge_dct):
    """ Taken from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    import collections
    for k in merge_dct.keys():
        if (k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


@calcfunction
def aiida_dict_merge(to_dict, from_dict):
    """Merge two aiida Dict objects."""
    to_dict = to_dict.get_dict()

    if isinstance(from_dict, Dict):
        from_dict = from_dict.get_dict()

    dict_merge(to_dict, from_dict)

    return Dict(dict=to_dict)


@calcfunction
def aiida_cif_merge(aiida_cif_a, aiida_cif_b):
    """Merge the coordinates of two CifData into a sigle one. Note: the two unit cells must be the same."""
    import ase
    ase_a = aiida_cif_a.get_ase()
    ase_b = aiida_cif_b.get_ase()
    if not (ase_a.cell == ase_b.cell).all():
        raise ValueError('Attempting to merge two CifData with different unit cells.')
    ase_ab = ase.Atoms(  #Maybe there is a more direct way...
        symbols=list(ase_a.symbols) + list(ase_b.symbols),
        cell=ase_a.cell,
        positions=list(ase_a.positions) + list(ase_b.positions),
        pbc=True)
    return CifData(ase=ase_ab)
