"""Other utilities"""

from aiida.orm import Dict, CifData, StructureData
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
    cif_ab = CifData(ase=ase_ab, filename='fragments_a_b.cif')  #TODO: check why the filename is not assigned. # pylint: disable=fixme
    cif_ab.label = 'Loaded structure'
    cif_ab.description = 'Fragment A: {} atoms, fragment B: {} atoms.'.format(len(ase_a), len(ase_b))
    return cif_ab


@calcfunction
def aiida_structure_merge(aiida_structure_a, aiida_structure_b):
    """Merge the coordinates of two StructureData into a sigle one. Note: the two unit cells must be the same."""
    import ase
    ase_a = aiida_structure_a.get_ase()
    ase_b = aiida_structure_b.get_ase()
    if not (ase_a.cell == ase_b.cell).all():
        raise ValueError('Attempting to merge two StructureData with different unit cells.')
    ase_ab = ase.Atoms(  #Maybe there is a more direct way...
        symbols=list(ase_a.symbols) + list(ase_b.symbols),
        cell=ase_a.cell,
        positions=list(ase_a.positions) + list(ase_b.positions),
        pbc=True)
    return StructureData(ase=ase_ab)


@calcfunction
def get_structure_from_cif(cifdata):
    """Convert StructureData to CifData maintaining the provenance."""
    return cifdata.get_structure()


@calcfunction
def get_cif_from_structure(structuredata):
    """Convert CifData to StructureData maintaining the provenance."""
    return structuredata.get_cif()
