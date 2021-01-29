# -*- coding: utf-8 -*-
"""Utilities related to CP2K."""
from aiida_lsmo.utils import HARTREE2EV


def get_kinds_info(atoms):
    """Get kinds information from ASE atoms

    :param atoms: ASE atoms instance
    :returns: list of kind_info dictionaries (keys: 'kind', 'element', 'magnetization')
    """
    symbols = sorted(set(atoms.get_chemical_symbols()))

    kinds_info = []
    for symbol in symbols:
        ats = [atom for atom in atoms if atom.symbol == symbol]

        tags = {}
        for atom in ats:
            # we assume that atoms are already tagged properly (atoms with the same tag have the same properties)
            tags[atom.tag] = {'element': atom.symbol, 'magnetization': {atom.magmom}, 'tag': {atom.tag}}

        if len(tags) == 1:
            kind = tags[0]
            kind['kind'] = kind['element']
            kinds_info.append(kind)
        else:
            for tag, kind in tags.items():
                kind['kind'] = f"{kind['element']}{str(tag)}"
                kinds_info.append(kind)

    return kinds_info


def get_multiplicity_section(atoms, protocol):
    """ Compute the total multiplicity of the structure by summing the atomic magnetizations.

        multiplicity = 1 + sum_i ( natoms_i * magnetization_i ), for each atom_type i
                     = 1 + sum_i magnetization_j, for each atomic site j

    :param atoms: ASE atoms instance
    :param protocol: protocol dict
    :returns: dict (for cp2k input)
    """
    if protocol['initial_magnetization'] == 'zero':
        # base multiplicity on number of electrons
        # (even atomic number <=> even number of valence electrons)
        is_even = sum(atoms.get_atomic_numbers()) % 2 == 0
        if is_even:
            multiplicity = 1
        else:
            multiplicity = 2
    else:
        # base multiplicity on starting magnetization
        multiplicity = 1 + sum([atom.magmom for atom in atoms])
        multiplicity = int(round(multiplicity))

    multiplicity_dict = {'FORCE_EVAL': {'DFT': {'MULTIPLICITY': multiplicity}}}
    if multiplicity != 1:
        multiplicity_dict['FORCE_EVAL']['DFT']['UKS'] = True
    return multiplicity_dict


def get_kinds_section(atoms, protocol):
    """ Write the &KIND sections given the structure and the settings_dict"""
    kinds_info = get_kinds_info(atoms)

    kinds = [{
        '_': k['kind'],
        'ELEMENT': k['element'],
        'BASIS_SET': protocol['basis_set'][k['element']],
        'POTENTIAL': protocol['pseudopotential'][k['element']],
        'MAGNETIZATION': k['magnetization'],
    } for k in kinds_info]

    return {'FORCE_EVAL': {'SUBSYS': {'KIND': kinds}}}


def get_kinds_with_ghost_section(atoms, protocol):
    """Write the &KIND sections given the structure and the settings_dict, and add also GHOST atoms"""
    kinds_info = get_kinds_info(atoms)

    kinds = []
    for kind_info in kinds_info:
        kinds.append({
            '_': kind_info['kind'],
            'ELEMENT': kind_info['element'],
            'BASIS_SET': protocol['basis_set'][kind_info['element']],
            'POTENTIAL': protocol['pseudopotential'][kind_info['element']],
            'MAGNETIZATION': kind_info['magnetization'],
        })
        kinds.append({
            '_': kind_info['kind'] + '_ghost',
            'ELEMENT': kind_info['element'],
            'BASIS_SET': protocol['basis_set'][kind_info['element']],
            'GHOST': True
        })

    return {'FORCE_EVAL': {'SUBSYS': {'KIND': kinds}}}


def get_bsse_section(natoms_a, natoms_b, mult_a=1, mult_b=1, charge_a=0, charge_b=0):  # pylint: disable=too-many-arguments
    """Get the &FORCE_EVAL/&BSSE section."""
    bsse_section = {
        'FORCE_EVAL': {
            'BSSE' : {
                'FRAGMENT': [{
                    'LIST': '1..{}'.format(natoms_a)
                },
                {
                    'LIST': '{}..{}'.format(natoms_a + 1, natoms_a + natoms_b)
                }],
                'CONFIGURATION': [
                    { # A fragment with basis set A
                        'MULTIPLICITY': mult_a,
                        'CHARGE': charge_a,
                        'GLB_CONF': '1 0',
                        'SUB_CONF': '1 0',
                        },
                    { # B fragment with basis set B
                        'MULTIPLICITY': mult_b,
                        'CHARGE': charge_b,
                        'GLB_CONF': '0 1',
                        'SUB_CONF': '0 1',
                        },
                    { # A fragment with basis set A+B
                        'MULTIPLICITY': mult_a,
                        'CHARGE': charge_a,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '1 0',
                        },
                    { # B fragment with basis set A+B
                        'MULTIPLICITY': mult_b,
                        'CHARGE': charge_b,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '0 1',
                        },
                    { # A+B fragments with basis set A+B
                        'MULTIPLICITY': mult_a + mult_b - 1,
                        'CHARGE': charge_a + charge_b,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '1 1',
                        }
                ]
            }
        }
    }
    return bsse_section


# Functions to parse results


def ot_has_small_bandgap(cp2k_input, cp2k_output, bandgap_thr_ev):
    """ Returns True if the calculation used OT and had a smaller bandgap then the guess needed for the OT.
    (NOTE: It has been observed also negative bandgap with OT in CP2K!)
    cp2k_input: dict
    cp2k_output: dict
    bandgap_thr_ev: float [eV]
    """
    list_true = [True, 'T', 't', '.TRUE.', 'True', 'true']  #add more?
    try:
        ot_settings = cp2k_input['FORCE_EVAL']['DFT']['SCF']['OT']
        if '_' not in ot_settings.keys() or ot_settings['_'] in list_true:  #pylint: disable=simplifiable-if-statement
            using_ot = True
        else:
            using_ot = False
    except KeyError:
        using_ot = False
    min_bandgap_ev = min(cp2k_output['bandgap_spin1_au'], cp2k_output['bandgap_spin2_au']) * HARTREE2EV
    is_bandgap_small = (min_bandgap_ev < bandgap_thr_ev)
    return using_ot and is_bandgap_small
