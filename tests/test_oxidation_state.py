# -*- coding: utf-8 -*-
"""Test oxidation state prediction."""
import os
from collections import OrderedDict

import pytest
from numpy import array
from ase.io import read
from aiida.orm import CifData

from aiida_lsmo.calcfunctions.oxidation_state import compute_oxidation_states
from aiida_lsmo.workchains.cp2k_multistage_protocols import set_initial_conditions, load_isotherm_protocol
from aiida_lsmo.utils.cp2k_utils import get_kinds_section, get_multiplicity_section
from . import DATA_DIR

CIF_FILES = ['Fe-MOF-74.cif', 'Mg_MOF_74.cif', 'Cu-I-II-HKUST-1.cif']
OXIDATION_STATES = {
    'Fe-MOF-74.cif': [2, 2, 2, 2, 2, 2],
    'Mg_MOF_74.cif': [2, 2, 2, 2, 2, 2],
    'Cu-I-II-HKUST-1.cif': [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
}


@pytest.fixture()
def cu_hkust_oxidation():
    """Oxidation state prediction for Cu-HKUST-1

    Cached to speed up tests.
    """
    return OrderedDict([
        ('metal_indices', [112, 113, 114, 115, 116, 117, 118, 119, 120, 124, 128, 132, 136, 140, 144, 148]),
        ('metal_symbols',
         ['Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu', 'Cu']),
        ('prediction', [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2]),
        ('max_probas', [
            0.9944271144788622, 0.9944271144788622, 0.9944271144788622, 0.9944271144788622, 0.9944271144788622,
            0.9944271144788622, 0.9944271144788622, 0.9944271144788622, 0.9935515666748629, 0.9935515666748629,
            0.9935515666748629, 0.9935515666748629, 0.9935515666748629, 0.9935515666748629, 0.9935515666748629,
            0.9935515666748629
        ]),
        ('base_predictions', [
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([0, 0, 0, 0]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1]),
            array([1, 1, 1, 1])
        ])
    ])


@pytest.mark.parametrize('cif_filename', CIF_FILES)
def test_oxidation_state(cif_filename):
    """Check that oxidation states are as expected."""
    with open(os.path.join(DATA_DIR, cif_filename), 'rb') as handle:
        cif = CifData(file=handle)

    results_dict = compute_oxidation_states(cif).get_dict()
    assert results_dict['prediction'] == OXIDATION_STATES[cif_filename]


def test_initial_conditions(cu_hkust_oxidation):  # pylint: disable=redefined-outer-name
    """Test that atoms are correctly tagged and assigned charges and initial magnetic moments."""

    # test oxidation_state
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    atoms = set_initial_conditions(atoms=atoms,
                                   initial_magnetization='oxidation_state',
                                   oxidation_states=cu_hkust_oxidation)
    assert atoms[113].charge == -1  # Cu(I)
    assert atoms[113].tag == 1  # Cu(I)
    assert atoms[113].magmom == 0.0  # Cu(I)
    assert atoms[140].charge == -2  # Cu(II)
    assert atoms[140].tag == 2  # Cu(II)
    assert atoms[140].magmom == 1.0  # Cu(I)

    # test element
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    atoms = set_initial_conditions(atoms=atoms, initial_magnetization='element')
    assert atoms[113].charge == -2  # Cu(I)
    assert atoms[113].tag == 0  # Cu(I)
    assert atoms[113].magmom == 1.0  # Cu(I)
    assert atoms[140].charge == -2  # Cu(II)
    assert atoms[140].tag == 0  # Cu(II)
    assert atoms[140].magmom == 1.0  # Cu(I)

    # test zero
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    atoms = set_initial_conditions(atoms=atoms, initial_magnetization='zero')
    assert atoms[113].charge == 0  # Cu(I)
    assert atoms[113].tag == 0  # Cu(I)
    assert atoms[113].magmom == 0.0  # Cu(I)
    assert atoms[140].charge == 0  # Cu(II)
    assert atoms[140].tag == 0  # Cu(II)
    assert atoms[140].magmom == 0.0  # Cu(I)

    # test manual (magnetization only)
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    mags = {'Cu': 2.0}
    atoms = set_initial_conditions(atoms=atoms, initial_magnetization=mags)
    assert atoms[113].charge == 0  # Cu(I)
    assert atoms[113].tag == 0  # Cu(I)
    assert atoms[113].magmom == 2.0  # Cu(I)
    assert atoms[140].charge == 0  # Cu(II)
    assert atoms[140].tag == 0  # Cu(II)
    assert atoms[140].magmom == 2.0  # Cu(I)

    # test manual (oxidation state and magnetization)
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    mags = {
        'Cu': {
            'magnetization': {
                2: 1.0
            },
            'default_oxidation': 2,
        }
    }
    atoms = set_initial_conditions(atoms=atoms, initial_magnetization=mags)
    assert atoms[113].charge == -2  # Cu(I)
    assert atoms[113].tag == 0  # Cu(I)
    assert atoms[113].magmom == 1.0  # Cu(I)
    assert atoms[140].charge == -2  # Cu(II)
    assert atoms[140].tag == 0  # Cu(II)
    assert atoms[140].magmom == 1.0  # Cu(I)


def test_cp2k_kinds(cu_hkust_oxidation):  # pylint: disable=redefined-outer-name
    """"Test the CP2K kinds"""

    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    protocol = load_isotherm_protocol(tag='standard')
    atoms = set_initial_conditions(atoms=atoms,
                                   initial_magnetization=protocol['initial_magnetization'],
                                   oxidation_states=cu_hkust_oxidation)

    kinds_section = get_kinds_section(atoms=atoms, protocol=protocol)
    assert len(kinds_section['FORCE_EVAL']['SUBSYS']['KIND']) == len(set(
        atoms.get_chemical_symbols())) + 1, kinds_section['FORCE_EVAL']['SUBSYS']['KIND']

    multiplicity_section = get_multiplicity_section(atoms=atoms)
    assert multiplicity_section['FORCE_EVAL']['DFT']['MULTIPLICITY'] == 8 * 1 + 1, multiplicity_section
