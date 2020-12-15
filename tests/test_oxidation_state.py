# -*- coding: utf-8 -*-
"""Test oxidation state prediction."""
import os
import pytest
from aiida.orm import CifData
from aiida_lsmo.calcfunctions.oxidation_state import compute_oxidation_states, tag_ase_atoms, OXIMACHINE_RUNNER

from . import DATA_DIR

CIF_FILES = ['Fe-MOF-74.cif', 'Mg_MOF_74.cif', 'Cu-I-II-HKUST-1.cif']
OXIDATION_STATES = {
    'Fe-MOF-74.cif': [2, 2, 2, 2, 2, 2],
    'Mg_MOF_74.cif': [2, 2, 2, 2, 2, 2],
    'Cu-I-II-HKUST-1.cif': [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
}


@pytest.mark.parametrize('cif_filename', CIF_FILES)
def test_oxidation_state(cif_filename):
    """Check that oxidation states are as expected."""
    with open(os.path.join(DATA_DIR, cif_filename), 'rb') as handle:
        cif = CifData(file=handle)

    results_dict = compute_oxidation_states(cif).get_dict()
    assert results_dict['prediction'] == OXIDATION_STATES[cif_filename]


def test_tag_atoms():
    """Test that atoms are correctly tagged and assigned charges."""
    from ase.io import read
    atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
    atoms = tag_ase_atoms(atoms=atoms, oxi_results=OXIMACHINE_RUNNER.run_oximachine(atoms))

    assert atoms[113].charge == -1  # Cu(I)
    assert atoms[113].tag == 1  # Cu(I)
    assert atoms[140].charge == -2  # Cu(II)
    assert atoms[140].tag == 2  # Cu(II)
