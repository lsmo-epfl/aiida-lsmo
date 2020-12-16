# -*- coding: utf-8 -*-
"""Test oxidation state prediction."""
import os
import pytest
from aiida.orm import CifData, StructureData
from aiida_lsmo.calcfunctions.oxidation_state import compute_oxidation_states, OXIMACHINE_RUNNER

from aiida_lsmo.utils.cp2k_utils import Cp2kSubsys
from aiida_lsmo.workchains.cp2k_multistage_protocols import load_isotherm_protocol
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


# def test_tag_atoms():
#     """Test that atoms are correctly tagged and assigned charges."""
#     from ase.io import read
#     atoms = read(os.path.join(DATA_DIR, 'Cu-I-II-HKUST-1.cif'))
#     atoms = tag_ase_atoms(atoms=atoms, oxi_results=OXIMACHINE_RUNNER.run_oximachine(atoms))
#
#     assert atoms[113].charge == -1  # Cu(I)
#     assert atoms[113].tag == 1  # Cu(I)
#     assert atoms[140].charge == -2  # Cu(II)
#     assert atoms[140].tag == 2  # Cu(II)


def test_cp2k_subsys():
    """"Test the CP2K subsys features"""
    from ase.io import read
    atoms = read(os.path.join(DATA_DIR, 'Fe-MOF-74.cif'))
    structure = StructureData(ase=atoms)

    subsys = Cp2kSubsys(structure=structure,
                        oxidation_states=OXIMACHINE_RUNNER.run_oximachine(atoms),
                        protocol=load_isotherm_protocol(tag='standard'))

    subsys.tag_atoms()
    assert len(subsys.get_kinds_section()['FORCE_EVAL']['SUBSYS']['KIND']) == 4
    assert subsys.get_multiplicity_section()['FORCE_EVAL']['DFT']['MULTIPLICITY'] == 25
