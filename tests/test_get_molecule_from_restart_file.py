#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""test calcfunction in SimAnnealing WC to extract cif from restart file"""

from pathlib import Path
from aiida.plugins import DataFactory
from aiida.orm import Dict
from aiida_lsmo.workchains.sim_annealing import get_molecule_from_restart_file

THIS_DIR = Path(__file__).resolve().parent
DATA_DIR = THIS_DIR / 'data'

# Data objects
CifData = DataFactory('core.cif')  # pylint: disable=invalid-name
FolderData = DataFactory('core.folder')


def test_get_molecule_from_restart_file():
    """Check that the molecule cif is extracted correctly from the RestartFile"""

    with open(DATA_DIR / 'HKUST-1.cif', 'rb') as handle:
        structure_cif = CifData(file=handle)

    molecule_folderdata = FolderData(tree=(DATA_DIR / 'RASPA_3N2_in_HKUST1'))

    input_dict = Dict({'number_of_molecules': int(3)})

    molecule_dict = Dict({'name': 'N2', 'forcefield': 'TraPPE'})

    cif_out = get_molecule_from_restart_file(structure_cif, molecule_folderdata, input_dict, molecule_dict)

    atoms_output = cif_out.get_ase().symbols.get_chemical_formula()
    assert atoms_output == str('N6')
