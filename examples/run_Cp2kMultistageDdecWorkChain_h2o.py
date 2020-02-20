#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
""" Test/example for the DdecCp2kChargesWorkChain"""

import click
import ase.build

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain object
MultistageDdecWorkChain = WorkflowFactory('lsmo.multistageddec')  # pylint: disable=invalid-name

# Data objects
StructureData = DataFactory('structure')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('cp2k_code_string')
@click.argument('ddec_code_string')
@click.argument('ddec_atdens_path')
def main(cp2k_code_string, ddec_code_string, ddec_atdens_path):
    """Example usage:
    ATDENS_PATH='/home/daniele/Programs/aiida-database/data/chargemol_09_26_2017/atomic_densities/'
    verdi run run_Cp2kMultistageDdecWorkChain_h2o.py cp2k@localhost ddec@localhost $ATDENS_PATH
    """
    print('Testing CP2K-Multistage calculation + DDEC on H2O...')

    cp2k_code = Code.get_from_string(cp2k_code_string)
    ddec_code = Code.get_from_string(ddec_code_string)

    atoms = ase.build.molecule('H2O')
    atoms.center(vacuum=2.0)
    structure = StructureData(ase=atoms)

    cp2k_options = {
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 10 * 60,
        'withmpi': True,
    }

    ddec_options = {
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 10 * 60,
        'withmpi': False,
    }

    ddec_params = Dict(
        dict={
            'net charge': 0.0,
            'charge type': 'DDEC6',
            'periodicity along A, B, and C vectors': [True, True, True],
            'compute BOs': False,
            'atomic densities directory complete path': ddec_atdens_path,
            'input filename': 'valence_density',
        })

    inputs = {
        'structure': structure,
        'metadata': {
            'label': 'test-h2o'
        },
        'protocol_tag': Str('test'),
        'cp2k_base': {
            'cp2k': {
                'code': cp2k_code,
                'metadata': {
                    'options': cp2k_options,
                }
            }
        },
        'ddec': {
            'parameters': ddec_params,
            'code': ddec_code,
            'metadata': {
                'options': ddec_options,
            }
        }
    }

    run(MultistageDdecWorkChain, **inputs)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
