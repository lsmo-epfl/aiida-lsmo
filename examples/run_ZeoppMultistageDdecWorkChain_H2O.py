#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Submit ZeoppMultistageDdecWorkChain for H2O"""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain objects
ZeoppMultistageDdecWorkChain = WorkflowFactory('lsmo.zeopp_multistage_ddec')  # pylint: disable=invalid-name

#Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('zeopp_code_string')
@click.argument('cp2k_code_string')
@click.argument('ddec_code_string')
def main(zeopp_code_string, cp2k_code_string, ddec_code_string, ddec_atdens_path):
    """Example usage:
    ATDENS_PATH='/home/daniele/Programs/aiida-database/data/chargemol_09_26_2017/atomic_densities/'
    verdi run run_ZeoppMultistageDdecWorkChain_H2O.py zeopp@localhost cp2k@localhost ddec@localhost $ATDENS_PATH
    """
    zeopp_code = Code.get_from_string(zeopp_code_string)
    cp2k_code = Code.get_from_string(cp2k_code_string)
    ddec_code = Code.get_from_string(ddec_code_string)

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

    zeopp_options = {
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

    zeopp_params = NetworkParameters(
        dict={
            'ha': 'DEF',  # Using high accuracy (mandatory!)
            'res': True,  # Max included, free and incl in free sphere
            'sa': [1.86, 1.86, 1000],  # Nitrogen probe to compute surface
            'vol': [0.0, 0.0, 1000],  # Geometric pore volume
        })

    structure = CifData(file=os.path.join(os.getcwd(), 'data/H2O.cif')).store()
    structure.label = 'H2O'

    inputs = {
        'structure': structure,
        'protocol_tag': Str('test'),
        'metadata': {
            'label': 'test',
        },
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
        },
        'zeopp': {
            'parameters': zeopp_params,
            'code': zeopp_code,
            'metadata': {
                'options': zeopp_options,
            }
        }
    }

    run(ZeoppMultistageDdecWorkChain, **inputs)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
