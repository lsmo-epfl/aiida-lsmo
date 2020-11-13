#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example 3-components GCMC in HKUST-1, at 3 different T/P conditions. Skipping calculation of blocking spheres."""
import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict

# Workchain objects
MulticompGcmcWorkChain = WorkflowFactory('lsmo.multicomp_gcmc')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name
SinglefileData = DataFactory('singlefile')


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the workchain.
    Usage: verdi run run_thisworkchainexample.py raspa@localhost zeopp@localhost"""

    builder = MulticompGcmcWorkChain.get_builder()

    builder.metadata.label = 'test'

    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.zeopp.code = Code.get_from_string(zeopp_code_label)

    options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
        'withmpi': False,
    }
    builder.raspa_base.raspa.metadata.options = options
    builder.zeopp.metadata.options = options
    builder.structure = CifData(file=os.path.abspath('data/HKUST-1.cif'), label='HKUST-1')
    builder.conditions = Dict(dict={
        'molfraction': {
            'co': 0.2,
            'ethene': 0.3,
            'ethane': 0.5,
        },
        'tp_gcmc': [
            [200, 0.1],
            [300, 0.5],
            [400, 0.7],
        ]
    })

    builder.parameters = Dict(
        dict={
            'zeopp_probe_scaling': 0.0,  # NOTE: excluding blocking sphere calculations
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_gcmc_init_cycles': 100,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
