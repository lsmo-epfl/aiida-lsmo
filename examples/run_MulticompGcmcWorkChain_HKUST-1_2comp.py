#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example 2-components GCMC in HKUST-1, at 3 different T/P conditions. Computing blocking spheres."""
import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict
from aiida import cmdline

# Workchain objects
MulticompGcmcWorkChain = WorkflowFactory('lsmo.multicomp_gcmc')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('core.cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name
SinglefileData = DataFactory('core.singlefile')


@click.command('cli')
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa_code', type=cmdline.params.types.CodeParamType())
@click.option('--zeopp_code', type=cmdline.params.types.CodeParamType())
def main(raspa_code, zeopp_code):
    """Prepare inputs and submit the workchain.
    Usage: verdi run run_thisworkchainexample.py raspa@localhost zeopp@localhost"""

    builder = MulticompGcmcWorkChain.get_builder()

    builder.metadata.label = 'test'

    builder.raspa_base.raspa.code = raspa_code
    builder.zeopp.code = zeopp_code

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
    builder.conditions = Dict({
        'molfraction': {
            'xenon': 0.2,
            'krypton': 0.8,
        },
        'temp_press': [
            [200, 0.1],
            [300, 0.5],
            [400, 1.0],
        ]
    })

    builder.parameters = Dict(
        {
            'zeopp_probe_scaling': 0.95,
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_gcmc_init_cycles': 100,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
