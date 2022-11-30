#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example sim annealing of 3 CO2 molecules in HKUST-1 framework."""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict
from aiida import cmdline

# Workchain objects
SimAnnealingWorkChain = WorkflowFactory('lsmo.sim_annealing')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


@click.command('cli')
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa_code', type=cmdline.params.types.CodeParamType())
def main(raspa_code):
    """Prepare inputs and submit the work chain."""

    builder = SimAnnealingWorkChain.get_builder()
    builder.metadata.label = 'test'
    builder.raspa_base.raspa.code = raspa_code
    builder.raspa_base.raspa.metadata.options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
    }
    builder.structure = CifData(file=os.path.abspath('data/HKUST-1.cif'), label='HKUST-1')
    builder.molecule = Dict(dict={
        'name': 'DMA',
        'forcefield': 'BOYDDDEC',
        'proberad': 1.7365,
        'singlebead': False,
        'charged': True,
    })
    builder.parameters = Dict(
        dict={
            'ff_framework': 'UFF',  # (str) Forcefield of the structure.
            'temperature_list': [300, 150],  # (list) List of decreasing temperatures for the annealing.
            'mc_steps': int(10),  # (int) Number of MC cycles.
            'number_of_molecules': 1  # (int) Number of molecules loaded in the framework.
            # 'reinsertion_probability': float(0.0),
            # 'randomtranslation_probability': float(1.0),
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
