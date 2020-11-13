#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example sim annealing of CO2 in Zn-MOF-74 framework."""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain objects
SimAnnealingWorkChain = WorkflowFactory('lsmo.sim_annealing')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
def main(raspa_code_label):
    """Prepare inputs and submit the work chain."""

    builder = SimAnnealingWorkChain.get_builder()
    builder.metadata.label = 'test'
    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.raspa_base.raspa.metadata.options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
    }
    builder.structure = CifData(file=os.path.abspath('data/Zn-MOF-74.cif'), label='Zn-MOF-74')
    builder.molecule = Str('co2')
    builder.parameters = Dict(dict={
        'ff_framework': 'UFF',  # (str) Forcefield of the structure.
        'mc_steps': int(10),  # (int) Number of MC cycles.
    })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
