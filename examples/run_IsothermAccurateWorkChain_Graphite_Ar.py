#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example isotherm-accurate work chain with Graphite-7A framework."""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain objects
IsothermAccurateWorkChain = WorkflowFactory('lsmo.isotherm_accurate')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the Isotherm workchain.
    Usage: verdi run run_isotherm_hkust1.py raspa@localhost network@localhost"""

    builder = IsothermAccurateWorkChain.get_builder()

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

    builder.structure = CifData(file=os.path.abspath('data/Graphite_7A.cif'), label='Graphite_7A')
    builder.molecule = Dict(
        dict={
            'name': 'Ar',
            'forcefield': 'HIRSCHFELDER',
            'molsatdens': 35.4,
            'proberad': 1.7,
            'singlebead': True,
            'charged': False,
        })

    builder.parameters = Dict(
        dict={
            'ff_framework': 'UFF',  # Default: UFF
            'ff_cutoff': 10,  # To speed up
            'temperature': 400,  # (K) Note: higher temperature will have less adsorbate and it is faster
            'zeopp_probe_scaling': 0.8,
            'zeopp_volpo_samples': 1000,  # Default: 1e5 *NOTE: default is good for standard real-case!
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 10,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
