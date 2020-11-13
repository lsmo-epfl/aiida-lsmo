#!/usr/bin/env python  # pylint: disable=invalid-name
"""Run example IsothermInflection for Ar in Graphite."""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict

# Workchain objects
IsothermInflectionWorkChain = WorkflowFactory('lsmo.isotherm_inflection')

# Data objects
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the workchain.
    Usage: verdi run run_thisworkchainexample.py raspa@localhost zeopp@localhost"""

    builder = IsothermInflectionWorkChain.get_builder()

    builder.metadata.label = "test"

    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.zeopp.code = Code.get_from_string(zeopp_code_label)

    options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 60 * 60,
        "withmpi": False,
    }
    builder.raspa_base.raspa.metadata.options = options
    builder.zeopp.metadata.options = options
    builder.structure = CifData(file=os.path.abspath('data/Graphite_20A.cif'), label="Graphite_20A")

    builder.molecule = Dict(
        dict={
            'name': 'Ar',
            'forcefield': 'HIRSCHFELDER',
            "ff_cutoff": 8,
            'molsatdens': 35.4,
            'proberad': 1.7,
            'singlebead': True,
            'charged': False,
            'rosenbluth': 1,
            'pressure_zero': 1,
        })

    builder.parameters = Dict(
        dict={
            'ff_framework': 'DREIDING',
            'temperature': 87,  # Tsat Ar
            "ff_cutoff": 8.0,  # NOTE: Low to have cheap testing
            'box_length': 16.0,
            "zeopp_probe_scaling": 1.0,
            'zeopp_volpo_samples': 10000,
            'zeopp_block_samples': 100,
            'raspa_widom_cycles': 1000,
            'raspa_gcmc_init_cycles': 200,
            'raspa_gcmc_prod_cycles': 200,
            "pressure_num": 4,
            "raspa_verbosity": 10
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
