#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example PE calculation calculation with HKUST1 framework."""

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict

# Workchain objects
IsothermCalcPEWorkChain = WorkflowFactory('lsmo.isotherm_calc_pe')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the Isotherm workchain.
    Usage: verdi run run_IsothremCalcPE_HKUST-1.py raspa@localhost zeopp@localhost"""

    builder = IsothermCalcPEWorkChain.get_builder()

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

    builder.structure = CifData(file=os.path.abspath('data/HKUST-1.cif'), label="HKUST-1")
    # NOTE1: parameters are chosen for speed purpose. Use the default to be consistent to 10.1021/acscentsci.9b00619
    # NOTE2: calc_pe can fail due to this raw sampling of the isotherm. It happens in the ca. 20% of the cases.
    builder.parameters = Dict(
        dict={
            'ff_framework': 'UFF',
            'temperature': 400,
            'zeopp_volpo_samples': 1000,
            'zeopp_block_samples': 10,
            'raspa_widom_cycles': 100,
            'raspa_gcmc_init_cycles': 100,
            'raspa_gcmc_prod_cycles': 2000,
            'pressure_min': 0.001,
            'pressure_max': 10,
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
