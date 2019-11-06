#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example PE calculation calculation with HKUST1 framework."""

from __future__ import absolute_import
from __future__ import print_function

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
    builder.parameters = Dict(
        dict={
            'forcefield': 'UFF',
            'temperature': 400,  # WARNING: use 300 K for real PE calculation!
            'zeopp_volpo_samples': 1000,  # Default: 1e5
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 100,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 1000,  # Default: 1e4
            'pressure_min': 0.001,  # Default: 0.001 (bar)
            'pressure_max': 10,  # WARNING: use 30 bar for real PE calculation!
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
