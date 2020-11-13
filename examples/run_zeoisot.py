# -*- coding: utf-8 -*-
"""Run example isotherm calculation with COF-1 framework."""

from __future__ import absolute_import
from __future__ import print_function
import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory
from aiida.orm import Code, Dict, Str, Float, Int

# Workchain objects
#IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  # pylint: disable=invalid-name
from aiida_lsmo.workchains.zeoisotherm import ZeoIsothermWorkChain

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the Isotherm workchain.
    Usage: verdi run run_IsothermWorkChain_COF-1.py raspa@localhost network@localhost"""

    builder = ZeoIsothermWorkChain.get_builder()

    builder.metadata.label = "test"  # pylint: disable=no-member

    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)  # pylint: disable=no-member
    builder.zeopp.code = Code.get_from_string(zeopp_code_label)  # pylint: disable=no-member

    options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 60 * 60,
        "withmpi": False,
    }

    builder.raspa_base.raspa.metadata.options = options  # pylint: disable=no-member
    builder.zeopp.metadata.options = options  # pylint: disable=no-member
    builder.raspa_base.max_iterations = Int(1)  # pylint: disable=no-member

    builder.structure = CifData(file=os.path.abspath('data/zeo_LTA_p91.cif'), label="zeo_LTA_p91")
    builder.sial_ratio = Float(1.101)
    builder.cation = Str('Na')

    builder.molecule = Str('co2')
    builder.parameters = Dict(
        dict={
            'ff_framework': 'Zeo-test',  # Default: UFF
            'temperature': 300,  # (K) Note: higher temperature will have less adsorbate and it is faster
            'zeopp_volpo_samples': 1000,  # Default: 1e5 *NOTE: default is good for standard real-case!
            'zeopp_block_samples': 10,  # Default: 100
            "raspa_nvt_cations_cycles": 1000,
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 10,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
            'pressure_min': 0.001,  # Default: 0.001 (bar)
            'pressure_max': 3,  # Default: 10 (bar)
        })

    run(builder)


if __name__ == '__main__':
    print('Testing Zeolite Isotherm work chain.')
    main()  # pylint: disable=no-value-for-parameter

# EOF
