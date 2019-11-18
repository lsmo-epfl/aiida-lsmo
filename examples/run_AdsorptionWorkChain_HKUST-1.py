#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example two-component isotherm calculation with HKUST1 framework."""

from __future__ import absolute_import
from __future__ import print_function

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory  #WorkflowFactory
from aiida.orm import Code, Dict
from aiida_lsmo.workchains import AdsorptionWorkChain

# Workchain objects
#IsothermMultiCompWorkChain = WorkflowFactory('lsmo.isotherm_multi_comp')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name
SinglefileData = DataFactory('singlefile')


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the Isotherm workchain.
    Usage: verdi run run_IsothermMultiCompWorkChain_HKUST-1.py raspa@localhost network@localhost"""

    builder = AdsorptionWorkChain.get_builder()

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
    builder.mixture = Dict(dict={
        'comp1': {
            'name': 'xenon',
            'molfraction': 0.2
        },
        'comp2': {
            'name': 'krypton',
            'molfraction': 0.8
        }
    })

    builder.parameters = Dict(
        dict={
            'forcefield': 'UFF',  # Default: UFF
            'temperature': 298,  # (K) Note: higher temperature will have less adsorbate and it is faster
            'zeopp_volpo_samples': 10,  # Default: 1e5 *NOTE: default is good for standard real-case!
            'zeopp_sa_samples': 10,  # Default: 1e5 *NOTE: default is good for standard real-case!
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 100,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
            'pressure_list': [
                1.0, 2.0
            ],  # Comment this line and uncomment the following three to have linear range of pressure points.
            # "pressure_precision": 0.5,
            # "pressure_min": 0.1,
            # "pressure_max": 2.0
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
