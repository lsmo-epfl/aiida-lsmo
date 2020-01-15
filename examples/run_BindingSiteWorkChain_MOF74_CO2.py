#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example binding site work chain for CO2 in Zn-MOF-74 framework."""

from __future__ import absolute_import
from __future__ import print_function

import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain objects
BindingSiteWorkChain = WorkflowFactory('lsmo.binding_site')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
@click.argument('cp2k_code_label')
def main(raspa_code_label, cp2k_code_label):
    """Prepare inputs and submit the work chain."""

    print("Testing BindingSite work chain (FF + DFT) for CO2 in Zn-MOF-74 ...")
    print("[NOTE: this test will run on 4 cpus and take ca. 10 minutes]")

    builder = BindingSiteWorkChain.get_builder()
    builder.metadata.label = "test"
    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.cp2k_base.cp2k.code = Code.get_from_string(cp2k_code_label)
    builder.raspa_base.raspa.metadata.options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 10 * 60,
    }
    builder.cp2k_base.cp2k.metadata.options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 4,
        },
        "max_wallclock_seconds": 1 * 10 * 60,
    }
    builder.structure = CifData(file=os.path.abspath('data/Zn-MOF-74.cif'), label="Zn-MOF-74")
    builder.molecule = Str('co2')
    builder.parameters = Dict(
        dict={
            "ff_framework": "UFF",  # (str) Forcefield of the structure.
            "mc_steps": int(10),  # (int) Number of MC cycles.
            "temperature_list": [300, 150],
        })
    builder.protocol_tag = Str('test')
    builder.cp2k_base.cp2k.parameters = Dict(dict={ # Lowering CP2K default setting for a faster test calculation
        'FORCE_EVAL': {
            'DFT': {
                'SCF': {
                    'EPS_SCF': 1.0E-4,
                    'OUTER_SCF': {
                        'EPS_SCF': 1.0E-4,
                    },
                },
            },
        },
        'MOTION': {
            'GEO_OPT': {
                'MAX_ITER': 5
            }
        },
    })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
