# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
""" Test/example for the BindingEnergyWorkChain"""

import os
import click
import ase.io

from aiida.engine import run
from aiida.orm import Code, Dict, StructureData, Str
from aiida.plugins import WorkflowFactory

BindingEnergyWorkChain = WorkflowFactory('lsmo.cp2k_binding_energy')


@click.command('cli')
@click.argument('cp2k_code_label')
def main(cp2k_code_label):
    """Example usage: verdi run thistest.py cp2k@localhost"""

    cp2k_code = Code.get_from_string(cp2k_code_label)

    print("Testing CP2K BindingEnergy work chain for CO2 in Zn-MOF-74 ...")
    print("[NOTE: this test will run on 4 cpus and take ca. 10 minutes]")

    thisdir = os.path.dirname(os.path.abspath(__file__))

    # Construct process builder
    builder = BindingEnergyWorkChain.get_builder()
    builder.structure = StructureData(ase=ase.io.read(os.path.join(thisdir, 'data/Zn-MOF-74.cif')))
    builder.molecule = StructureData(ase=ase.io.read(os.path.join(thisdir, 'data/CO2_in_Zn-MOF-74.cif')))
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
    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 4,
    }
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 5 * 60

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
