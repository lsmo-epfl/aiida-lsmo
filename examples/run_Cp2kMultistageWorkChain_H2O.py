# -*- coding: utf-8 -*-
""" Test/example for the Cp2kMultistageWorkChain"""

import click
import ase.build

from aiida.engine import run
from aiida.orm import Dict, StructureData, Str
from aiida.plugins import WorkflowFactory
from aiida import cmdline

Cp2kMultistageWorkChain = WorkflowFactory('lsmo.cp2k_multistage')  # pylint: disable=invalid-name


def run_multistage_h2o(cp2k_code):
    """Example usage: verdi run thistest.py cp2k@localhost"""

    print('Testing CP2K multistage workchain on H2O (RKS, no need for smearing)...')
    print(">>> Using 'test' tag")

    atoms = ase.build.molecule('H2O')
    atoms.center(vacuum=2.0)
    structure = StructureData(ase=atoms)

    # Construct process builder
    builder = Cp2kMultistageWorkChain.get_builder()
    builder.structure = structure
    builder.protocol_tag = Str('test')
    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 1,
    }
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 3 * 60
    builder.cp2k_base.cp2k.parameters = Dict(dict={
        'GLOBAL': {
            'PREFERRED_DIAG_LIBRARY': 'SL'
        },
    })

    run(builder)


@click.command('cli')
@cmdline.utils.decorators.with_dbenv()
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
def cli(cp2k_code):
    """Click interface"""
    run_multistage_h2o(cp2k_code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
