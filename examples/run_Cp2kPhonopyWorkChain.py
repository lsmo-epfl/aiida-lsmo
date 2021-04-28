# -*- coding: utf-8 -*-
""" Test/example for the Cp2kPhonopyWorkChain"""

import click

from aiida.engine import run
from aiida.orm import load_node
from aiida.plugins import WorkflowFactory, DataFactory
from aiida import cmdline

Cp2kPhonopyWorkChain = WorkflowFactory('lsmo.cp2k_phonopy')
Str = DataFactory('str')  # pylint: disable=invalid-name
Int = DataFactory('int')  # pylint: disable=invalid-name


def run_cp2k_phonopy(cp2k_code, structure_pk):
    """Test the workchain in both serial and parallel mode"""
    for mode in ['serial', 'parallel']:
        print(f'>>> Compute forces + 3 displacements for water - MODE: {mode}')
        builder = Cp2kPhonopyWorkChain.get_builder()
        builder.structure = load_node(structure_pk)
        builder.mode = Str(mode)
        builder.max_displacements = Int(3)  # Compute only few displacements (instead of 6N) for sake of time
        builder.cp2k_base.cp2k.code = cp2k_code
        builder.cp2k_base.cp2k.metadata.options.resources = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
        }
        builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 3 * 60

        run(builder)


@click.command('cli')
@cmdline.utils.decorators.with_dbenv()
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
@click.option('--structure-pk')  # use the pk of a CifData/StructureData which is the descendant of a CP2K calculation
def cli(cp2k_code, structure_pk):
    """Click interface"""
    run_cp2k_phonopy(cp2k_code, structure_pk)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
