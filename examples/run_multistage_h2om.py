""" Test/example for the Cp2kMultistageWorkChain"""

import sys
import click
import ase.build

from aiida.engine import run
from aiida.orm import Code, Dict, StructureData, Float, Str
from aiida.common import NotExistent
from aiida.plugins import WorkflowFactory

Cp2kMultistageWorkChain = WorkflowFactory('lsmo.cp2k_multistage')  # pylint: disable=invalid-name


def run_multistage_h2om(cp2k_code):
    """Example usage: verdi run thistest.py cp2k@localhost"""

    print("Testing CP2K multistage workchain on 2xH2O- (UKS, no need for smearing)...")
    print("This is checking:")
    print(" > unit cell resizing")
    print(" > protocol modification")
    print(" > cp2k calc modification")

    atoms = ase.build.molecule('H2O')
    atoms.center(vacuum=2.0)
    structure = StructureData(ase=atoms)

    protocol_mod = Dict(dict={'settings_0': {'FORCE_EVAL': {'DFT': {'MGRID': {'CUTOFF': 300,}}}}})
    parameters = Dict(dict={'FORCE_EVAL': {
        'DFT': {
            'UKS': True,
            'MULTIPLICITY': 3,
            'CHARGE': -2,
        }
    }})

    # Construct process builder
    builder = Cp2kMultistageWorkChain.get_builder()
    builder.structure = structure
    builder.min_cell_size = Float(4.1)
    builder.protocol_tag = Str('test')
    builder.protocol_modify = protocol_mod

    builder.cp2k_base.cp2k.parameters = parameters
    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 3 * 60

    run(builder)


@click.command('cli')
@click.argument('codelabel')
def cli(codelabel):
    """Click interface"""
    try:
        code = Code.get_from_string(codelabel)
    except NotExistent:
        print("The code '{}' does not exist".format(codelabel))
        sys.exit(1)
    run_multistage_h2om(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
