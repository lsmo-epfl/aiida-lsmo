#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run Cp2kMultistageWorkChain workchain on Aluminum."""

from pathlib import Path
import pytest
import click
import ase.build

from aiida.plugins import DataFactory, WorkflowFactory
from aiida import cmdline
from aiida import engine
from aiida.orm import Dict, StructureData, Str, Int, SinglefileData

THIS_DIR = Path(__file__).resolve().parent
DATA_DIR = THIS_DIR / 'data'

# Workchain objects
Cp2kMultistageWorkChain = WorkflowFactory('lsmo.cp2k_multistage')

# Data objects
StructureData = DataFactory('core.structure')


@pytest.fixture(scope='function')
def al_structuredata():
    """StructureData for Aluminum."""
    return StructureData(ase=ase.io.read(DATA_DIR / 'Al.cif'))


def run_multistage_al(cp2k_code, al_structuredata):  # pylint: disable=redefined-outer-name
    """Run Cp2kMultistageWorkChain workchain on Aluminum."""

    print('Testing CP2K multistage workchain on Al (RKS, needs smearing)...')
    print('EXPECTED: the OT (settings_0) will converge to a negative bandgap, then use SMEARING (settings_1)')

    # testing user change of parameters and protocol
    parameters = Dict({'FORCE_EVAL': {'DFT': {'MGRID': {'CUTOFF': 250,}}}})
    protocol_mod = Dict({
        'initial_magnetization': {
            'Al': 0
        },
        'settings_0': {
            'FORCE_EVAL': {
                'DFT': {
                    'SCF': {
                        'OUTER_SCF': {
                            'MAX_SCF': 5,
                        }
                    }
                }
            }
        }
    })

    # Construct process builder
    builder = Cp2kMultistageWorkChain.get_builder()
    builder.structure = al_structuredata
    builder.protocol_tag = Str('test')
    builder.starting_settings_idx = Int(0)
    builder.protocol_modify = protocol_mod
    builder.cp2k_base.cp2k.parameters = parameters
    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 1,
    }
    builder.cp2k_base.cp2k.metadata.options.withmpi = False  # comment this for parallel cp2k executable
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 3 * 60
    # The following is not needed, if the files are available in the data directory of your CP2K executable
    cp2k_dir = DATA_DIR / 'cp2k'
    builder.cp2k_base.cp2k.file = {
        'basis': SinglefileData(file=str(cp2k_dir / 'BASIS_MOLOPT')),
        'pseudo': SinglefileData(file=str(cp2k_dir / 'GTH_POTENTIALS')),
        'dftd3': SinglefileData(file=str(cp2k_dir / 'dftd3.dat')),
    }

    results, node = engine.run_get_node(builder)

    assert node.is_finished_ok, results
    output_parameters = results['output_parameters'].get_dict()
    assert output_parameters['final_bandgap_spin1_au'] == pytest.approx(0, abs=1e-3)
    assert output_parameters['final_bandgap_spin2_au'] == pytest.approx(0, abs=1e-3)
    assert output_parameters['step_info']['scf_converged'][-1]


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
def cli(cp2k_code):
    """Run example.

    Example usage: $ ./test_Cp2kMultistageWorkChain_Al.py --cp2k-code my-cp2k@myhost
    """
    run_multistage_al(cp2k_code, al_structuredata=StructureData(ase=ase.io.read(DATA_DIR / 'Al.cif')))


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
