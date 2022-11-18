#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Test/example for the BindingEnergyWorkChain"""

from pathlib import Path
import pytest
import click
import ase.build

from aiida.plugins import DataFactory, WorkflowFactory
from aiida import cmdline
from aiida import engine
from aiida.orm import Dict, StructureData, Str, SinglefileData, WorkChainNode, CalcJobNode

THIS_DIR = Path(__file__).resolve().parent
DATA_DIR = THIS_DIR / 'data'

# Workchain objects
BindingEnergyWorkChain = WorkflowFactory('lsmo.cp2k_binding_energy')

# Data objects
StructureData = DataFactory('structure')


@pytest.fixture(scope='function')
def zn_mof74():
    """StructureData for MOF 74."""
    return StructureData(ase=ase.io.read(DATA_DIR / 'Zn-MOF-74.cif'))


@pytest.fixture(scope='function')
def co2_in_mof74():
    """StructureData for CO2 in Zn-MOF-74."""
    return StructureData(ase=ase.io.read(DATA_DIR / 'CO2_in_Zn-MOF-74.cif'))


def run_binding_energy_co2_mof74(cp2k_code, zn_mof74, co2_in_mof74):  # pylint: disable=redefined-outer-name
    """Compute binding energy of CO2 in MOF 74"""

    print('Testing CP2K BindingEnergy work chain for CO2 in Zn-MOF-74 ...')

    # Construct process builder
    builder = BindingEnergyWorkChain.get_builder()
    builder.structure = zn_mof74
    builder.molecule = co2_in_mof74
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
        'num_machines': 1,
        'num_mpiprocs_per_machine': 1,  # increase this to 4 in order to speed up the calculation
    }
    builder.cp2k_base.cp2k.metadata.options.withmpi = False  # comment this for parallel cp2k executable
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 60 * 60

    # The following is not needed, if the files are available in the data directory of your CP2K executable
    cp2k_dir = DATA_DIR / 'cp2k'
    builder.cp2k_base.cp2k.file = {
        'basis': SinglefileData(file=str(cp2k_dir / 'BASIS_MOLOPT')),
        'pseudo': SinglefileData(file=str(cp2k_dir / 'GTH_POTENTIALS')),
        'dftd3': SinglefileData(file=str(cp2k_dir / 'dftd3.dat')),
    }

    node: WorkChainNode
    results, node = engine.run_get_node(builder)

    print(node.attributes)

    for called in node.called_descendants:
        print()
        print(called.attributes)
        if isinstance(called, CalcJobNode):
            wdir = Path(called.get_remote_workdir())
            print('Remote work directory:', wdir)
            if wdir.exists():
                print([str(p.relative_to(wdir)) for p in wdir.glob('**/*')])
                if wdir.joinpath('aiida.coords.xyz').exists():
                    print("aiida.coords.xyz:")
                    print(wdir.joinpath('aiida.coords.xyz').read_text())
            else:
                print('Remote work directory does not exist')

    raise

    assert node.is_finished_ok, results

    params = results['output_parameters'].get_dict()
    assert params['binding_energy_raw'] == pytest.approx(-24.86, abs=0.1)
    assert params['binding_energy_corr'] == pytest.approx(-23.07, abs=0.1)
    assert params['motion_step_info']['scf_converged'][-1]


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
def cli(cp2k_code):
    """Run example.

    Example usage: $ ./test_Cp2kBindingEnergy_CO2_MOF74.py --cp2k-code my-cp2k@myhost
    """
    run_binding_energy_co2_mof74(cp2k_code,
                                 zn_mof74=StructureData(ase=ase.io.read(DATA_DIR / 'Zn-MOF-74.cif')),
                                 co2_in_mof74=StructureData(ase=ase.io.read(DATA_DIR / 'CO2_in_Zn-MOF-74.cif')))


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
