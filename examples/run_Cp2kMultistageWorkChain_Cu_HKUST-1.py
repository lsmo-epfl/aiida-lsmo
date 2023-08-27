#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run Cp2kMultistageWorkChain workchain on Cu-HKUST-1."""

import pytest
import click
import ase.build
import data

from importlib.resources import files

from aiida.plugins import DataFactory, WorkflowFactory
from aiida import cmdline
from aiida import engine
from aiida.orm import Dict, StructureData, Str, SinglefileData


# Workchain objects
Cp2kMultistageWorkChain = WorkflowFactory('lsmo.cp2k_multistage')

# Data objects
StructureData = DataFactory('core.structure')


@pytest.fixture(scope='function')
def cu_hkust1_structuredata():
    """StructureData for Aluminum."""
    return StructureData(ase=ase.io.read(files(data) / 'Cu-I-II-HKUST-1.cif'))


def run_multistage_cu_hkust1(cp2k_code, cu_hkust1_structuredata):  # pylint: disable=redefined-outer-name
    """Run Cp2kMultistageWorkChain on Cu-HKUST-1."""

    # testing user change of parameters and protocol
    parameters = Dict({'FORCE_EVAL': {'DFT': {'MGRID': {'CUTOFF': 250,}}}})

    # Construct process builder
    builder = Cp2kMultistageWorkChain.get_builder()
    builder.structure = cu_hkust1_structuredata
    builder.protocol_tag = Str('standard')
    builder.cp2k_base.cp2k.parameters = parameters
    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        'num_machines': 1,
        #'num_mpiprocs_per_machine': 1,
    }
    #builder.cp2k_base.cp2k.metadata.options.withmpi = False  # comment this for parallel cp2k executable
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 1 * 60 * 60
    # The following is not needed, if the files are available in the data directory of your CP2K executable
    cp2k_dir = files(data) / 'cp2k'
    builder.cp2k_base.cp2k.file = {
        'basis': SinglefileData(file=str(cp2k_dir / 'BASIS_MOLOPT')),
        'pseudo': SinglefileData(file=str(cp2k_dir / 'GTH_POTENTIALS')),
        'dftd3': SinglefileData(file=str(cp2k_dir / 'dftd3.dat')),
    }
    results, node = engine.run_get_node(builder)

    import pdb
    pdb.set_trace()

    assert node.is_finished_ok, results
    output_parameters = results['output_parameters'].get_dict()
    assert output_parameters['step_info']['scf_converged'][-1]


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
def cli(cp2k_code):
    """Run example.

    Example usage: $ ./test_multistage_aluminum.py --cp2k-code my-cp2k@myhost
    """
    run_multistage_cu_hkust1(cp2k_code,
                             cu_hkust1_structuredata=StructureData(ase=ase.io.read(files(data) / 'Cu-I-II-HKUST-1.cif')))


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
