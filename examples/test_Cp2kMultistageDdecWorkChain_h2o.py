#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Test/example for the DdecMultistageDdec WorkChain."""

import click
import ase.build

from aiida.plugins import DataFactory, WorkflowFactory
from aiida import cmdline
from aiida import engine
from aiida.orm import Dict, StructureData, Str, SinglefileData
from . import DATA_DIR

# Workchain object
MultistageDdecWorkChain = WorkflowFactory('lsmo.cp2k_multistage_ddec')  # pylint: disable=invalid-name

# Data objects
StructureData = DataFactory('structure')  # pylint: disable=invalid-name


def run_cp2k_multistage_ddec_h2o(cp2k_code, ddec_code):  # pylint: disable=redefined-outer-name
    """Run DdecMultistageDdec WorkChainExample on H2O.
    """
    print('Testing CP2K-Multistage calculation + DDEC on H2O...')

    builder = MultistageDdecWorkChain.get_builder()
    atoms = ase.build.molecule('H2O')
    atoms.center(vacuum=2.0)
    builder.structure = StructureData(ase=atoms)
    builder.protocol_tag = Str('test')
    builder.metadata.label = 'test-h2o'

    builder.cp2k_base.cp2k.code = cp2k_code
    builder.cp2k_base.cp2k.metadata.options.resources = {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 1,
    }
    builder.cp2k_base.cp2k.metadata.options.max_wallclock_seconds = 10 * 60
    builder.cp2k_base.cp2k.metadata.options.withmpi = False  # set to True for production
    # The following is not needed, if the files are available in the data directory of your CP2K executable
    cp2k_dir = DATA_DIR / 'cp2k'
    builder.cp2k_base.cp2k.file = {
        'basis': SinglefileData(file=str(cp2k_dir / 'BASIS_MOLOPT')),
        'pseudo': SinglefileData(file=str(cp2k_dir / 'GTH_POTENTIALS')),
        'dftd3': SinglefileData(file=str(cp2k_dir / 'dftd3.dat')),
    }

    builder.ddec.code = ddec_code
    builder.ddec.parameters = Dict(
        dict={
            'net charge': 0.0,
            'charge type': 'DDEC6',
            'periodicity along A, B, and C vectors': [True, True, True],
            'compute BOs': False,
            'input filename': 'valence_density',
        })
    builder.ddec.metadata.options.max_wallclock_seconds = 10 * 60

    results, node = engine.run_get_node(builder)

    assert node.is_finished_ok, results
    assert 'structure_ddec' in results, results
    parameters = results['output_parameters'].get_dict()
    assert 'final_bandgap_spin1_au' in parameters, parameters


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--ddec-code', type=cmdline.params.types.CodeParamType())
@click.option('--cp2k-code', type=cmdline.params.types.CodeParamType())
def cli(ddec_code, cp2k_code):
    """Run DdecMultistageDdec WorkChainExample on H2O.
    """
    run_cp2k_multistage_ddec_h2o(ddec_code=ddec_code, cp2k_code=cp2k_code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
