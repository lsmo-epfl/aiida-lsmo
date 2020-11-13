#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example isotherm calculation on Mg MOF 74."""

import os
import click
import pytest

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida import cmdline

from . import DATA_DIR

# Workchain objects
IsothermWorkChain = WorkflowFactory('lsmo.isotherm')

# Data objects
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')


@pytest.fixture(scope='function')
def mg_mof74_cifdata():
    """CifData for Mg MOF74 CIF."""
    with open(os.path.join(DATA_DIR, 'Mg_MOF_74.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    return cif


def run_isotherm_mg_mof74(raspa_code, zeopp_code, mg_mof74_cifdata):  # pylint: disable=redefined-outer-name
    """Test Isotherm workchain on MOF 74."""

    builder = IsothermWorkChain.get_builder()

    builder.zeopp.code = zeopp_code
    builder.zeopp.metadata.options = {'max_wallclock_seconds': 60}

    builder.raspa_base.raspa.code = raspa_code
    builder.raspa_base.raspa.metadata.options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
        'withmpi': False,
    }

    builder.structure = mg_mof74_cifdata
    builder.molecule = Str('co2')
    builder.parameters = Dict(
        dict={
            'ff_framework': 'UFF',  # Default: UFF
            'temperature': 400,  # (K) Note: higher temperature will have less adsorbate and it is faster
            'zeopp_volpo_samples': 1000,  # Default: 1e5 *NOTE: default is good for standard real-case!
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 10,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
            'pressure_min': 0.001,  # Default: 0.001 (bar)
            'pressure_max': 3,  # Default: 10 (bar)
        })

    results = run(builder)

    params = results['output_parameters'].get_dict()

    assert params['Estimated_saturation_loading'] == pytest.approx(13.49433, 0.1), params


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa-code', type=cmdline.params.types.CodeParamType())
@click.option('--zeopp-code', type=cmdline.params.types.CodeParamType())
def cli(raspa_code, zeopp_code, cif):
    """Run example.

    Example usage: $ ./test_isotherm_workchain.py --raspa-code raspa-4467e14@fidis --zeopp-code zeopp-46ce745@fidis

    Help: $ ./test_isotherm_workchain.py --help
    """
    with open(os.path.join(DATA_DIR, 'Mg_MOF_74.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    run_isotherm_mg_mof74(raspa_code, zeopp_code, cif)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
