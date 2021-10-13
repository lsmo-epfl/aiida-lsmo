#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example IsothermInflection for Ar in Graphite."""

from pathlib import Path
import os
import click
import pytest

from aiida import engine
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict
from aiida import cmdline

THIS_DIR = Path(__file__).resolve().parent
DATA_DIR = THIS_DIR / 'data'

# Workchain objects
IsothermInflectionWorkChain = WorkflowFactory('lsmo.isotherm_inflection')

# Data objects
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')


@pytest.fixture(scope='function')
def graphite_20a():
    """CifData for graphite."""
    with open(os.path.join(DATA_DIR, 'Graphite_20A.cif'), 'rb') as handle:
        cif = CifData(file=handle, label='Graphite_20A')

    return cif


def run_isotherm_inflection_ar_graphite(raspa_code, zeopp_code, graphite_20a):  # pylint: disable=redefined-outer-name
    """Prepare inputs and submit the workchain.
    Usage: verdi run run_thisworkchainexample.py raspa@localhost zeopp@localhost"""

    builder = IsothermInflectionWorkChain.get_builder()

    builder.metadata.label = 'test'

    builder.raspa_base.raspa.code = raspa_code
    builder.zeopp.code = zeopp_code

    options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
        'withmpi': False,
    }
    builder.raspa_base.raspa.metadata.options = options
    builder.zeopp.metadata.options = options
    builder.structure = graphite_20a

    builder.molecule = Dict(
        dict={
            'name': 'Ar',
            'forcefield': 'HIRSCHFELDER',
            'ff_cutoff': 8,
            'molsatdens': 35.4,
            'proberad': 1.7,
            'singlebead': True,
            'charged': False,
            'pressure_zero': 1,
        })

    builder.parameters = Dict(
        dict={
            'ff_framework': 'DREIDING',
            'temperature': 87,  # T_sat Ar
            'ff_cutoff': 8.0,  # NOTE: Low to have cheap testing
            'box_length': 16.0,
            'zeopp_probe_scaling': 1.0,
            'zeopp_volpo_samples': 10000,
            'zeopp_block_samples': 100,
            'raspa_widom_cycles': 1000,
            'raspa_gcmc_init_cycles': 300,
            'raspa_gcmc_prod_cycles': 300,
            'pressure_num': 4,
            'raspa_verbosity': 10
        })

    results, node = engine.run_get_node(builder)

    assert node.is_finished_ok, results

    params = results['output_parameters'].get_dict()
    assert 'loading_absolute_dev_from_dil' in params['isotherm']


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa-code', type=cmdline.params.types.CodeParamType())
@click.option('--zeopp-code', type=cmdline.params.types.CodeParamType())
def cli(raspa_code, zeopp_code):
    """Run example.

    Example usage: $ ./test_isotherm_inflection_workchain_graphite_ar.py --raspa-code ... --zeopp-code ...

    Help: $ ./test_isotherm_inflection_workchain_graphite.py --help
    """
    with open(os.path.join(DATA_DIR, 'Graphite_20A.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    run_isotherm_inflection_ar_graphite(raspa_code, zeopp_code, cif)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
