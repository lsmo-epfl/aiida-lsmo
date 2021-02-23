#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example two-component isotherm calculation with HKUST1 framework."""

import os
import click
import pytest

from aiida import engine
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict
from aiida import cmdline

from . import DATA_DIR

# Workchain objects
MulticompAdsDesWorkChain = WorkflowFactory('lsmo.multicomp_ads_des')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name
SinglefileData = DataFactory('singlefile')


@pytest.fixture(scope='function')
def hkust_1_cifdata():
    """CifData for HKUST-1."""
    with open(os.path.join(DATA_DIR, 'HKUST-1.cif'), 'rb') as handle:
        cif = CifData(file=handle, label='HKUST-1')

    return cif


def run_multicomp_ads_des_hkust_1(raspa_code, zeopp_code, hkust_1_cifdata):  # pylint: disable=redefined-outer-name
    """Prepare inputs and submit the Isotherm workchain."""

    builder = MulticompAdsDesWorkChain.get_builder()

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
    builder.structure = hkust_1_cifdata
    builder.conditions = Dict(
        dict={
            'molfraction': {
                'xenon': 0.2,
                'krypton': 0.8
            },
            'adsorption': {
                'temperature': 298,  #K
                'pressure': 1,  #bar
            },
            'desorption': {
                'temperature': 308,
                'pressure': 0.1,
            },
        })

    builder.parameters = Dict(
        dict={
            'zeopp_probe_scaling': 0.9,  # Default: 1.0
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'raspa_gcmc_init_cycles': 100,  # Default: 1e3
            'raspa_gcmc_prod_cycles': 100,  # Default: 1e4
        })

    results, node = engine.run_get_node(builder)

    assert node.is_finished_ok, results

    params = results['output_parameters'].get_dict()

    # checking results
    assert params['loading_absolute_average']['Kr'][1] == pytest.approx(0.02, abs=0.02)
    assert params['loading_absolute_average']['Xe'][1] == pytest.approx(0.74, abs=0.2)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa-code', type=cmdline.params.types.CodeParamType())
@click.option('--zeopp-code', type=cmdline.params.types.CodeParamType())
def cli(raspa_code, zeopp_code):
    """Run example.

    Example usage: $ ./test_MulticompAdsDesWorkChain_HKUST-1.py --raspa-code ... --zeopp-code ...

    Help: $ ./test_MulticompAdsDesWorkChain_HKUST.py --help
    """
    with open(os.path.join(DATA_DIR, 'HKUST-1.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    run_multicomp_ads_des_hkust_1(raspa_code, zeopp_code, cif)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
