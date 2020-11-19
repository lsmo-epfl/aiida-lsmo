#!/usr/bin/env python  # pylint: disable=invalid-name
# -*- coding: utf-8 -*-
"""Run example sim annealing of CO2 in Zn-MOF-74 framework."""

import os
import click
import pytest

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida import cmdline

from . import DATA_DIR

# Workchain objects
SimAnnealingWorkChain = WorkflowFactory('lsmo.sim_annealing')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


@pytest.fixture(scope='function')
def zn_mof74_cifdata():
    """CifData for Zn MOF74 CIF."""
    with open(os.path.join(DATA_DIR, 'Zn-MOF-74.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    return cif


def run_sim_annealing_zn_mof74(raspa_code, zn_mof74_cifdata):  # pylint: disable=redefined-outer-name
    """Prepare inputs and submit the work chain."""

    builder = SimAnnealingWorkChain.get_builder()
    builder.metadata.label = 'test'
    builder.raspa_base.raspa.code = raspa_code
    builder.raspa_base.raspa.metadata.options = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 1 * 60 * 60,
    }
    builder.structure = zn_mof74_cifdata
    builder.molecule = Str('co2')
    builder.parameters = Dict(dict={
        'ff_framework': 'UFF',  # (str) Forcefield of the structure.
        'mc_steps': int(10),  # (int) Number of MC cycles.
    })

    results = run(builder)
    params = results['output_parameters'].get_dict()
    assert params['energy_host/ads_tot_final'][-1] == pytest.approx(-36, abs=5)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa-code', type=cmdline.params.types.CodeParamType())
def cli(raspa_code):
    """Run example.

    Example usage: $ ./test_SimAnnealingWorkChain_MOF74_CO2.py --raspa-code raspa-4467e14@fidis

    Help: $ ./test_SimAnnealingWorkChain_MOF74_CO2.py --help
    """
    with open(os.path.join(DATA_DIR, 'Zn-MOF-74.cif.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    run_sim_annealing_zn_mof74(raspa_code, cif)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
