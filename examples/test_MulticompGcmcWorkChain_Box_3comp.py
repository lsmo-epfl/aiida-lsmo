#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example 3-components GCMC in a box, at 3 different T/P conditions."""
import click

from aiida import engine
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict
from aiida import cmdline

# Workchain objects
MulticompGcmcWorkChain = WorkflowFactory('lsmo.multicomp_gcmc')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


def run_multicomp_gcmc_box(raspa_code, zeopp_code):  # pylint: disable=redefined-outer-name
    """Prepare inputs and submit the workchain."""

    builder = MulticompGcmcWorkChain.get_builder()

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
    builder.conditions = Dict(dict={
        'molfraction': {
            'co': 0.2,
            'ethene': 0.3,
            'ethane': 0.5,
        },
        'temp_press': [
            [200, 0.1],
            [300, 0.5],
            [400, 0.7],
        ]
    })

    builder.parameters = Dict(dict={
        'raspa_gcmc_init_cycles': 1000,  # Default: 1e3
        'raspa_gcmc_prod_cycles': 1000,  # Default: 1e4
    })

    results, node = engine.run_get_node(builder)

    assert node.is_finished_ok, results

    params = results['output_parameters'].get_dict()
    for molecule in ['CO', 'C2H4', 'C2H6']:
        assert molecule in params['loading_absolute_average']


@click.command()
@cmdline.utils.decorators.with_dbenv()
@click.option('--raspa-code', type=cmdline.params.types.CodeParamType())
@click.option('--zeopp-code', type=cmdline.params.types.CodeParamType())
def cli(raspa_code, zeopp_code):
    """Run example.

    Example usage: $ ./test_MulticompGcmcWorkChain_Box_3comp.py --raspa-code ... --zeopp-code ...

    Help: $ ./test_MulticompGcmcWorkChain_Box_3comp.py --help
    """
    run_multicomp_gcmc_box(raspa_code, zeopp_code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
