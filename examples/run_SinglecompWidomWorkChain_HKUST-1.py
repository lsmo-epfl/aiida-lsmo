#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run example Widom calculation in HKUST1 framework, using 0.7 scaling for the probe radius to compute blocks."""
import os
import click

from aiida.engine import run
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Str, Dict
from aiida import cmdline

# Workchain objects
SinglecompWidomWorkChain = WorkflowFactory('lsmo.singlecomp_widom')

# Data objects
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')


@click.command('cli')
@cmdline.utils.decorators.with_dbenv()
@click.argument('raspa_code_label')
@click.argument('zeopp_code_label')
def main(raspa_code_label, zeopp_code_label):
    """Prepare inputs and submit the workchain.
    Usage: verdi run run_thisworkchainexample.py raspa@localhost zeopp@localhost"""

    builder = SinglecompWidomWorkChain.get_builder()

    builder.metadata.label = 'test'

    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.zeopp.code = Code.get_from_string(zeopp_code_label)

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
    builder.structure = CifData(file=os.path.abspath('data/HKUST-1.cif'), label='HKUST-1')
    builder.molecule = Str('h2o')

    builder.parameters = Dict(
        dict={
            'zeopp_probe_scaling': 0.7,
            'zeopp_block_samples': 10,  # Default: 100
            'raspa_widom_cycles': 100,  # Default: 1e5
            'temperatures': [200, 300]
        })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
