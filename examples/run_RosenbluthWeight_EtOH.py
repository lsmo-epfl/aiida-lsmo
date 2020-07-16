#!/usr/bin/env python  # pylint: disable=invalid-name
"""Run example RosenbluthWeight work chain for EtOH."""

import click

from aiida.engine import run
from aiida.plugins import WorkflowFactory
from aiida.orm import Code, Dict

# Workchain objects
RosenbluthWeightWorkChain = WorkflowFactory('lsmo.rosenbluth_weight')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('raspa_code_label')
def main(raspa_code_label):
    """Prepare inputs and submit the work chain."""

    builder = RosenbluthWeightWorkChain.get_builder()
    builder.metadata.label = "test"
    builder.raspa_base.raspa.code = Code.get_from_string(raspa_code_label)
    builder.raspa_base.raspa.metadata.options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 60 * 60,
        "withmpi": False,
    }

    builder.molecule = Dict(dict={
        'name': 'EtOH',
        'forcefield': 'TraPPE',
        'charged': True,
    })

    builder.parameters = Dict(dict={
        'temperature': 300,
        'raspa_cycles': 1000,
    })

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
