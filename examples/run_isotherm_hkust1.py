# -*- coding: utf-8 -*-
"""Run example isotherm calculation with HKUST1 framework."""
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import click

from aiida.engine import run
from aiida.common import NotExistent
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Float, Int, SinglefileData

# Workchain object
IsothermWorkChain = WorkflowFactory('aiida_lsmo.isotherm')  # pylint: disable=invalid-name

# Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


@click.command('cli')
@click.argument('label_raspa_code')
@click.argument('label_zeopp_code')
def main(label_raspa_code, label_zeopp_code):
    """Prepare inputs and submit the Isotherm workchain."""
    # Test the codes and specify the nodes and walltime
    try:
        raspa_code = Code.get_from_string(label_raspa_code)
    except NotExistent:
        print("The code '{}' does not exist".format(label_raspa_code))
        sys.exit(1)

    try:
        zeopp_code = Code.get_from_string(label_zeopp_code)
    except NotExistent:
        print("The code '{}' does not exist".format(label_zeopp_code))
        sys.exit(1)

    builder = IsothermWorkChain.get_builder()

    builder.metadata.label = "Volpo-Kh-Isotherm"

    # Import and attach the structure
    builder.structure = CifData(file=os.path.abspath("HKUST1.cif"),
                                label="HKUST1")

    builder.raspa_molsatdens = Float(2.5e4)

    # Raspa parameters
    builder.raspa_base.raspa.parameters = Dict(
        dict={
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfCycles": 1000,
                "PrintPropertiesEvery": 100,  # info on henry coeff
                "NumberOfInitializationCycles": 2000,
                "PrintEvery": 200,
                "Forcefield": "GenericMOFs",
                "RemoveAtomNumberCodeFromLabel": True,
                "EwaldPrecision": 1e-6,
                "CutOff": 12.0,
            },
            "System": {
                "hkust1": {
                    "type": "Framework",
                    "UnitCells": "1 1 1",
                    "ExternalTemperature": 300.0,
                    "ExternalPressure": 0,  # Will be set by the workchain
                }
            },
            "Component": {
                "CO2": {
                    "MoleculeDefinition": "TraPPE",
                },
            },
        })

    # Raspa code
    builder.raspa_base.raspa.code = raspa_code

    # Raspa options
    builder.raspa_base.raspa.metadata.options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 60 * 60,
        "withmpi": False,
    }

    # Radius file for the framework
    builder.zeopp.atomic_radii = SinglefileData(
        file=os.path.abspath("UFF.rad"))

    # Zeopp parameters
    builder.zeopp.parameters = NetworkParameters(dict={
        #'ha': 'DEF',
        'volpo': [1.82, 1.82, 1000],
        'block': [1.82, 100],
    })

    # Zeopp code
    builder.zeopp.code = zeopp_code

    builder.raspa_widom_cycle_mult = Int(2)

    # Zeopp options
    builder.zeopp.metadata.options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
        },
        "max_wallclock_seconds": 1 * 60 * 60,
        "withmpi": False,
    }

    run(builder)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

# EOF
