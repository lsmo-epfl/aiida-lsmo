import os
import numpy as np
from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.data.base import Float
from aiida.work.run import submit

from ase.io import read
from aiida_lsmo_workflows.isotherm import Isotherm

# data objects
ArrayData = DataFactory('array')
ParameterData = DataFactory('parameter')
CifData = DataFactory('cif')

#structure = load_node(6952)
structure = CifData(file=os.getcwd()+'/Cu-MOF-74.cif')

# option for zeo++ and raspa
zr_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 1 * 60 * 60,
    "withmpi": False,
    }

raspa_parameters = ParameterData(dict={
        "GeneralSettings":
        {
        "SimulationType"                   : "MonteCarlo",
        #"NumberOfCycles"                   : 100000,
        "NumberOfCycles"                   : 1000,
        #"NumberOfInitializationCycles"     : 10000,   # 20000
        "NumberOfInitializationCycles"     : 2000,   # 20000

        #"PrintEvery"                       : 10000, "PrintEvery"                       : 1000,

        "ChargeMethod"                     : "Ewald",
        "CutOff"                           : 12.0,
        "Forcefield"                       : "LSMO_UFF-TraPPE",
        "EwaldPrecision"                   : 1e-6,

        "Framework"                        : 0,
        "UnitCells"                        : "1 1 1",
        "HeliumVoidFraction"               : 0.0,

        "ExternalTemperature"              : 298.0,
        "ExternalPressure"                 : 58e4,
        },
        "Component":
        [{
        "MoleculeName"                     : "CO2",
        "MoleculeDefinition"               : "TraPPE",
        "TranslationProbability"           : 0.5,
        "RotationProbability"              : 0.5,
        "ReinsertionProbability"           : 0.5,
        "SwapProbability"                  : 1.0,
        "CreateNumberOfMolecules"          : 0,
        }],
        })

zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa@deneb', expected_code_type='raspa')


pressures = ArrayData()
pressures.set_array("pressures", np.array([0.01e5]))
submit(Isotherm,
        structure=structure,
        probe_radius=Float(1.525),
        pressures=pressures,
        min_cell_size=Float(10.0),
        zeopp_code=zeopp_code,
        _zeopp_options=zr_options,
        raspa_code=raspa_code,
        raspa_parameters=raspa_parameters,
        _raspa_options=zr_options,
        _usecharges=True,
        _label='Isotherm',
        )
