import os
import numpy as np
from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.data.base import Float, Str
from aiida.work.run import submit

from ase.io import read
from workflows.isotherm import Isotherm

# data objects
ArrayData = DataFactory('array')
ParameterData = DataFactory('parameter')
CifData = DataFactory('cif')

structure = CifData(file=os.getcwd()+'/Cu-MOF-74.cif')

cp2k_options = {
    "resources": {
        "num_machines": 2,
    },
    "max_wallclock_seconds": 1 * 60 * 60,
    }

ddec_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 1 * 60 * 60 / 2,
    "withmpi": False,
    }

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

        #"PrintEvery"                       : 10000,
        "PrintEvery"                       : 1000,

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

cp2k_code = test_and_get_code('cp2k@fidis-debug', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis-debug', expected_code_type='ddec')
zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa@deneb', expected_code_type='raspa')


pressures = ArrayData()
pressures.set_array("pressures", np.array([0.01e5, 0.05e5, 0.1e5, 0.15e5, 0.2e5]))
submit(Isotherm,
        structure=structure,
        probe_molecule=ParameterData(dict={"sigma":1.525}),
        pressures=pressures,
        min_cell_size=Float(10.0),
        cp2k_code=cp2k_code,
        _cp2k_options=cp2k_options,
        ddec_code=ddec_code,
        _ddec_options=ddec_options,
        zeopp_code=zeopp_code,
        _zeopp_options=zr_options,
        raspa_code=raspa_code,
        raspa_parameters=raspa_parameters,
        _raspa_options=zr_options,
        _usecharges=True,
        _guess_multiplicity=True,
        _label='Isotherm',
        )
