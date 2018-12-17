import os
from aiida.common.example_helpers import test_and_get_code
from aiida.orm import DataFactory
from aiida.orm.data.base import Float
from aiida.orm.calculation.work import WorkCalculation
from aiida.work.run import submit
from aiida_lsmo_workflows.volpo_Kh import VolpoKh
ParameterData = DataFactory('parameter')
SinglefileData = DataFactory('singlefile')
CifData = DataFactory('cif')

def dict_merge_ez(dict1, dict2):
    sumdicts = dict1.copy()
    sumdicts.update(dict2)
    return sumdicts

# Test the codes and specify the nodes and walltime
zeopp_code = test_and_get_code('zeopp@localhost', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa@localhost', expected_code_type='raspa')

zeopp_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 3 * 60 * 60,
    "withmpi": False,
    }
raspa_options = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 72 * 60 * 60,
    "withmpi": False,
    }

# Settings for Zeopp (Block-pockets and VolPO) and Raspa (Widom)
zeopp_probe_radius_co2_trappe = Float(2.0) # It will create 8 pore blocks for test purpose
zeopp_atomic_radii_file = SinglefileData(file=os.path.abspath("./UFF.rad")) # Radius file for the framework
raspa_params_dict = {
        "GeneralSettings":
        {
        "NumberOfCycles"                   : 1000,
        "PrintPropertiesEvery"             : 100,  # info on henry coeff
        "Forcefield"                       : "LSMO_UFF-TraPPE",
        "CutOff"                           : 12.0,
        "ExternalTemperature"              : 300.0,
        },
}
raspa_co2_dict = {
        "Component":
        [{
        "MoleculeName"                     : "CO2",
        "MoleculeDefinition"               : "TraPPE",
        "WidomProbability"                 : 1.0,
        }],
}

# Import the structure
structure = CifData(file=os.path.abspath("./HKUST1.cif"))
structure.label = "HKUST1"

# Run for CO2, using UFF-TraPPE force field
submit(VolpoKh,
    structure=structure,
    zeopp_code=zeopp_code,
    _zeopp_options=zeopp_options,
    zeopp_probe_radius=zeopp_probe_radius_co2_trappe,
    zeopp_atomic_radii=zeopp_atomic_radii_file,
    raspa_code=raspa_code,
    raspa_parameters=ParameterData(dict=dict_merge_ez(raspa_params_dict,raspa_co2_dict)),
    _raspa_options=raspa_options,
    _raspa_usecharges=True,
    _label='volpo-Kh-test',
    )
