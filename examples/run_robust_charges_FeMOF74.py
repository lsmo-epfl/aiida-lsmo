from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.structure import StructureData  # noqa
from aiida.orm.data.parameter import ParameterData  # noqa
from aiida.work.run import submit

from ase.io import read
from aiida_lsmo_workflows.geoopt_charges import Cp2kGeoOptDdecWorkChain

atoms = read('Fe-MOF-74.cif')
structure = StructureData(ase=atoms)
structure.label='Fe-MOF-74'
structure.store()

cp2k_options = {
    "resources": {
        "num_machines": 2,
    },
    "max_wallclock_seconds": 1 * 60 * 60,
    }

params_dict = {
        'MOTION':{
            'MD':{
                'STEPS': 5,
                },
            'GEO_OPT': {
                'MAX_ITER': 5,
            },
            'CELL_OPT': {
                'MAX_ITER': 5,
            },
        },
}
ddec_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 1 * 60 * 60 / 2,
    "withmpi": False,
    }
cp2k_parameters = ParameterData(dict=params_dict)
cp2k_code = test_and_get_code('cp2k@fidis-debug', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis-debug', expected_code_type='ddec')
submit(Cp2kGeoOptDdecWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        cp2k_parameters=cp2k_parameters,
        _cp2k_options=cp2k_options,
        ddec_code=ddec_code,
        _ddec_options=ddec_options,
        _label='MyFirstWokchain',
        _guess_multiplicity=True,
        )
