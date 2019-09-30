from __future__ import print_function
from __future__ import absolute_import

import sys, os
from ase.io import read

from aiida.engine import submit
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Float, Str

# Workchain objects
ZeoppMultistageDdecWorkChain = WorkflowFactory('lsmo.zeoppmultistageddec')  # pylint: disable=invalid-name

#Data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory("zeopp.parameters")  # pylint: disable=invalid-name


cp2k_code = Code.get_from_string('cp2k@localhost')
ddec_code = Code.get_from_string('ddec@localhost')
zeopp_code = Code.get_from_string('zeopp@localhost')

cp2k_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 10 * 60,
    'withmpi': True,
}

ddec_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 10 * 60,
    'withmpi': False,
}

zeopp_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 10 * 60,
    'withmpi': False,
}

ddec_params = Dict(
    dict={
        'net charge': 0.0,
        'charge type': 'DDEC6',
        'periodicity along A, B, and C vectors': [True, True, True],
        'compute BOs': False,
        'atomic densities directory complete path':
            '/home/daniele/Programs/aiida-database/data/chargemol_09_26_2017/atomic_densities/',
        'input filename': 'valence_density',
    })

zeopp_params = NetworkParameters(dict={
    'ha': 'DEF',                    # Using high accuracy (mandatory!)
    'res': True,                   # Max included, free and incl in free sphere
    'sa': [1.86, 1.86, 1000],    # Nitrogen probe to compute surface
    'vol': [0.0, 0.0, 1000],    # Geometric pore volume
})

structure = CifData(file=os.path.join(os.getcwd(), 'data/H2O.cif')).store()
structure.label = 'H2O'

inputs = {
    'structure': structure,
    'protocol_tag': Str('test'),
    'metadata': {
        'label': 'test',
    },
    'cp2k_base': {
        'cp2k': {
            'code': cp2k_code,
            'metadata': {
                'options': cp2k_options,
            }
        }
    },
    'ddec': {
        'parameters': ddec_params,
        'code': ddec_code,
        'metadata': {
            'options': ddec_options,
        }
    },
    'zeopp': {
        'parameters': zeopp_params,
        'code': zeopp_code,
        'metadata': {
            'options': zeopp_options,
        }
    }
}
wc = submit(ZeoppMultistageDdecWorkChain, **inputs)

print("Submitted CifData<{}> to ZeoppMultistageDdecWorkChain<{}>".format(
      structure.pk,
      wc.pk
      ))
