#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Submit ZeoppMultistageDdecWorkChain for COF-5"""

import os

from aiida.engine import submit
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Code, Dict, Str

# Workchain objects
ZeoppMultistageDdecWorkChain = WorkflowFactory('lsmo.zeopp_multistage_ddec')  # pylint: disable=invalid-name

#Data objects
CifData = DataFactory('core.cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name

print('NOTE: this test will perform a real-case calculation on Fidis, taking a couple of hours!')

cp2k_code = Code.get_from_string('cp2k-5.1@fidis')
ddec_code = Code.get_from_string('ddec@fidis')
zeopp_code = Code.get_from_string('zeopp@fidis')

cp2k_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 72 * 60 * 60,
    'withmpi': True,
}

ddec_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 10 * 60 * 60,
    'withmpi': False,
}

zeopp_options = {
    'resources': {
        'num_machines': 1
    },
    'max_wallclock_seconds': 10 * 60 * 60,
    'withmpi': False,
}

ddec_params = Dict(
    {
        'net charge':
            0.0,
        'charge type':
            'DDEC6',
        'periodicity along A, B, and C vectors': [True, True, True],
        'compute BOs':
            False,
        'atomic densities directory complete path':
            '/home/ongari/aiida-database/data/chargemol_09_26_2017/atomic_densities/',
        'input filename':
            'valence_density',
    })

structure = CifData(file=os.path.join(os.getcwd(), 'data/COF-5.cif')).store()
structure.label = 'COF-5'

inputs = {
    'structure': structure,
    'protocol_tag': Str('standard'),
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
        #'parameters': zeopp_params, #using standard
        'code': zeopp_code,
        'metadata': {
            'options': zeopp_options,
        }
    }
}
wc = submit(ZeoppMultistageDdecWorkChain, **inputs)

print('Submitted CifData<{}> to ZeoppMultistageDdecWorkChain<{}>'.format(structure.pk, wc.pk))
