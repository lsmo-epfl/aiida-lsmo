# -*- coding: utf-8 -*-
"""IsothermMultiTemp workchain."""
from __future__ import absolute_import

import os

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Float, Int, Str, List, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, append_, while_, if_
from aiida_lsmo.utils import check_resize_unit_cell, aiida_dict_merge

# Workchain objects
IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  # pylint: disable=invalid-name

def get_parameters_singletemp(i, parameters):
    parameters_singletemp = parameters.get_dict()
    parameters_singletemp['temperature'] = parameters_singletemp['temperature_list'][i]
    parameters_singletemp['temperature_list'] = None
    return Dict(dict=parameters_singletemp)

@calcfunction
def get_isotherms_output(parameters, **isotherm_dict):
    """Gather together all the results, returning lists for the multi temperature values"""

    multi_temp_labels = [
        'temperature',
        'henry_coefficient_average',
        'henry_coefficient_dev',
        'adsorption_energy_widom_average',
        'adsorption_energy_widom_dev',
        'is_kh_enough',
        'isotherm'
    ]

    single_temp_labels = [
        'temperature_unit',
        'henry_coefficient_unit',
        'adsorption_energy_widom_unit',
        'conversion_factor_molec_uc_to_cm3stp_cm3',
        'conversion_factor_molec_uc_to_gr_gr',
        'conversion_factor_molec_uc_to_mol_kg'
    ]

    isotherms_output = {}
    for label in multi_temp_labels:
        isotherms_output[label] = []

    for i in range(len(isotherm_dict)):
        for label in multi_temp_labels:
            isotherm_out_i = isotherm_dict['isotherm_out_{}'.format(i)]
            isotherms_output[label].append(isotherm_out_i[label])

    # Same for all, take the last for convenience
        for label in single_temp_labels:
            isotherms_output[label] = isotherm_out_i[label]

    return Dict(dict=isotherms_output)

class IsothermMultiTempWorkChain(WorkChain):
    """ Run IsothermWorkChain for multiple temperatures: first compute geometric properties
    and then submit Widom+GCMC at different temperatures in parallel
    """

    @classmethod
    def define(cls, spec):
        super(IsothermMultiTempWorkChain, cls).define(spec)

        spec.expose_inputs(IsothermWorkChain)

        spec.outline(
            cls.run_geometric,
            if_(cls.should_continue)(  # if porous
                cls.run_isotherms,  # run raspa widom calculation
                ),
                cls.collect_isotherms
            )

        spec.expose_outputs(IsothermWorkChain, include=['geometric_output','blocking_spheres'])

        spec.outputs.dynamic = True  # any outputs are accepted

    def run_geometric(self):
        """Perform Zeo++ block and VOLPO calculation with IsothermWC."""

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(IsothermWorkChain)

        # Set inputs for zeopp
        inputs.update({
            'metadata': {
                'label': "IsothermGeometric",
                'call_link_label': 'run_geometric',
            },
        })

        running = self.submit(IsothermWorkChain, **inputs)
        self.report("Computing common gemetric properties")
        return ToContext(geometric=running)

    def should_continue(self):
        """Continue if porous"""

        self.out('geometric_output', self.ctx.geometric.outputs['geometric_output']) #TODO: use expose instead
        # Why not working?:
        #self.out_many(self.exposed_outputs(self.ctx.geometric, IsothermWorkChain))
        #TODO: expose also blocking spheres if exposed from IsothermWC
        return self.outputs['geometric_output']['is_porous']

    def run_isotherms(self):
        """ Compuite isotherms at different temperatures """

        self.ctx.ntemp = len(self.inputs.parameters['temperature_list'])

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(IsothermWorkChain)
        inputs['geometric'] = self.outputs['geometric_output']

        # Update the parameters with only one temperature and submit
        for i in range(self.ctx.ntemp):
            self.ctx.parameters_singletemp = get_parameters_singletemp(i, self.inputs.parameters)

            inputs.update({
                'metadata': {
                    'label': "Isotherm_{}".format(i),
                    'call_link_label': 'run_isotherm_{}'.format(i),

                },
                'parameters': self.ctx.parameters_singletemp
            })

            running = self.submit(IsothermWorkChain, **inputs)
            self.to_context(**{'isotherm_{}'.format(i):running})

    def collect_isotherms(self):
        """ Collect all the results in one Dict """

        output_dict = {}
        for i in range(self.ctx.ntemp):
            output_dict['isotherm_out_{}'.format(i)] =  self.ctx['isotherm_{}'.format(i)].outputs['isotherm_output']

        self.out("isotherms_output", get_isotherms_output(self.inputs.parameters,**output_dict))

        self.report("All the isotherms computed: geom Dict<{}>, isotherms Dict<{}>".format(
            self.outputs['geometric_output'].pk, self.outputs['isotherms_output'].pk))
