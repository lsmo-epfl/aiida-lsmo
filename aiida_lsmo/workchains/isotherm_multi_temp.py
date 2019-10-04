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
    parameters_singletemp'temperature_list'] = None
    return Dict(dict=paramterers_singletemp)

@calcfunction
def get_isotherms_output(**dict):
    """Gather together all the results"""
    temperature_list = 


class IsothermMultiTempWorkChain(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient.
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

        spec.expose_outputs(IsothermWorkChain)

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

        self.out_many(self.exposed_outputs(self.ctx.geometric, IsothermWorkChain))

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
                'geometric': self.outputs['geometric_output']['is_porous'],
                'parameters': self.ctx.parameters_singletemp
            })

            running = self.submit(IsothermWorkChain, **inputs)
            assert self.to_context(**{'isotherm_{}'.format(i):running})

    def collect_isotherms(self):
        """ Collect all the results in one Dict """
        output_dict = {}
        for i in range(self.ctx.ntemp):
            output_dict['widom_out_{}'] = self.ctx['isotherm_{}'.format(i)].outputs['widom_output']
            output_dict['isotherm_out_{}'] =  self.ctx['isotherm_{}'.format(i)].outputs['widom_output']

        self.out("isotherms_output", get_isotherms_output(self.inputs.parameters,**output_dict))
