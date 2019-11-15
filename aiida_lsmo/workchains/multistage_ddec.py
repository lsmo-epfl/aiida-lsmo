# -*- coding: utf-8 -*-
"""MultistageDdecWorkChain workchain"""

from __future__ import absolute_import

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext

# import sub-workchains
Cp2kMultistageWorkChain = WorkflowFactory('cp2k.multistage')  # pylint: disable=invalid-name
Cp2kDdecWorkChain = WorkflowFactory('ddec.cp2k_ddec')  # pylint: disable=invalid-name

# import calculations
DdecCalculation = CalculationFactory('ddec')  # pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  # pylint: disable=invalid-name


class MultistageDdecWorkChain(WorkChain):
    """A workchain that combines: Cp2kMultistageWorkChain + Cp2kDdecWorkChain"""

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super(MultistageDdecWorkChain, cls).define(spec)

        spec.expose_inputs(Cp2kMultistageWorkChain)
        spec.expose_inputs(Cp2kDdecWorkChain, exclude=['cp2k_base'])

        # specify the chain of calculations to be performed
        spec.outline(cls.run_cp2kmultistage, cls.run_cp2kddec, cls.return_results)

        spec.expose_outputs(Cp2kMultistageWorkChain, exclude=['output_structure'])
        spec.expose_outputs(Cp2kDdecWorkChain, include=['structure_ddec'])

    def run_cp2kmultistage(self):
        """Run CP2K-Multistage"""
        cp2k_ms_inputs = AttributeDict(self.exposed_inputs(Cp2kMultistageWorkChain))
        cp2k_ms_inputs['metadata']['call_link_label'] = 'call_cp2kmultistage'
        running = self.submit(Cp2kMultistageWorkChain, **cp2k_ms_inputs)
        self.report('Running Cp2MultistageWorkChain to move the structure')
        return ToContext(ms_wc=running)

    def run_cp2kddec(self):
        """Pass the Cp2kMultistageWorkChain outputs as inputs for
        Cp2kDdecWorkChain: cp2k_base (metadata), cp2k_params, structure and WFN.
        """
        cp2k_ddec_inputs = AttributeDict(self.exposed_inputs(Cp2kDdecWorkChain))
        cp2k_ddec_inputs['cp2k_base'] = self.exposed_inputs(Cp2kMultistageWorkChain)['cp2k_base']
        cp2k_ddec_inputs['cp2k_base']['cp2k']['parameters'] = self.ctx.ms_wc.outputs.last_input_parameters
        cp2k_ddec_inputs['cp2k_base']['cp2k']['structure'] = self.ctx.ms_wc.outputs.output_structure
        cp2k_ddec_inputs['cp2k_base']['cp2k']['parent_calc_folder'] = self.ctx.ms_wc.outputs.remote_folder
        cp2k_ddec_inputs['metadata']['call_link_label'] = 'call_cp2kddec'

        running = self.submit(Cp2kDdecWorkChain, **cp2k_ddec_inputs)
        return ToContext(cp2k_ddec_wc=running)

    def return_results(self):
        """Return exposed outputs and print the pk of the CifData w/DDEC"""
        self.out_many(self.exposed_outputs(self.ctx.ms_wc, Cp2kMultistageWorkChain))
        self.out_many(self.exposed_outputs(self.ctx.cp2k_ddec_wc, Cp2kDdecWorkChain))
