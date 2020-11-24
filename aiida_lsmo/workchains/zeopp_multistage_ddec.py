# -*- coding: utf-8 -*-
"""ZeoppMultistageDdecWorkChain work chain"""
import functools

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext
from aiida_lsmo.utils import get_structure_from_cif, validate_dict

from .parameters_schemas import NUMBER, Required, Schema

# import sub-workchains
Cp2kMultistageDdecWorkChain = WorkflowFactory('lsmo.cp2k_multistage_ddec')  # pylint: disable=invalid-name

# import calculations
DdecCalculation = CalculationFactory('ddec')  # pylint: disable=invalid-name
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  # pylint: disable=invalid-name
NetworkParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


class ZeoppMultistageDdecWorkChain(WorkChain):
    """A workchain that combines: Zeopp + Cp2kMultistageWorkChain + Cp2kDdecWorkChain + Zeopp"""

    parameters_schema = Schema({
        Required('ha', default='DEF', description='Using high accuracy (mandatory!)'):
            str,
        Required('res', default=True, description='Max included, free and incl in free sphere'):
            bool,
        Required('sa', default=[1.86, 1.86, 100000], description='Nitrogen probe to compute surface'): [
            NUMBER, NUMBER, int
        ],
        Required('vol', default=[0.0, 0.0, 1000000], description='Geometric pore volume'): [NUMBER, NUMBER, int],
        Required('volpo', default=[1.86, 1.86, 100000], description='Nitrogen probe to compute PO pore volume'): [
            NUMBER, NUMBER, int
        ],
        Required('psd', default=[1.2, 1.2, 10000], description='Small probe to compute the pore size distr'): [
            NUMBER, NUMBER, int
        ],
    })
    parameters_info = parameters_schema.schema  # shorthand for printing

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super().define(spec)

        spec.input('structure', valid_type=CifData, help='input structure')
        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', exclude=['structure'])
        spec.inputs['zeopp']['parameters'].default = lambda: NetworkParameters(dict=cls.parameters_schema({}))
        spec.inputs['zeopp']['parameters'].validator = functools.partial(validate_dict, schema=cls.parameters_schema)
        spec.inputs['zeopp']['parameters'].required = False

        spec.expose_inputs(Cp2kMultistageDdecWorkChain, exclude=['structure'])

        spec.outline(cls.run_zeopp_before, cls.run_multistageddec, cls.run_zeopp_after, cls.return_results)

        spec.expose_outputs(ZeoppCalculation, namespace='zeopp_before_opt', include=['output_parameters'])
        spec.expose_outputs(ZeoppCalculation, namespace='zeopp_after_opt', include=['output_parameters'])
        spec.expose_outputs(Cp2kMultistageDdecWorkChain)

    def run_zeopp_before(self):
        """Run Zeo++ for the original structure"""
        #Merging all inputs
        zeopp_inp = AttributeDict(self.exposed_inputs(ZeoppCalculation, 'zeopp'))
        zeopp_inp['parameters'] = self.inputs.zeopp.parameters
        zeopp_inp['structure'] = self.inputs.structure
        zeopp_inp['metadata']['label'] = 'zeopp_before_opt'
        zeopp_inp['metadata']['call_link_label'] = 'call_zeopp_before_opt'

        running_zeopp = self.submit(ZeoppCalculation, **zeopp_inp)
        self.report('Running Zeo++ calculation <{}>'.format(running_zeopp.pk))
        return ToContext(zeopp_before=running_zeopp)

    def run_multistageddec(self):
        """Run MultistageDdec work chain"""
        msddec_inputs = AttributeDict(self.exposed_inputs(Cp2kMultistageDdecWorkChain))
        msddec_inputs['structure'] = get_structure_from_cif(self.inputs.structure)
        msddec_inputs['metadata']['call_link_label'] = 'call_multistageddec'

        running = self.submit(Cp2kMultistageDdecWorkChain, **msddec_inputs)
        return ToContext(msddec_wc=running)

    def run_zeopp_after(self):
        """Run Zeo++ for the oprimized structure"""
        zeopp_inp = AttributeDict(self.exposed_inputs(ZeoppCalculation, 'zeopp'))
        zeopp_inp['parameters'] = self.inputs.zeopp.parameters
        zeopp_inp['structure'] = self.ctx.msddec_wc.outputs.structure_ddec
        zeopp_inp['metadata']['label'] = 'zeopp_after_opt'
        zeopp_inp['metadata']['call_link_label'] = 'call_zeopp_after_opt'

        running_zeopp = self.submit(ZeoppCalculation, **zeopp_inp)
        self.report('Running Zeo++ calculation <{}>'.format(running_zeopp.pk))
        return ToContext(zeopp_after=running_zeopp)

    def return_results(self):
        """Return exposed outputs"""
        self.out_many(self.exposed_outputs(self.ctx.zeopp_before, ZeoppCalculation, namespace='zeopp_before_opt'))
        self.out_many(self.exposed_outputs(self.ctx.msddec_wc, Cp2kMultistageDdecWorkChain))
        self.out_many(self.exposed_outputs(self.ctx.zeopp_after, ZeoppCalculation, namespace='zeopp_after_opt'))
        self.report('WorkChain terminated correctly.')
