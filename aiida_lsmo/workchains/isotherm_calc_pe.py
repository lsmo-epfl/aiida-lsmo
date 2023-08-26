# -*- coding: utf-8 -*-
"""IsothermCalcPE work chain."""
import functools

from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida.engine import WorkChain
from aiida_lsmo.utils import dict_merge, validate_dict
from aiida_lsmo.calcfunctions import PE_PARAMETERS_DEFAULT, calc_co2_parasitic_energy

from .parameters_schemas import FF_PARAMETERS_VALIDATOR, NUMBER, Required

# import sub-workchains
IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  #pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('core.cif')  #pylint: disable=invalid-name


class IsothermCalcPEWorkChain(WorkChain):
    """ Compute CO2 parassitic energy (PE) after running IsothermWorkChain for CO2 and N2 at 300K."""

    parameters_schema = FF_PARAMETERS_VALIDATOR.extend({
        Required('zeopp_probe_scaling', default=1.0, description="scaling probe's diameter: molecular_rad * scaling"):
            NUMBER,
        Required('zeopp_volpo_samples', default=int(1e5)):
            int,  # Number of samples for VOLPO calculation (per UC volume).
        Required('zeopp_block_samples',
                 default=int(100),
                 description='Number of samples for BLOCK calculation (per A^3).'):
            int,
        Required('raspa_verbosity', default=10, description='Print stats every: number of cycles / raspa_verbosity.'):
            int,
        Required('raspa_widom_cycles', default=int(1e5), description='Number of Widom cycles.'):
            int,
        Required('raspa_gcmc_init_cycles', default=int(1e3), description='Number of GCMC initialization cycles.'):
            int,
        Required('raspa_gcmc_prod_cycles', default=int(1e4), description='Number of GCMC production cycles.'):
            int,
        Required('raspa_minKh',
                 default=1e-10,
                 description='If Henry coefficient < raspa_minKh do not run the isotherm (mol/kg/Pa).'):
            NUMBER,
        Required('temperature', default=300, description='Temperature of the simulation.'):
            NUMBER,
        Required('pressure_min', default=0.001, description='Lower pressure to sample (bar).'):
            NUMBER,
        Required('pressure_max', default=30, description='Upper pressure to sample (bar).'):
            NUMBER,
        Required('pressure_maxstep', default=5.0, description='(float) Max distance between pressure points (bar).'):
            NUMBER,
        Required('pressure_precision',
                 default=0.1,
                 description='Precision in the sampling of the isotherm: 0.1 ok, 0.05 for high resolution.'):
            NUMBER,
    })
    parameters_info = parameters_schema.schema  # shorthand for printing

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(IsothermWorkChain, exclude=['molecule', 'parameters'])
        spec.input('parameters',
                   valid_type=Dict,
                   default=lambda: Dict(dict=cls.parameters_schema({})),
                   validator=functools.partial(validate_dict, schema=cls.parameters_schema),
                   help='Parameters for Isotherm work chain')

        spec.input('pe_parameters',
                   valid_type=Dict,
                   default=lambda: Dict(dict=PE_PARAMETERS_DEFAULT),
                   help='Parameters for PE process modelling')

        spec.outline(
            cls.run_isotherms,
            cls.run_calcpe,
        )

        spec.expose_outputs(IsothermWorkChain, namespace='co2')
        spec.expose_outputs(IsothermWorkChain, namespace='n2')
        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='Output parmaters of a calc_PE calculations')

    def run_isotherms(self):
        """Run Isotherm work chain for CO2 and N2."""

        inputs = self.exposed_inputs(IsothermWorkChain)
        self.report('Run Isotherm work chain for CO2 and N2, in CifData<{}> (label: {} )'.format(
            self.inputs.structure.pk, self.inputs.structure.label))

        for mol in ['co2', 'n2']:

            dict_merge(
                inputs, {
                    'metadata': {
                        'call_link_label': 'run_isotherm_for_{}'.format(mol),
                    },
                    'molecule': Str(mol),
                    'parameters': self.inputs.parameters
                })

            running = self.submit(IsothermWorkChain, **inputs)
            self.to_context(**{'isotherm_{}'.format(mol): running})

    def run_calcpe(self):
        """Expose isotherm outputs, prepare calc_pe, run it and return the output."""

        self.out_many(self.exposed_outputs(self.ctx.isotherm_co2, IsothermWorkChain, namespace='co2'))
        self.out_many(self.exposed_outputs(self.ctx.isotherm_n2, IsothermWorkChain, namespace='n2'))

        calcpe = calc_co2_parasitic_energy(isot_co2=self.ctx.isotherm_co2.outputs['output_parameters'],
                                           isot_n2=self.ctx.isotherm_n2.outputs['output_parameters'],
                                           pe_parameters=self.inputs['pe_parameters'])

        self.out('output_parameters', calcpe)
        self.report('calc_pe Dict<{}>'.format(calcpe.pk))
