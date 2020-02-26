"""IsothermCalcPE work chain."""

from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida.engine import WorkChain
from aiida_lsmo.utils import dict_merge
from aiida_lsmo.calcfunctions import PE_PARAMETERS_DEFAULT, calc_co2_parasitic_energy

# import sub-workchains
IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  #pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  #pylint: disable=invalid-name

ISOTHERM_PARAMETERS_DEFAULT = {  # Parameters used in 10.1021/acscentsci.9b00619
    "ff_framework": "UFF",  # valid_type=Str, help='Forcefield of the structure.'
    "ff_shifted": True,  # shift or truncate at cutoff
    "ff_tail_corrections": False,  # apply tail corrections
    "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
    "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
    "ff_cutoff": 12.0,  # valid_type=Float, help='CutOff truncation for the VdW interactions (Angstrom)'
    "temperature": 300,  # valid_type=Float, help='Temperature of the simulation'
    "zeopp_volpo_samples": 1e5,  # valid_type=Int,help='Number of samples for VOLPO calculation (per UC volume)'
    "zeopp_block_samples": 100,  # valid_type=Int, help='Number of samples for BLOCK calculation (per A^3)'
    "raspa_minKh": 1e-10,
    "raspa_verbosity": 10,  # valid_type=Int,help='Print stats every: number of cycles / raspa_verbosity'
    "raspa_widom_cycles": 1e5,  # valid_type=Int, help='Number of widom cycles'
    "raspa_gcmc_init_cycles": 1e3,  # valid_type=Int, help='Number of GCMC initialization cycles'
    "raspa_gcmc_prod_cycles": 1e4,  # valid_type=Int, help='Number of GCMC production cycles'
    "pressure_precision": 0.1,
    "pressure_maxstep": 5,  # valid_type=Float, help='Max distance between pressure points (bar)'
    "pressure_min": 0.001,  # valid_type=Float, help='Lower pressure to sample (bar)'
    "pressure_max": 30
}


class IsothermCalcPEWorkChain(WorkChain):
    """ Compute CO2 parassitic energy (PE) after running IsothermWorkChain for CO2 and N2 at 300K."""

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(IsothermWorkChain, exclude=['molecule', 'parameters'])
        spec.input('parameters',
                   valid_type=Dict,
                   default=lambda: Dict(dict=ISOTHERM_PARAMETERS_DEFAULT),
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
        self.report("Run Isotherm work chain for CO2 and N2, in CifData<{}> (label: {} )".format(
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

        self.out("output_parameters", calcpe)
        self.report('calc_pe Dict<{}>'.format(calcpe.pk))
