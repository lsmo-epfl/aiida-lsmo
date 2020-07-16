"""RosenbluthWeight work chain."""

import os

from aiida.plugins import WorkflowFactory, CalculationFactory
from aiida.orm import Dict, Str, Float
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext
from aiida_lsmo.utils import aiida_dict_merge

# import sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# import calculations
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name

# Deafault parameters
PARAMETERS_DEFAULT = {  #TODO: create IsothermParameters instead of Dict # pylint: disable=fixme
    "temperature": 300,  # (float) Temperature of the simulation.
    "raspa_cycles": int(1e6),  # (int) Number of Widom cycles.
    "raspa_verbosity": 10,  # (int) Print stats every: number of cycles / raspa_verbosity.
}


# calcfunctions (in order of appearence)
@calcfunction
def get_molecule_dict(molecule_name):
    """Get a Dict from the isotherm_molecules.yaml"""
    import ruamel.yaml as yaml
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, "isotherm_data", "isotherm_molecules.yaml")
    with open(yamlfile, 'r') as stream:
        yaml_dict = yaml.safe_load(stream)
    molecule_dict = yaml_dict[molecule_name.value]
    return Dict(dict=molecule_dict)


@calcfunction
def get_ff_parameters(molecule_dict):
    """Get the parameters for ff_builder."""
    ff_params = {}
    ff_params['ff_framework'] = 'NONE'
    ff_params['ff_molecules'] = {molecule_dict['name']: molecule_dict['forcefield']}
    ff_params['shifted'] = False
    ff_params['tail_corrections'] = False
    ff_params['mixing_rule'] = 'Lorentz-Berthelot'
    ff_params['separate_interactions'] = False
    return Dict(dict=ff_params)


@calcfunction
def get_rosenbluth_weight(out_dict):
    """Get the parameters for ff_builder."""
    comp_dict = list(out_dict['big_box']['components'].values())[0]
    rbw = Float(comp_dict['widom_rosenbluth_factor_average'])
    rbw.label = 'Rosenbluth Weight'
    rbw.description = 'Relative error: {:.3f}%'.format(100 * comp_dict['widom_rosenbluth_factor_dev'] / rbw.value)
    return rbw


class RosenbluthWeightWorkChain(WorkChain):
    """Workchain that computes the Rosenbluth Weight at a given temperature, for a flexible molecule in vacuum."""

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        spec.input("molecule",
                   valid_type=(Str, Dict),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')

        spec.input("parameters",
                   valid_type=Dict,
                   help='Parameters for the Isotherm workchain: will be merged with IsothermParameters_defaults.')

        spec.outline(
            cls.setup,
            cls.run_raspa_widom,
            cls.return_rosenbluth_weight,
        )

        spec.output('rosenbluth_weight', valid_type=Float, required=True, help='Rosenbluth weight value.')

        spec.output('rosenbluth_weight_rel_err',
                    valid_type=Float,
                    required=True,
                    help='Rosenbluth weight relative error.')

    def setup(self):
        """Initialize the parameters"""

        # Get the molecule Dict from the yaml or directly as an input
        if isinstance(self.inputs.molecule, Str):
            self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
        elif isinstance(self.inputs.molecule, Dict):
            self.ctx.molecule = self.inputs.molecule

        # Get the parameters Dict, merging defaults with user settings
        self.ctx.parameters = aiida_dict_merge(Dict(dict=PARAMETERS_DEFAULT), self.inputs.parameters)

    def _get_widom_param(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        param = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfInitializationCycles": 0,
                "NumberOfCycles": self.ctx.parameters['raspa_cycles'],
                "PrintPropertiesEvery": self.ctx.parameters['raspa_cycles'] / self.ctx.parameters['raspa_verbosity'],
                "PrintEvery": int(1e10),
                "Forcefield": "Local",
                "CutOff": 20.00,
            },
            "System": {
                "big_box": {
                    "type": "Box",
                    "BoxLengths": "100 100 100",
                    "ExternalTemperature": self.ctx.parameters['temperature'],
                },
            },
            "Component": {
                self.ctx.molecule['name']: {
                    "MoleculeDefinition": "Local",
                    "WidomProbability": 1.0,
                },
            },
        }

        if self.ctx.molecule['charged']:  # NOTE: `Chargemethod Ewald` is the default in Raspa!
            param["GeneralSettings"].update({"ChargeMethod": "CoulombTruncated", "CutOffChargeCharge": 20})
        else:
            param["GeneralSettings"].update({"ChargeMethod": "None"})

        return param

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa, with the molecule in vacuum."""

        # Initialize the input for raspa_base
        self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.inp['metadata']['label'] = "RaspaWidom"
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_widom"
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self._get_widom_param())

        # Generate the force field with the ff_builder
        ff_params = get_ff_parameters(self.ctx.molecule)
        files_dict = FFBuilder(ff_params)
        self.ctx.inp['raspa']['file'] = files_dict

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)

        return ToContext(raspa_widom=running)

    def return_rosenbluth_weight(self):
        """Extract the Rosenbluth Weight from the Raspa calculation and return it as output."""

        self.out("rosenbluth_weight", get_rosenbluth_weight(self.ctx.raspa_widom.outputs.output_parameters))

        self.report("Rosenbluth weight: {} Float<{}>".format(self.outputs['rosenbluth_weight'].value,
                                                             self.outputs['rosenbluth_weight'].pk))

        self.report(self.outputs['rosenbluth_weight'].description)
