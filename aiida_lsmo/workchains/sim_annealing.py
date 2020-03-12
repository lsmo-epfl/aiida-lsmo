"""Isotherm workchain"""

import os
import ruamel.yaml as yaml
import ase
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, append_, while_
from aiida_lsmo.utils import check_resize_unit_cell, aiida_dict_merge, aiida_cif_merge

# import sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# import calculations
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  # pylint: disable=invalid-name

# Deafault parameters
PARAMETERS_DEFAULT = {
    "ff_framework": "UFF",  # (str) Forcefield of the structure.
    "ff_separate_interactions": False,  # (bool) Use "separate_interactions" in the FF builder.
    "ff_mixing_rule": "Lorentz-Berthelot",  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
    "ff_tail_corrections": True,  # (bool) Apply tail corrections.
    "ff_shifted": False,  # (bool) Shift or truncate the potential at cutoff.
    "ff_cutoff": 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
    "temperature_list": [300, 250, 200, 250, 100, 50],  # (list) List of decreasing temperatures for the annealing.
    "mc_steps": int(1e3),  # (int) Number of MC cycles.
    "number_of_molecules": 1  # (int) Number of molecules loaded in the framework.
}


# calcfunctions (in order of appearence)
@calcfunction
def get_molecule_dict(molecule_name):
    """Get a Dict from the isotherm_molecules.yaml"""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, "isotherm_data", "isotherm_molecules.yaml")
    with open(yamlfile, 'r') as stream:
        yaml_dict = yaml.safe_load(stream)
    molecule_dict = yaml_dict[molecule_name.value]
    return Dict(dict=molecule_dict)


def get_ff_parameters(molecule_dict, isotparam):
    """Get the parameters for ff_builder."""
    ff_params = {}
    ff_params['ff_framework'] = isotparam['ff_framework']
    ff_params['ff_molecules'] = {molecule_dict['name']: molecule_dict['forcefield']}
    ff_params['shifted'] = isotparam['ff_shifted']
    ff_params['tail_corrections'] = isotparam['ff_tail_corrections']
    ff_params['mixing_rule'] = isotparam['ff_mixing_rule']
    ff_params['separate_interactions'] = isotparam['ff_separate_interactions']
    return Dict(dict=ff_params)


def load_yaml():
    """ Load the ff_data.yaml as a dict."""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfullpath = os.path.join(thisdir, '..', 'calcfunctions', 'ff_data.yaml')
    with open(yamlfullpath, 'r') as stream:
        ff_data = yaml.safe_load(stream)
    return ff_data


@calcfunction
def get_molecule_from_restart_file(structure_cif, molecule_folderdata, input_dict, molecule_dict):
    """Get a CifData file having the cell of the initial (unexpanded) structure and the geometry of the loaded molecule.
    TODO: this is source of error if there are more than one molecule, and the cell has been expandes,
    as you can not wrap them in the small cell.
    """

    # Get number of guest molecules
    if "number_of_molecules" in input_dict.get_dict():
        number_of_molecules = input_dict["number_of_molecules"]
    else:
        number_of_molecules = 1

    # Get the non-M (non-dummy) atom types of the molecule (ASE accepts only atomic elements as "symbols")
    ff_data = load_yaml()
    ff_data_molecule = ff_data[molecule_dict['name']][molecule_dict['forcefield']]
    symbols = [x[0].split("_")[0] for x in ff_data_molecule['atomic_positions']]
    symbols = [x for x in symbols if x != 'M']
    symbols *= number_of_molecules

    # Get the coordinates of the molecule in the extended uni cell
    restart_fname = molecule_folderdata._repository.list_object_names(os.path.join('Restart', 'System_0'))[0]  # pylint: disable=protected-access
    restart_abs_path = os.path.join(
        molecule_folderdata._repository._get_base_folder().abspath,  # pylint: disable=protected-access
        'Restart',
        'System_0',
        restart_fname)

    positions = []
    with open(restart_abs_path, "r") as fobj:
        for line in fobj:
            if 'Adsorbate-atom-position:' in line:
                positions.append(line.split()[3:6])

    # Get the cell, combine the ASE object and return the CifData
    cell = structure_cif.get_ase().get_cell()
    molecule_ase = ase.Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
    molecule_ase.wrap()  # wrap molecule in the unit cell
    return CifData(ase=molecule_ase)


@calcfunction
def get_output_parameters(input_dict, min_out_dict, **nvt_out_dict):
    """Merge energy info from the calculations."""

    if "number_of_molecules" in input_dict.get_dict():
        number_of_molecules = input_dict["number_of_molecules"]
    else:
        number_of_molecules = 1

    if "temperature_list" in input_dict.get_dict():
        temperature_list = input_dict["temperature_list"]
    else:
        temperature_list = [300, 250, 200, 250, 100, 50]

    out_dict = {'number_of_molecules': number_of_molecules, 'description': [], 'energy_unit': 'kJ/mol'}

    key_list = [
        'energy_host/ads_tot_final',
        'energy_host/ads_vdw_final',
        'energy_host/ads_coulomb_final',
        'energy_ads/ads_tot_final',
        'energy_ads/ads_vdw_final',
        'energy_ads/ads_coulomb_final',
    ]

    for key in key_list:
        out_dict[key] = []

    for i, temp in enumerate(temperature_list):
        out_dict['description'].append('NVT simulation at {} K'.format(temp))
        nvt_out_dict_i = nvt_out_dict['RaspaNVT_{}'.format(i + 1)]
        for key in key_list:
            out_dict[key].append(nvt_out_dict_i["framework_1"]['general'][key])

    out_dict['description'].append('Final energy minimization')
    for key in key_list:
        out_dict[key].append(min_out_dict["framework_1"]['general'][key])

    return Dict(dict=out_dict)


class SimAnnealingWorkChain(WorkChain):
    """A work chain to compute the minimum energy geometry of a molecule inside a framework, using simulated annealing,
    i.e., decreasing the temperature of a Monte Carlo simulation and finally running and energy minimization step.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        spec.input('structure', valid_type=CifData, help='Adsorbent framework CIF.')

        spec.input("molecule",
                   valid_type=(Str, Dict),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')

        spec.input("parameters",
                   valid_type=Dict,
                   help='Parameters for the SimAnnealing workchain: will be merged with default ones.')

        spec.outline(cls.setup, while_(cls.should_run_nvt)(cls.run_raspa_nvt,), cls.run_raspa_min, cls.return_results)

        spec.output('loaded_molecule', valid_type=CifData, help='CIF containing the final postition of the molecule.')

        spec.output('loaded_structure', valid_type=CifData, help='CIF containing the loaded structure.')

        spec.output(
            'output_parameters',
            valid_type=Dict,
            required=False,  #DO later: this needs to implement the parsing of final energies in aiida-raspa
            help='Information about the final configuration.')

    def _get_raspa_nvt_param(self):
        """Write Raspa input parameters from scratch, for an MC NVT calculation"""

        param = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfCycles": self.ctx.parameters['mc_steps'],
                "PrintPropertiesEvery": int(1e10),
                "PrintEvery": self.ctx.parameters['mc_steps'] / 100,
                "RemoveAtomNumberCodeFromLabel":
                    True,  # BE CAREFULL: needed in AiiDA-1.0.0 because of github.com/aiidateam/aiida-core/issues/3304
                "Forcefield": "Local",
                "UseChargesFromCIFFile": "yes",
                "CutOff": self.ctx.parameters['ff_cutoff'],
            },
            "System": {
                "framework_1": {
                    "type": "Framework",
                }
            },
            "Component": {
                self.ctx.molecule['name']: {
                    "MoleculeDefinition": "Local",
                    "TranslationProbability": 1.0,
                    "ReinsertionProbability": 1.0,
                    "CreateNumberOfMolecules": self.ctx.parameters['number_of_molecules'],
                },
            },
        }

        # Check particular conditions and settings
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
        param["System"]["framework_1"]["UnitCells"] = "{} {} {}".format(mult[0], mult[1], mult[2])

        if self.ctx.inp["raspa"]["block_pocket"]:  #TODO: check if working properly! # pylint: disable=fixme
            param["Component"][self.ctx.molecule['name']]["BlockPocketsFileName"] = "block_file"

        if not self.ctx.molecule['singlebead']:
            param["Component"][self.ctx.molecule['name']].update({"RotationProbability": 1.0})

        if self.ctx.molecule['charged']:  # NOTE: `Chargemethod Ewald` is the default in Raspa!
            param["GeneralSettings"].update({"ChargeMethod": "Ewald", "EwaldPrecision": 1e-6})
        else:
            param["GeneralSettings"].update({"ChargeMethod": "None"})
        return param

    def setup(self):
        """Initialize the parameters"""

        # Get the molecule Dict from the yaml or directly as an input
        if isinstance(self.inputs.molecule, Str):
            self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
        elif isinstance(self.inputs.molecule, Dict):
            self.ctx.molecule = self.inputs.molecule

        # Get the parameters Dict, merging defaults with user settings
        self.ctx.parameters = aiida_dict_merge(Dict(dict=PARAMETERS_DEFAULT), self.inputs.parameters)

        # Initialize the input for raspa_base, which later will need only minor updates
        self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.inp['raspa']['framework'] = {"framework_1": self.inputs.structure}
        self.ctx.raspa_param = self._get_raspa_nvt_param()

        # Generate the force field with the ff_builder
        ff_params = get_ff_parameters(self.ctx.molecule, self.ctx.parameters)
        files_dict = FFBuilder(ff_params)
        self.ctx.inp['raspa']['file'] = files_dict

        self.ctx.count = 0

    def should_run_nvt(self):
        """Update temperature untill the last of the list."""
        if self.ctx.count == 1:  #Placed here to work even with just one temperature in temperature_list
            self.ctx.raspa_param["Component"][self.ctx.molecule['name']]["CreateNumberOfMolecules"] = 0
        return self.ctx.count < len(self.ctx.parameters['temperature_list'])

    def run_raspa_nvt(self):
        """Run a NVT calculation in Raspa."""

        temperature = self.ctx.parameters['temperature_list'][self.ctx.count]
        self.ctx.raspa_param["System"]["framework_1"]["ExternalTemperature"] = temperature
        if self.ctx.count > 0:
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_nvt[self.ctx.count - 1].outputs.retrieved

        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
        self.ctx.inp['metadata']['label'] = "RaspaNVT_{}".format(self.ctx.count + 1)
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_nvt_{}".format(self.ctx.count + 1)

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa NVT ({} of {})".format(self.ctx.count + 1,
                                                          len(self.ctx.parameters['temperature_list'])))
        self.ctx.count += 1

        return ToContext(raspa_nvt=append_(running))

    def run_raspa_min(self):
        """Run a Energy Minimization in Raspa."""

        # Update parameters Dict and labels
        self.ctx.raspa_param['GeneralSettings'].update({
            'SimulationType': 'Minimization',
            'NumberOfInitializationCycles': 0,  # no NVT steps
            'NumberOfCycles': 1,  # do only one Minimization
            'PrintEvery': 1,  # not readen: it prints anyway all min iterations
            'MaximumNumberOfMinimizationSteps': 10000,  # it will stop before, if converged
            'RMSGradientTolerance': 1e-6,
            'MaxGradientTolerance': 1e-6,
        })
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
        self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_nvt[self.ctx.count - 1].outputs.retrieved
        self.ctx.inp['metadata']['label'] = "RaspaMin"
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_min"

        # Create the calculation process, launch it and update pressure index
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa final minimization")
        return ToContext(raspa_min=running)

    def return_results(self):
        """Return molecule position and energy info."""

        self.out(
            "loaded_molecule",
            get_molecule_from_restart_file(self.inputs.structure, self.ctx.raspa_min.outputs.retrieved,
                                           self.inputs.parameters, self.ctx.molecule))
        self.out("loaded_structure", aiida_cif_merge(self.inputs.structure, self.outputs['loaded_molecule']))

        nvt_out_dict = {}
        for calc in self.ctx.raspa_nvt:
            nvt_out_dict[calc.label] = calc.outputs.output_parameters

        self.out(
            "output_parameters",
            get_output_parameters(input_dict=self.inputs.parameters,
                                  min_out_dict=self.ctx.raspa_min.outputs.output_parameters,
                                  **nvt_out_dict))

        self.report(
            "Work chain competed! Molecule CifData<{}>, loaded structure CifData<{}>, output parameters Dict<{}>".
            format(self.outputs['loaded_molecule'].pk, self.outputs['loaded_structure'].pk,
                   self.outputs['output_parameters'].pk))
