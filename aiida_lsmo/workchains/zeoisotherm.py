# -*- coding: utf-8 -*-
"""ZeoIsotherm work chain"""

import os

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Float, Str, List, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, append_, while_, if_
from aiida_lsmo.utils import aiida_dict_merge
from aiida_lsmo.utils import dict_merge

# import sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# import calculations
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  # pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name

# Deafault parameters
ISOTHERMPARAMETERS_DEFAULT = {  #TODO: create IsothermParameters instead of Dict # pylint: disable=fixme
    'ff_framework': 'Zeo-test',  # (str) Forcefield of the structure.
    'ff_separate_interactions': False,  # (bool) Use "separate_interactions" in the FF builder.
    'ff_mixing_rule': 'Lorentz-Berthelot',  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
    'ff_tail_corrections': True,  # (bool) Apply tail corrections.
    'ff_shifted': False,  # (bool) Shift or truncate the potential at cutoff.
    'ff_cutoff': 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
    'cations_equilibration_temp_list': [
        1000, 800, 600, 500, 400
    ],  # (list) List of temperatures for the equilibrations of cations: last will be "temperature".
    'temperature': 300,  # (float) Temperature of the simulation.
    'temperature_list': None,  # (list) To be used by IsothermMultiTempWorkChain.
    'zeopp_volpo_samples': int(1e5),  # (int) Number of samples for VOLPO calculation (per UC volume).
    'zeopp_block_samples': int(100),  # (int) Number of samples for BLOCK calculation (per A^3).
    'raspa_minKh': 1e-10,  # (float) If Henry coefficient < raspa_minKh do not run the isotherm (mol/kg/Pa).
    'raspa_verbosity': 10,  # (int) Print stats every: number of cycles / raspa_verbosity.
    'raspa_widom_cycles': int(1e5),  # (int) Number of Widom cycles.
    'raspa_nvt_cations_cycles': int(1e4),  # (int) Number of GCMC production cycles.
    'raspa_gcmc_init_cycles': int(1e3),  # (int) Number of GCMC initialization cycles.
    'raspa_gcmc_prod_cycles': int(1e4),  # (int) Number of GCMC production cycles.
    'pressure_list': None,  # (list) Pressure list for the isotherm (bar): if given it will skip to guess it.
    'pressure_precision': 0.1,  # (float) Precision in the sampling of the isotherm: 0.1 ok, 0.05 for high resolution.
    'pressure_maxstep': 5,  # (float) Max distance between pressure points (bar).
    'pressure_min': 0.001,  # (float) Lower pressure to sample (bar).
    'pressure_max': 10  # (float) Upper pressure to sample (bar).
}


# calcfunctions (in order of appearence)
@calcfunction
def get_molecule_dict(molecule_name):
    """Get a Dict from the isotherm_molecules.yaml"""
    import ruamel.yaml as yaml
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, 'isotherm_data', 'isotherm_molecules.yaml')
    with open(yamlfile, 'r') as stream:
        yaml_dict = yaml.safe_load(stream)
    molecule_dict = yaml_dict[molecule_name.value]
    return Dict(dict=molecule_dict)


@calcfunction
def get_zeopp_parameters(molecule_dict, isotparam):
    """Get the ZeoppParameters from the inputs of the workchain"""
    probe_rad = molecule_dict['proberad']
    param_dict = {
        'ha': 'DEF',
        'volpo': [probe_rad, probe_rad, isotparam['zeopp_volpo_samples']],
        'block': [probe_rad, isotparam['zeopp_block_samples']],
    }
    return ZeoppParameters(dict=param_dict)


@calcfunction
def get_ff_parameters(molecule_dict, cation_name, isotparam):
    """Get the parameters for ff_builder.
    NOTE: for the cation consider a ff with the same name as the ff_framework.
    """
    ff_params = {}
    ff_params['ff_framework'] = isotparam['ff_framework']
    ff_params['ff_molecules'] = {
        molecule_dict['name']: molecule_dict['forcefield'],
        cation_name.value: isotparam['ff_framework']
    }
    ff_params['shifted'] = isotparam['ff_shifted']
    ff_params['tail_corrections'] = isotparam['ff_tail_corrections']
    ff_params['mixing_rule'] = isotparam['ff_mixing_rule']
    ff_params['separate_interactions'] = isotparam['ff_separate_interactions']
    return Dict(dict=ff_params)


@calcfunction
def choose_pressure_points(inp_param, geom, raspa_widom_out):
    """If 'presure_list' is not provide, model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm, in a List.
    """
    if inp_param['pressure_list']:
        pressure_points = inp_param['pressure_list']
    else:
        khenry = list(raspa_widom_out['framework_1']['components'].values())[1]['henry_coefficient_average']  #mol/kg/Pa
        b_value = khenry / geom['Estimated_saturation_loading'] * 1e5  #(1/bar)
        pressure_points = [inp_param['pressure_min']]
        while True:
            pold = pressure_points[-1]
            delta_p = min(inp_param['pressure_maxstep'],
                          inp_param['pressure_precision'] * (b_value * pold**2 + 2 * pold + 1 / b_value))
            pnew = pold + delta_p
            if pnew <= inp_param['pressure_max']:
                pressure_points.append(pnew)
            else:
                pressure_points.append(inp_param['pressure_max'])
                break
    return List(list=pressure_points)


@calcfunction
def get_geometric_dict(zeopp_out, molecule):
    """Return the geometric Dict from Zeopp results, including Qsat and is_porous"""
    geometric_dict = zeopp_out.get_dict()
    geometric_dict.update({
        'Estimated_saturation_loading': zeopp_out['POAV_cm^3/g'] * molecule['molsatdens'],
        'Estimated_saturation_loading_unit': 'mol/kg',
        'is_porous': geometric_dict['POAV_A^3'] > 0.000
    })
    return Dict(dict=geometric_dict)


@calcfunction
def get_output_parameters(geom_out, inp_params, widom_out=None, pressures=None, **gcmc_out_dict):
    """Merge results from all the steps of the work chain."""

    out_dict = geom_out.get_dict()

    if out_dict['is_porous'] and widom_out:
        widom_out_mol = list(widom_out['framework_1']['components'].values())[1]

        out_dict.update({
            'temperature': inp_params['temperature'],
            'temperature_unit': 'K',
            'is_kh_enough': widom_out_mol['henry_coefficient_average'] > inp_params['raspa_minKh']
        })

        widom_labels = [
            'henry_coefficient_average',
            'henry_coefficient_dev',
            'henry_coefficient_unit',
            'adsorption_energy_widom_average',
            'adsorption_energy_widom_dev',
            'adsorption_energy_widom_unit',
        ]

        for label in widom_labels:
            out_dict.update({label: widom_out_mol[label]})

        if out_dict['is_kh_enough']:

            isotherm = {
                'pressure': pressures,
                'pressure_unit': 'bar',
                'loading_absolute_average': [],
                'loading_absolute_dev': [],
                'loading_absolute_unit': 'mol/kg',
                'enthalpy_of_adsorption_average': [],
                'enthalpy_of_adsorption_dev': [],
                'enthalpy_of_adsorption_unit': 'kJ/mol'
            }

            conv_ener = 1.0 / 120.273  # K to kJ/mol
            for i in range(len(pressures)):
                gcmc_out = gcmc_out_dict['RaspaGCMC_{}'.format(i + 1)]['framework_1']
                gcmc_out_mol = list(gcmc_out['components'].values())[1]
                conv_load = gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg']

                for label in ['loading_absolute_average', 'loading_absolute_dev']:
                    isotherm[label].append(conv_load * gcmc_out_mol[label])

                for label in ['enthalpy_of_adsorption_average', 'enthalpy_of_adsorption_dev']:
                    if gcmc_out['general'][label]:
                        isotherm[label].append(conv_ener * gcmc_out['general'][label])
                    else:  # when there are no particles and Raspa return Null enthalpy
                        isotherm[label].append(None)

            out_dict.update({
                'isotherm': isotherm,
                'conversion_factor_molec_uc_to_cm3stp_cm3': gcmc_out_mol['conversion_factor_molec_uc_to_cm3stp_cm3'],
                'conversion_factor_molec_uc_to_mg_g': gcmc_out_mol['conversion_factor_molec_uc_to_mg_g'],
                'conversion_factor_molec_uc_to_mol_kg': gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg'],
            })

    return Dict(dict=out_dict)


class ZeoIsothermWorkChain(WorkChain):
    """Isotherm work chain designed specifically for zeolites: both all-silica and Al-substituted.
    """

    @classmethod
    def define(cls, spec):
        super(ZeoIsothermWorkChain, cls).define(spec)

        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])

        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        spec.input('structure', valid_type=CifData, help='Adsorbent framework CIF.')

        spec.input('sial_ratio', valid_type=Float, default=Float(1), help='Get the Si/Al ratio in the zeolite.')

        spec.input('molecule',
                   valid_type=(Str, Dict),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')

        spec.input('cation', valid_type=Str, default=Str('Na'), help='Cation to be used for charged zeolites (Na/Ca)')
        #TODO: need to implement Ca+2 # pylint: disable=fixme

        spec.input('parameters',
                   valid_type=Dict,
                   help='Parameters for the Isotherm workchain: will be merged with IsothermParameters_defaults.')

        spec.input('geometric',
                   valid_type=Dict,
                   required=False,
                   help='[Only used by IsothermMultiTempWorkChain] Already computed geometric properties')

        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            while_(cls.should_equilibrate_cations)(cls.run_equilibrate_cations,),
            if_(cls.should_run_widom)(  # run Widom only if porous
                cls.run_raspa_widom,  # run raspa widom calculation
                if_(cls.should_run_gcmc)(  # Kh is high enough
                    cls.init_raspa_gcmc,  # initializate setting for GCMC
                    while_(cls.should_run_another_gcmc)(  # new pressure
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation
                    ),
                ),
            ),
            cls.return_output_parameters,
        )

        spec.expose_outputs(ZeoppCalculation, include=['block'])  #only if porous

        spec.output(
            'output_parameters',
            valid_type=Dict,
            required=True,
            help='Results of the single temperature wc: keys can vay depending on is_porous and is_kh_enough booleans.')

    def setup(self):
        """Initialize the parameters"""

        # Get the molecule Dict from the yaml or directly as an input
        if isinstance(self.inputs.molecule, Str):
            self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
        elif isinstance(self.inputs.molecule, Dict):
            self.ctx.molecule = self.inputs.molecule

        # Get the parameters Dict, merging defaults with user settings
        self.ctx.parameters = aiida_dict_merge(Dict(dict=ISOTHERMPARAMETERS_DEFAULT), self.inputs.parameters)

        # Get integer temperature in context for easy reports
        self.ctx.temperature = int(round(self.ctx.parameters['temperature']))

        # Get the number of cations needed and other data for the cation equilibration
        natoms = self.inputs.structure.get_ase().get_global_number_of_atoms()
        ntatoms = natoms / 3
        self.ctx.ncations = int(ntatoms / (self.inputs.sial_ratio + 1))
        self.ctx.current_equilibration_id = 0
        self.ctx.cations_equilibration_temp_list = self.ctx.parameters['cations_equilibration_temp_list'] + [
            self.ctx.parameters['temperature']
        ]

        # Understand if IsothermMultiTempWorkChain is calling this work chain
        if 'geometric' in self.inputs:
            self.ctx.multitemp_mode = 'run_single_temp'
        elif self.ctx.parameters['temperature_list']:
            self.ctx.multitemp_mode = 'run_geom_only'
        else:
            self.ctx.multitemp_mode = None

    def run_zeopp(self):
        """Perform Zeo++ block and VOLPO calculations."""

        # Skip zeopp calculation if the geometric properties are already provided by IsothermMultiTemp
        if self.ctx.multitemp_mode == 'run_single_temp':
            return None

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')

        # Set inputs for zeopp
        dict_merge(
            inputs, {
                'metadata': {
                    'label': 'ZeoppVolpoBlock',
                    'call_link_label': 'run_zeopp_block_and_volpo',
                },
                'structure':
                    self.inputs.structure,
                'atomic_radii':
                    SinglefileData(
                        file=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'isotherm_data', 'DEFAULT.rad')),
                'parameters':
                    get_zeopp_parameters(self.ctx.molecule, self.ctx.parameters)
            })

        running = self.submit(ZeoppCalculation, **inputs)
        self.report('Running zeo++ block and volpo Calculation<{}>'.format(running.id))
        return ToContext(zeopp=running)

    def _get_nvt_cations(self):
        """Write Raspa input parameters from scratch, for a NVT calculation"""

        param = {
            'GeneralSettings': {
                'SimulationType': 'MonteCarlo',
                'NumberOfInitializationCycles': self.ctx.parameters['raspa_nvt_cations_cycles'] / 2,
                'NumberOfCycles': self.ctx.parameters['raspa_nvt_cations_cycles'] / 2,
                'PrintPropertiesEvery': int(1e10),
                'PrintEvery': self.ctx.parameters['raspa_nvt_cations_cycles'] / self.ctx.parameters['raspa_verbosity'],
                'RemoveAtomNumberCodeFromLabel': True,
                'ModifyOxgensConnectedToAluminium': True,
                'Forcefield': 'Local',
                'CutOff': self.ctx.parameters['ff_cutoff'],
                'ChargeMethod': 'Ewald',
                'EwaldPrecision': 1e-6,
            },
            'System': {
                'framework_1': {
                    'type': 'Framework',
                    'UnitCells': '1 1 1',  #NOTE: it is expanded before!
                    'ExternalTemperature': self.ctx.cations_equilibration_temp_list[0],
                }
            },
            'Component': {
                self.inputs.cation.value: {
                    'MoleculeDefinition': 'Local',
                    'ExtraFrameworkMolecule': True,
                    'TranslationProbability': 1.0,
                    'RandomTranslationProbability': 1.0,  #Raspa warns to use this instead of ReinsertionProbability
                    'CreateNumberOfMolecules': self.ctx.ncations,
                },
            },
        }
        return param

    def should_equilibrate_cations(self):
        """Decide whether to run the cations equilibration."""
        if self.ctx.ncations == 0:
            return False
        return self.ctx.current_equilibration_id < len(self.ctx.cations_equilibration_temp_list)

    def run_equilibrate_cations(self):
        """Run NVT for the equilibration of the cations as a list of calculations."""

        if self.ctx.current_equilibration_id == 0:

            # Initialize the input for raspa_base, which later will need only minor updates for GCMC
            self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')

            # Create input parameters and pass the framework
            self.ctx.raspa_param = self._get_nvt_cations()
            self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
            self.ctx.inp['raspa']['framework'] = {'framework_1': self.inputs.structure}

            # Generate the force field with the ff_builder
            ff_params = get_ff_parameters(self.ctx.molecule, self.inputs.cation, self.ctx.parameters)
            files_dict = FFBuilder(ff_params)
            self.ctx.inp['raspa']['file'] = files_dict

        elif self.ctx.current_equilibration_id > 0:

            # Update Parameters
            self.ctx.raspa_param['System']['framework_1'][
                'ExternalTemperature'] = self.ctx.cations_equilibration_temp_list[self.ctx.current_equilibration_id]
            self.ctx.raspa_param['Component'][self.inputs.cation.value]['CreateNumberOfMolecules'] = 0
            self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

            # Update restart (if present, i.e., if current_p_index>0)
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_nvt_cations[
                self.ctx.current_equilibration_id - 1].outputs.retrieved

        # Update metadata and submit
        self.ctx.inp['metadata']['label'] = 'RaspaNVTCations{}'.format(self.ctx.current_equilibration_id)
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_cations_{}'.format(self.ctx.current_equilibration_id)
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report('Running Raspa NVT @ {}K for equilibrating cations ({}/{})'.format(
            int(self.ctx.cations_equilibration_temp_list[self.ctx.current_equilibration_id]),
            self.ctx.current_equilibration_id + 1, len(self.ctx.cations_equilibration_temp_list)))
        self.ctx.current_equilibration_id += 1
        return ToContext(raspa_nvt_cations=append_(running))

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume,
        also check the number of blocking spheres and estimate the saturation loading.
        Also, stop if called by IsothermMultiTemp for geometric results only."""

        # Get geometric properties and consider if IsothermMultiTempWorkChain is calling this workchain
        if self.ctx.multitemp_mode == 'run_single_temp':
            self.ctx.geom = self.inputs.geometric
            return True
        self.ctx.geom = get_geometric_dict(self.ctx.zeopp.outputs.output_parameters, self.ctx.molecule)

        if self.ctx.geom['is_porous']:
            self.report('Found accessible pore volume: continue')
            self.report('Found {} blocking spheres'.format(self.ctx.geom['Number_of_blocking_spheres']))
            # Return block file only if blocking spheres are present
            if self.ctx.geom['Number_of_blocking_spheres'] > 0:
                self.out_many(self.exposed_outputs(self.ctx.zeopp, ZeoppCalculation))
        else:
            self.report('No accessible pore volume: stop')

        return self.ctx.geom['is_porous'] and not self.ctx.multitemp_mode == 'run_geom_only'

    def _update_param_for_widom(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        param = self.ctx.raspa_param
        param['GeneralSettings'].update({
            'NumberOfInitializationCycles': 0,
            'NumberOfCycles': self.ctx.parameters['raspa_widom_cycles'],
            'PrintPropertiesEvery': self.ctx.parameters['raspa_widom_cycles'] / self.ctx.parameters['raspa_verbosity'],
            'PrintEvery': int(1e10),
            'HeliumVoidFraction': self.ctx.geom['POAV_Volume_fraction'],
        })
        param['Component'][self.inputs.cation.value].update({
            'CreateNumberOfMolecules': 0,
        })
        param['Component'][self.ctx.molecule['name']] = {
            'MoleculeDefinition': 'Local',
            'ExtraFrameworkMolecule': False,
            'WidomProbability': 1.0,
        }

        # Check particular conditions and settings
        if self.ctx.geom['Number_of_blocking_spheres'] > 0:
            param['Component'][self.ctx.molecule['name']]['BlockPocketsFileName'] = 'block_file'

        return param

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa."""

        # Update labels
        self.ctx.inp['metadata']['label'] = 'RaspaWidom'
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_widom'

        # Update parameters dict, add block file, add restart
        self.ctx.raspa_param = self._update_param_for_widom()
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
        if self.ctx.geom['Number_of_blocking_spheres'] > 0 and self.ctx.multitemp_mode != 'run_single_temp':
            self.ctx.inp['raspa']['block_pocket'] = {'block_file': self.ctx.zeopp.outputs.block}
        self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_nvt_cations[-1].outputs.retrieved

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report('Running Raspa Widom @ {}K for the Henry coefficient'.format(self.ctx.temperature))
        return ToContext(raspa_widom=running)

    def should_run_gcmc(self):
        """Output the widom results and decide to compute the isotherm if kH > kHmin, as defined by the user"""

        self.ctx.is_kh_enough = list(self.ctx.raspa_widom.outputs['output_parameters']['framework_1']['components'].
                                     values())[1]['henry_coefficient_average'] > self.ctx.parameters['raspa_minKh']

        if self.ctx.is_kh_enough:
            self.report('kH larger than the threshold: continue')
            return True

        self.report('kHh lower than the threshold: stop')
        return False

    def _update_param_for_gcmc(self):
        """Update Raspa input parameter, from Widom to GCMC"""

        param = self.ctx.raspa_param
        param['GeneralSettings'].update({
            'NumberOfInitializationCycles': self.ctx.parameters['raspa_gcmc_init_cycles'],
            'NumberOfCycles': self.ctx.parameters['raspa_gcmc_prod_cycles'],
            'PrintPropertiesEvery': int(1e6),
            'PrintEvery': self.ctx.parameters['raspa_gcmc_prod_cycles'] / self.ctx.parameters['raspa_verbosity']
        })
        param['Component'][self.ctx.molecule['name']].update({
            'WidomProbability': 0.0,
            'TranslationProbability': 1.0,
            'ReinsertionProbability': 1.0,
            'SwapProbability': 2.0,
        })
        # Check particular conditions
        if not self.ctx.molecule['singlebead']:
            param['Component'][self.ctx.molecule['name']].update({'RotationProbability': 1.0})
        return param

    def init_raspa_gcmc(self):
        """Choose the pressures we want to sample, report some details, and update settings for GCMC"""

        self.ctx.current_p_index = 0
        self.ctx.pressures = choose_pressure_points(self.ctx.parameters, self.ctx.geom,
                                                    self.ctx.raspa_widom.outputs.output_parameters)

        self.report('Computed Kh(mol/kg/Pa)={:.2e} POAV(cm3/g)={:.3f} Qsat(mol/kg)={:.2f}'.format(
            list(self.ctx.raspa_widom.outputs['output_parameters']['framework_1']['components'].values())[1]
            ['henry_coefficient_average'], self.ctx.geom['POAV_cm^3/g'], self.ctx.geom['Estimated_saturation_loading']))
        self.report('Now evaluating the isotherm @ {}K for {} pressure points'.format(
            self.ctx.temperature, len(self.ctx.pressures)))

        self.ctx.raspa_param = self._update_param_for_gcmc()

    def should_run_another_gcmc(self):
        """We run another raspa calculation only if the current iteration is
        smaller than the total number of pressures we want to compute.
        """
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_raspa_gcmc(self):
        """Run a GCMC calculation in Raspa @ T,P. """

        # Update labels
        self.ctx.inp['metadata']['label'] = 'RaspaGCMC_{}'.format(self.ctx.current_p_index + 1)
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_gcmc_{}'.format(self.ctx.current_p_index + 1)

        # Update pressure (NOTE: need to convert from bar to Pa)
        self.ctx.raspa_param['System']['framework_1']['ExternalPressure'] = \
            self.ctx.pressures[self.ctx.current_p_index] * 1e5

        # Update parameters Dict
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

        # Update restart (if present, i.e., if current_p_index>0)
        if self.ctx.current_p_index > 0:
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_gcmc[self.ctx.current_p_index -
                                                                                   1].outputs.retrieved

        # Create the calculation process, launch it and update pressure index
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report('Running Raspa GCMC @ {}K/{:.3f}bar (pressure {} of {})'.format(
            self.ctx.temperature, self.ctx.pressures[self.ctx.current_p_index], self.ctx.current_p_index + 1,
            len(self.ctx.pressures)))
        self.ctx.current_p_index += 1
        return ToContext(raspa_gcmc=append_(running))

    def return_output_parameters(self):
        """Merge all the parameters into output_parameters, depending on is_porous and is_kh_ehough."""

        gcmc_out_dict = {}
        if self.ctx.geom['is_porous'] and not self.ctx.multitemp_mode == 'run_geom_only':
            widom_out = self.ctx.raspa_widom.outputs.output_parameters
            if self.ctx.is_kh_enough:
                for calc in self.ctx.raspa_gcmc:
                    gcmc_out_dict[calc.label] = calc.outputs.output_parameters
            else:
                self.ctx.pressures = None
        else:
            widom_out = None
            self.ctx.pressures = None

        self.out(
            'output_parameters',
            get_output_parameters(geom_out=self.ctx.geom,
                                  inp_params=self.ctx.parameters,
                                  widom_out=widom_out,
                                  pressures=self.ctx.pressures,
                                  **gcmc_out_dict))

        if not self.ctx.multitemp_mode == 'run_geom_only':
            self.report('Isotherm @ {}K computed: ouput Dict<{}>'.format(self.ctx.temperature,
                                                                         self.outputs['output_parameters'].pk))
