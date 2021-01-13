# -*- coding: utf-8 -*-
"""IsothermAccurate work chain."""

import os
import functools
import ruamel.yaml as yaml

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Str, List, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, append_, while_, if_
from aiida_lsmo.utils import check_resize_unit_cell, dict_merge, validate_dict
from aiida_lsmo.utils.isotherm_molecules_schema import ISOTHERM_MOLECULES_SCHEMA
from .parameters_schemas import FF_PARAMETERS_VALIDATOR, Required, Optional, NUMBER
# import sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# import calculations
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name

# import aiida data
CifData = DataFactory('cif')  # pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


# calcfunctions (in order of appearence)
@calcfunction
def get_molecule_dict(molecule_name):
    """Get a Dict from the isotherm_molecules.yaml"""

    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, 'isotherm_data', 'isotherm_molecules.yaml')
    with open(yamlfile, 'r') as stream:
        yaml_dict = yaml.safe_load(stream)
        ISOTHERM_MOLECULES_SCHEMA(yaml_dict)
    molecule_dict = yaml_dict[molecule_name.value]
    return Dict(dict=molecule_dict)


@calcfunction
def get_atomic_radii(isotparam):
    """Get {ff_framework}.rad as SinglefileData form workchain/isotherm_data. If not existing use DEFAULT.rad."""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    filename = isotparam['ff_framework'] + '.rad'
    filepath = os.path.join(thisdir, 'isotherm_data', filename)
    if not os.path.isfile(filepath):
        filepath = os.path.join(thisdir, 'isotherm_data', 'DEFAULT.rad')
    return SinglefileData(file=filepath)


@calcfunction
def get_zeopp_parameters(molecule_dict, isotparam):
    """Get the ZeoppParameters from the inputs of the workchain"""
    probe_rad = molecule_dict['proberad'] * isotparam['zeopp_probe_scaling']
    param_dict = {
        'ha': 'DEF',
        'volpo': [probe_rad, probe_rad, isotparam['zeopp_volpo_samples']],
        'block': [probe_rad, isotparam['zeopp_block_samples']],
    }
    return ZeoppParameters(dict=param_dict)


@calcfunction
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
def get_output_parameters(geom_out, inp_params, widom_out=None, **gcmc_out_dict):
    """Merge results from all the steps of the work chain."""

    out_dict = geom_out.get_dict()

    if out_dict['is_porous'] and widom_out:
        widom_out_mol = list(widom_out['framework_1']['components'].values())[0]

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
                #'pressure': pressures, # TODO: I need to think how to pass the pressures!
                'pressure_unit': 'bar',
                'loading_absolute_average': [],
                'loading_absolute_dev': [],
                'loading_absolute_unit': 'mol/kg',
                'enthalpy_of_adsorption_average': [],
                'enthalpy_of_adsorption_dev': [],
                'enthalpy_of_adsorption_unit': 'kJ/mol'
            }

            conv_ener = 1.0 / 120.273  # K to kJ/mol
            for i, _ in range(len(gcmc_out_dict)):  # I want to keep the order imposed with the index
                gcmc_out = gcmc_out_dict['RaspaGCMC_{}'.format(i + 1)]['framework_1']
                gcmc_out_mol = list(gcmc_out['components'].values())[0]
                conv_load = gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg']

                for label in ['loading_absolute_average', 'loading_absolute_dev']:
                    isotherm[label].append(conv_load * gcmc_out_mol[label])

                for label in ['enthalpy_of_adsorption_average', 'enthalpy_of_adsorption_dev']:
                    if gcmc_out['general'][label]:
                        isotherm[label].append(conv_ener * gcmc_out['general'][label])
                    else:  # when there are no particles and Raspa return Null enthalpy
                        isotherm[label].append(None)

            # TODO: we could finally sort all the value by increasing pressure

            out_dict.update({
                'isotherm': isotherm,
                'conversion_factor_molec_uc_to_cm3stp_cm3': gcmc_out_mol['conversion_factor_molec_uc_to_cm3stp_cm3'],
                'conversion_factor_molec_uc_to_mg_g': gcmc_out_mol['conversion_factor_molec_uc_to_mg_g'],
                'conversion_factor_molec_uc_to_mol_kg': gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg'],
            })

    return Dict(dict=out_dict)


class IsothermAccurateWorkChain(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient.
    """

    parameters_schema = FF_PARAMETERS_VALIDATOR.extend({
        Required('zeopp_probe_scaling', default=1.0, description="scaling probe's diameter: molecular_rad * scaling"):
            NUMBER,
        Required('zeopp_volpo_samples',
                 default=int(1e5),
                 description='Number of samples for VOLPO calculation (per UC volume).'):
            int,
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
        Optional('temperature_list', description='To be used by IsothermMultiTempWorkChain.'):
            list,
        Required('loading_lowp_epsilon', default=0.1, description='Epsilon for convergence of low pressure loading.'):
            NUMBER,
        Required('loading_higp_sigma', default=0.95, description='Sigma fraction, to consider the sytem saturated.'):
            NUMBER,
        Required('ph0_reiteration_coeff', default=0.8, description='Coefficient for P_H0 to iterate GCMC at lower P.'):
            NUMBER,
    })
    parameters_info = parameters_schema.schema  # shorthand for printing

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])

        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        spec.input('structure', valid_type=CifData, help='Adsorbent framework CIF.')

        spec.input('molecule',
                   valid_type=(Str, Dict),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')

        spec.input('parameters',
                   valid_type=Dict,
                   validator=functools.partial(validate_dict, schema=cls.parameters_schema),
                   help='Parameters for the Isotherm workchain (see workchain.schema for default values).')

        spec.input('geometric',
                   valid_type=Dict,
                   required=False,
                   help='[Only used by IsothermMultiTempWorkChain] Already computed geometric properties')

        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            if_(cls.should_run_widom)(  # run Widom only if porous
                cls.run_raspa_widom,  # run raspa widom calculation
                if_(cls.should_run_gcmc)(  # Kh is high enough
                    cls.init_raspa_gcmc,  # initializate setting for GCMC
                    while_(cls.should_run_another_gcmc_lowp)(
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation at low pressure
                    ),
                    while_(cls.should_run_another_gcmc_highp)(
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation untill saturation
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
        @calcfunction
        def get_valid_dict(dict_node):
            return Dict(dict=self.parameters_schema(dict_node.get_dict()))

        self.ctx.parameters = get_valid_dict(self.inputs.parameters)

        # Get integer temperature in context for easy reports
        self.ctx.temperature = int(round(self.ctx.parameters['temperature']))

        # Understand if IsothermMultiTempWorkChain is calling this work chain
        if 'geometric' in self.inputs:
            self.ctx.multitemp_mode = 'run_single_temp'
        elif 'temperature_list' in self.ctx.parameters.attributes:
            self.ctx.multitemp_mode = 'run_geom_only'
        else:
            self.ctx.multitemp_mode = None

    def run_zeopp(self):
        """Step 1: run zeo++ to calculate the density, void fraction and POAV
        and check if blocking spheres are needed or not.
        """

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
                'structure': self.inputs.structure,
                'atomic_radii': get_atomic_radii(self.ctx.parameters),
                'parameters': get_zeopp_parameters(self.ctx.molecule, self.ctx.parameters)
            })

        running = self.submit(ZeoppCalculation, **inputs)
        self.report('Running zeo++ block and volpo for {} Calculation<{}>'.format(self.ctx.molecule['name'],
                                                                                  running.id))
        return ToContext(zeopp=running)

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume,
        also check the number of blocking spheres and estimate the saturation loading.
        Also, stop if called by IsothermMultiTemp for geometric results only.

        Step 2: Calculate the theoretical q_sat based on the liquid density of the molecule (in the function
        get_geometric dict).
        """

        # Get geometric properties and consider if IsothermMultiTempWorkChain is calling this workchain
        if self.ctx.multitemp_mode == 'run_single_temp':
            self.ctx.geom = self.inputs.geometric
            return True
        self.ctx.geom = get_geometric_dict(self.ctx.zeopp.outputs.output_parameters, self.ctx.molecule)

        if self.ctx.geom['is_porous']:
            self.report('Found accessible pore volume for {}: continue'.format(self.ctx.molecule['name']))
            self.report('Found {} blocking spheres'.format(self.ctx.geom['Number_of_blocking_spheres']))
            # Return block file only if blocking spheres are present
            if self.ctx.geom['Number_of_blocking_spheres'] > 0:
                self.out_many(self.exposed_outputs(self.ctx.zeopp, ZeoppCalculation))
        else:
            self.report('No accessible pore volume to {}: stop'.format(self.ctx.molecule['name']))

        return self.ctx.geom['is_porous'] and not self.ctx.multitemp_mode == 'run_geom_only'

    def _get_widom_param(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        param = {
            'GeneralSettings': {
                'SimulationType':
                    'MonteCarlo',
                'NumberOfInitializationCycles':
                    0,
                'NumberOfCycles':
                    self.ctx.parameters['raspa_widom_cycles'],
                'PrintPropertiesEvery':
                    self.ctx.parameters['raspa_widom_cycles'] / self.ctx.parameters['raspa_verbosity'],
                'PrintEvery':
                    int(1e10),
                'RemoveAtomNumberCodeFromLabel':
                    True,  # BE CAREFULL: needed in AiiDA-1.0.0 because of github.com/aiidateam/aiida-core/issues/3304
                'Forcefield':
                    'Local',
                'UseChargesFromCIFFile':
                    'yes',
                'CutOff':
                    self.ctx.parameters['ff_cutoff'],
            },
            'System': {
                'framework_1': {
                    'type': 'Framework',
                    'HeliumVoidFraction': self.ctx.geom['POAV_Volume_fraction'],
                    'ExternalTemperature': self.ctx.parameters['temperature'],
                }
            },
            'Component': {
                self.ctx.molecule['name']: {
                    'MoleculeDefinition': 'Local',
                    'WidomProbability': 1.0,
                },
            },
        }

        # Check particular conditions and settings
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
        param['System']['framework_1']['UnitCells'] = '{} {} {}'.format(mult[0], mult[1], mult[2])

        if self.ctx.geom['Number_of_blocking_spheres'] > 0:
            param['Component'][self.ctx.molecule['name']]['BlockPocketsFileName'] = 'block_file'

        if self.ctx.molecule['charged']:  # NOTE: `Chargemethod Ewald` is the default in Raspa!
            param['GeneralSettings'].update({'ChargeMethod': 'Ewald', 'EwaldPrecision': 1e-6})
        else:
            param['GeneralSettings'].update({'ChargeMethod': 'None'})

        if 'rosenbluth' in self.ctx.molecule.keys():  # flexible molecule which need a correction for the chem pot
            param['Component'][self.ctx.molecule['name']]['IdealGasRosenbluthWeight'] = self.ctx.molecule['rosenbluth']

        return param

    def run_raspa_widom(self):
        """Step 3.0: Run the “henry” function to obtain KH and PH."""

        # Initialize the input for raspa_base, which later will need only minor updates for GCMC
        self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.inp['metadata']['label'] = 'RaspaWidom'
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_widom'

        self.ctx.inp['raspa']['framework'] = {'framework_1': self.inputs.structure}
        if self.ctx.geom['Number_of_blocking_spheres'] > 0 and self.ctx.multitemp_mode != 'run_single_temp':
            self.ctx.inp['raspa']['block_pocket'] = {'block_file': self.ctx.zeopp.outputs.block}

        self.ctx.raspa_param = self._get_widom_param()
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

        # Generate the force field with the ff_builder
        ff_params = get_ff_parameters(self.ctx.molecule, self.ctx.parameters)

        files_dict = FFBuilder(ff_params)
        self.ctx.inp['raspa']['file'] = files_dict

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report('Running Raspa Widom {} @ {}K for the Henry coefficient'.format(self.ctx.molecule['name'],
                                                                                    self.ctx.temperature))

        return ToContext(raspa_widom=running)

    def should_run_gcmc(self):
        """Output the widom results and decide to compute the isotherm if kH > kHmin, as defined by the user.
        """

        self.ctx.kh = list(self.ctx.raspa_widom.outputs['output_parameters']['framework_1']
                           ['components'].values())[0]['henry_coefficient_average']
        self.ctx.is_kh_enough = self.ctx.kh > self.ctx.parameters['raspa_minKh']

        if self.ctx.is_kh_enough:
            self.report('kH larger than the threshold for {}: continue'.format(self.ctx.molecule['name']))
            return True

        self.report('kHh lower than the threshold for {}: stop'.format(self.ctx.molecule['name']))
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

        if 'rosenbluth' in self.ctx.molecule.keys():  # Flexible molecule needs ConfigurationalBias move
            param['Component'][self.ctx.molecule['name']].update({'CBMCProbability': 1.0})

        return param

    def init_raspa_gcmc(self):
        """Update settings for GCMC.
        """
        self.ctx.raspa_param = self._update_param_for_gcmc()
        self.ctx.first_iteration = True
        self.ctx.gcmc_i = 0

    def _get_last_loading_molkg(self):
        """Get the GCMC loading in mol/kg from the last GCMC calculation.
        """
        gcmc_out = self.ctx.raspa_gcmc[-1]['framework_1']
        gcmc_out_mol = list(gcmc_out['components'].values())[0]
        conv_load = gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg']
        load_molec_uc = gcmc_out_mol['loading_absolute_average']
        return conv_load * load_molec_uc

    def should_run_another_gcmc_lowp(self):
        """Step 3.1: I calculate the initial guess PH0
        Step 3.2: To check if P_H0 belongs to the Henry’s regime, the error between the obtained uptake and Henry’s
        uptake should be smaller than a precision value, epsilon. If the error is higher than epsilon,
        P_HO is multiplied by a factor of 0.8 and the second step is repeated until the error converges to a value
        smaller than epsilon.
        """
        if self.ctx.first_iteration:
            self.ctx.first_iteration = False
            self.ctx.ph0 = 10000 / self.ctx.geom['POAV_A^3'] / self.ctx.geom[
                'Density'] / 6.02 / self.ctx.kh  #TODO: check units!
            self.ctx.pressure = self.ctx.ph0
            return True
        else:
            self.ctx.gcmc_loading_average = self._get_last_loading_molkg()
            loading_convergence = (gcmc_loading_average - self.ctx.kh * self.ctx.pressure) / gcmc_loading_average
            if loading_convergence < self.ctx.parameters['loading_lowp_epsilon']:
                self.ctx.first_iteration = True
                return False  # Converged!
            else:
                self.ctx.pressure *= self.ctx.parameters['ph0_reiteration_coeff']
                self.ctx.gcmc_i += 1
                return True  # Not converged: lower P_H0 and recompute GCMC

        return self.ctx.current_p_index < len(self.ctx.pressures)

    def should_run_another_gcmc_higp(self):
        """Step 4: determine the pressure at which the saturation starts, Psat
        Step 5: After calculating PH and Psat, pressure values in between are generated based on the following
        sampling scheme: more pressure points are needed when the isotherm is steep compared to when the isotherm
        is closer to saturation. The maximum pressure step is determined so that in the smoothest case 20 points are
        generated.
        """
        if self.ctx.first_iteration:  # Proceed with first iteration
            self.ctx.first_iteration = False
            self.ctx.psat = self.ctx.parameters['loading_sat_sigma'] * self.ctx.geom['Estimated_saturation_loading'] /\
                                (1-self.ctx.parameters['loading_sat_sigma'])/self.ctx.kh  #TODO: check units!
            self.ctx.dpmax = (self.ctx.psat - self.ctx.ph0) / 20  # Comment: maybe more efficient to use the first ph0
            self.ctx.pressure_old = self.ctx.ph0
            self.ctx.pressure = self.ctx.pressure_old + self.ctx.dpmax
            return True
        else:
            gcmc_loading_average_old = self.ctx.gcmc_loading_average
            self.ctx.gcmc_loading_average = self._get_last_loading_molkg()
            loading_derivative = (self.ctx.gcmc_loading_average - gcmc_loading_average_old) / (self.ctx.pressure -
                                                                                               self.ctx.pressure_old)
            if loading_derivative < 0:  #Need to recompute GCMC. NOTE: very dangerous, as it can lead to infinite loop!
                del self.ctx.raspa_gcmc[-1]
                return True
            dp_computed = 2 * self.ctx.pressure / self.ctx.psat * self.ctx.geom[
                'Estimated_saturation_loading'] / loading_derivative
            self.ctx.pressure_old = self.ctx.pressure
            self.ctx.pressure = self.ctx.pressure_old + min(self.ctx.dpmax, dp_computed)
            if self.ctx.pressure > self.ctx.psat:
                return False  # Reached saturation!
            else:
                self.ctx.gcmc_i += 1
                return True  # Not converged: lower P_H0 and recompute GCMC

        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_raspa_gcmc(self):
        """Run a GCMC calculation in Raspa @ T,P."""

        # Update labels
        self.ctx.inp['metadata']['label'] = 'RaspaGCMC_{}'.format(self.ctx.gcmc_i)
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_gcmc_{}'.format(self.ctx.gcmc_i)

        # Update pressure (NOTE: need to convert from bar to Pa)
        self.ctx.raspa_param['System']['framework_1']['ExternalPressure'] = self.ctx.pressure * 1e5

        # Update parameters Dict
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

        # Update restart (if present, i.e., if current_p_index>0)
        if 'raspa_gcmc' in self.ctx:
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_gcmc[-1].outputs.retrieved

        # Create the calculation process, launch it and update pressure index
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report('Running Raspa GCMC {} @ {}K/{:.3f}bar'.format(self.ctx.molecule['name'], self.ctx.temperature,
                                                                   self.ctx.pressure))
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
                                  **gcmc_out_dict))

        if not self.ctx.multitemp_mode == 'run_geom_only':
            self.report('Isotherm {} @ {}K computed: ouput Dict<{}>'.format(self.ctx.molecule['name'],
                                                                            self.ctx.temperature,
                                                                            self.outputs['output_parameters'].pk))
