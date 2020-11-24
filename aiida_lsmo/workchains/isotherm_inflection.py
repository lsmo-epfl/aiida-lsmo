# -*- coding: utf-8 -*-
"""A work chain."""
import functools
import numpy as np

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Str, List
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, if_
from aiida_lsmo.utils import check_resize_unit_cell, dict_merge

from .isotherm import (get_molecule_dict, get_atomic_radii, get_ff_parameters, get_zeopp_parameters, validate_dict,
                       get_geometric_dict)
from .parameters_schemas import FF_PARAMETERS_VALIDATOR, NUMBER, Required, Optional

RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name
CifData = DataFactory('cif')  # pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name

NAVO = 6.022E+23  # molec/mol
A3TOL = 1E-27  # L/A3


# calcfunctions (in order of appearence)
@calcfunction
def get_pressure_points(molecule_dict, isotparam):
    """Multiply p/p0 with p0 to have pressure points in bar if pressure_list!=None,
    or choose them based on pressure_min/max/num, to be equispaced in a Log plot.
    """
    if 'pressure_list' in isotparam.attributes:
        pressure_points = [x * molecule_dict['pressure_zero'] for x in isotparam['pressure_list']]
    else:  # pressure_list==None
        exp_min = np.log10(isotparam['pressure_min'] * molecule_dict['pressure_zero'])
        exp_max = np.log10(isotparam['pressure_max'] * molecule_dict['pressure_zero'])
        exp_list = np.linspace(exp_min, exp_max, isotparam['pressure_num'])
        pressure_points = [10**x for x in exp_list]
    return List(list=pressure_points)


@calcfunction
def get_output_parameters(inp_params, pressures, geom_out, widom_out, **gcmc_dict):  # pylint: disable=too-many-locals
    """Merge results from all the steps of the work chain.
    geom_out (Dict) contains the output of Zeo++
    widom_out (Dict) contains the output of Raspa's Widom insertions calculation
    gcmc_dict (dict of Dicts) has the keys like: inp/out_RaspaGCMC/RaspaGCMCNew/RaspaGCMCSat_1..n
    """

    out_dict = geom_out.get_dict()

    if out_dict['is_porous'] and widom_out:
        widom_out_mol = list(widom_out['framework_1']['components'].values())[0]

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

        gcmc_out_mol = list(list(
            gcmc_dict.values())[0]['framework_1']['components'].values())[0]  #same for all, take just the first!

        isotherm = {
            'temperature': inp_params['temperature'],
            'temperature_unit': 'K',
            'pressure': pressures,
            'pressure_unit': 'bar',
            'loading_absolute_average_from_dil': [],
            'loading_absolute_average_from_sat': [],
            'loading_absolute_dev_from_dil': [],
            'loading_absolute_dev_from_sat': [],
            'loading_absolute_unit': 'mol/kg',
            'enthalpy_of_adsorption_average_from_dil': [],
            'enthalpy_of_adsorption_average_from_sat': [],
            'enthalpy_of_adsorption_dev_from_dil': [],
            'enthalpy_of_adsorption_dev_from_sat': [],
            'enthalpy_of_adsorption_unit': 'kJ/mol',
            'conversion_factor_molec_uc_to_cm3stp_cm3': gcmc_out_mol['conversion_factor_molec_uc_to_cm3stp_cm3'],
            'conversion_factor_molec_uc_to_mg_g': gcmc_out_mol['conversion_factor_molec_uc_to_mg_g'],
            'conversion_factor_molec_uc_to_mol_kg': gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg'],
        }

        conv_ener = 1.0 / 120.273  # K to kJ/mol
        conv_load = isotherm['conversion_factor_molec_uc_to_mol_kg']

        for i in range(len(pressures)):
            for inp, key_ext in [('dil', '_from_dil'), ('sat', '_from_sat')]:
                gcmc_out = gcmc_dict['raspa_gcmc_{}_{}'.format(inp, i + 1)]['framework_1']
                gcmc_out_mol = list(gcmc_out['components'].values())[0]

                for label in ['loading_absolute_average', 'loading_absolute_dev']:
                    isotherm[label + key_ext].append(conv_load * gcmc_out_mol[label])

                for label in ['enthalpy_of_adsorption_average', 'enthalpy_of_adsorption_dev']:
                    if gcmc_out['general'][label]:
                        isotherm[label + key_ext].append(conv_ener * gcmc_out['general'][label])
                    else:  # when there are no particles and Raspa return Null enthalpy
                        isotherm[label + key_ext].append(None)

        out_dict['isotherm'] = isotherm

    return Dict(dict=out_dict)


class IsothermInflectionWorkChain(WorkChain):
    """A work chain to compute single component isotherms at adsorption and desorption:
    GCMC calculations are run in parallell at all pressures, starting from the empty framework and the saturated system.
    This workchain is useful to spot adsorption hysteresis.
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
        Required('temperature', default=300, description='Temperature of the simulation.'):
            NUMBER,
        Required('pressure_min', default=0.001, description='Lower pressure to sample (bar).'):
            NUMBER,  # TODO: MIN selected from the henry coefficient!  # pylint: disable=fixme
        Required('pressure_max', default=1.0, description='Upper pressure to sample (bar).'):
            NUMBER,
        Optional('pressure_list',
                 description='Pressure list for the isotherm (bar): if given it will skip to guess it.'):
            list,
        Required('pressure_num',
                 default=20,
                 description='Number of pressure points considered, eqispaced in a log plot'):
            int,
        Required('box_length', description='length of simulation box for simulation without framework'):
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

        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            if_(cls.should_run_widom)(  # if is porous
                cls.run_raspa_widom_and_sat,  # run raspa (1) widom calculation (2) NVT at qsat
                cls.run_raspa_gcmc_from_dil_sat,  # run all calc from zero uptake and sat uptake
            ),
            cls.return_output_parameters,
        )

        spec.expose_outputs(ZeoppCalculation, include=['block'])  #only if porous

        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='Results of the single temperature wc: keys can vay depending on is_porous.')

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

        # Keeping the infrastructure for multitemperature, but excluding it
        self.ctx.multitemp_mode = None

    def run_zeopp(self):
        """Perform Zeo++ block and VOLPO calculations."""
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
        """

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
        self.ctx.number_uc = mult[0] * mult[1] * mult[2]

        if self.ctx.geom['Number_of_blocking_spheres'] > 0:
            param['Component'][self.ctx.molecule['name']]['BlockPocketsFileName'] = 'block_file'

        if self.ctx.molecule['charged']:  # NOTE: `Chargemethod Ewald` is the default in Raspa!
            param['GeneralSettings'].update({'ChargeMethod': 'Ewald', 'EwaldPrecision': 1e-6})
        else:
            param['GeneralSettings'].update({'ChargeMethod': 'None'})

        if 'rosenbluth' in self.ctx.molecule.keys():  # flexible molecule which need a correction for the chem pot
            param['Component'][self.ctx.molecule['name']]['IdealGasRosenbluthWeight'] = self.ctx.molecule['rosenbluth']

        return param

    def _update_param_for_gcmc(self, number_of_molecules=0, swap_prob=0.5):
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
            'ReinsertionProbability': 0.5,
            'SwapProbability': swap_prob,
            'CreateNumberOfMolecules': number_of_molecules,
        })
        # Check particular conditions
        if not self.ctx.molecule['singlebead']:
            param['Component'][self.ctx.molecule['name']].update({'RotationProbability': 1.0})

        if 'rosenbluth' in self.ctx.molecule.keys():  # Flexible molecule needs ConfigurationalBias move
            param['Component'][self.ctx.molecule['name']].update({'CBMCProbability': 1.0})

        return param

    def _get_saturation_molecules(self):
        """Compute the estimate of molecules at saturation by: pore_vol * lid_dens * number_uc."""
        molar_dens = self.ctx.molecule['molsatdens']  #mol/l
        pore_vol = self.ctx.zeopp.outputs.output_parameters['POAV_A^3']  #A3/UC
        molec_uc = molar_dens * NAVO * A3TOL * pore_vol  # molec/UC
        molec_tot = molec_uc * self.ctx.number_uc
        return int(molec_tot)

    def _get_mid_dens_molecules(self, raspa_calc_dil, raspa_calc_sat):
        """Given a calculation at diluted and saturation condition, compute the total molecules at mid density."""
        molec_uc_dil = list(raspa_calc_dil.outputs.output_parameters['framework_1']
                            ['components'].values())[0]['loading_absolute_average']
        molec_uc_sat = list(raspa_calc_sat.outputs.output_parameters['framework_1']
                            ['components'].values())[0]['loading_absolute_average']
        molec_uc_mid = max(0, (molec_uc_sat + molec_uc_dil) / 2)
        molec_tot_mid = molec_uc_mid * self.ctx.number_uc
        return int(molec_tot_mid)

    def run_raspa_widom_and_sat(self):
        """Run a Widom calculation in Raspa."""

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

        self.to_context(raspa_widom=running)

        # Run Raspa NVT at 90% Qsat
        self.ctx.inp['metadata']['label'] = 'RaspaNVT_sat'
        self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_nvt_sat'

        nmolecules = int(0.90 * self._get_saturation_molecules())
        self.ctx.raspa_param = self._update_param_for_gcmc(number_of_molecules=nmolecules, swap_prob=0.01)
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)

        self.report('Running Raspa NVT with {} molecules'.format(nmolecules))
        self.to_context(raspa_nvt_sat=running)

    def run_raspa_gcmc_from_dil_sat(self):
        """For each pressure point, run GCMC calculation from both diluted and saturated initial conditions."""

        self.report('Computed {} {}'.format(
            self.ctx.raspa_widom.outputs.output_parameters.pk,
            self.ctx.raspa_nvt_sat.outputs.output_parameters.pk,
        ))

        self.ctx.pressures = get_pressure_points(self.ctx.molecule, self.ctx.parameters)
        self.ctx.raspa_param = self._update_param_for_gcmc()

        for inp in ['dil', 'sat']:
            for i, _ in enumerate(self.ctx.pressures):
                # Update labels
                self.ctx.inp['metadata']['label'] = 'RaspaGCMC_{}_{}'.format(inp, i + 1)
                self.ctx.inp['metadata']['call_link_label'] = 'run_raspa_gcmc_{}_{}'.format(inp, i + 1)

                self.ctx.raspa_param['System']['framework_1']['ExternalPressure'] = self.ctx.pressures[i] * 1e5

                # Update parameters Dict
                self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

                if inp == 'dil':
                    pass
                elif inp == 'sat':
                    self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_nvt_sat.outputs.retrieved

                # Create the calculation process, launch it and update pressure index
                running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
                self.report('Running Raspa GCMC {} {}'.format(inp, i + 1))

                key_ctx = 'raspa_gcmc_{}_{}'.format(inp, i + 1)
                self.to_context(**{key_ctx: running})

    def return_output_parameters(self):
        """Merge all the parameters into output_parameters, depending on is_porous and is_kh_ehough."""

        gcmc_dict = {}
        if self.ctx.geom['is_porous'] and not self.ctx.multitemp_mode == 'run_geom_only':
            widom_out = self.ctx.raspa_widom.outputs.output_parameters
            for key in self.ctx.keys():
                if 'raspa_gcmc' in key:
                    gcmc_dict[key] = self.ctx[key].outputs.output_parameters
        else:
            widom_out = None

        self.out(
            'output_parameters',
            get_output_parameters(inp_params=self.ctx.parameters,
                                  pressures=self.ctx.pressures,
                                  geom_out=self.ctx.geom,
                                  widom_out=widom_out,
                                  **gcmc_dict))

        if not self.ctx.multitemp_mode == 'run_geom_only':
            self.report('IsothermInflection {} @ {}K computed: ouput Dict<{}>'.format(
                self.ctx.molecule['name'], self.ctx.temperature, self.outputs['output_parameters'].pk))
