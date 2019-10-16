# -*- coding: utf-8 -*-
"""Isotherm workchain"""
from __future__ import absolute_import

import os
from six.moves import range

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Str, List, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, append_, while_, if_
from aiida_lsmo.utils import check_resize_unit_cell, aiida_dict_merge

# sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# calculation objects
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name

# data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


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
def get_atomic_radii(isotparam):
    """Get {forcefield}.rad as SinglefileData form workchain/isotherm_data"""
    forcefield = isotparam['forcefield']
    thisdir = os.path.dirname(os.path.abspath(__file__))
    fullfilename = forcefield + ".rad"
    return SinglefileData(file=os.path.join(thisdir, "isotherm_data", fullfilename))


@calcfunction
def get_zeopp_parameters(molecule_dict, isotparam):
    """Get the ZeoppParameters from the inputs of the workchain"""
    probe_rad = molecule_dict["proberad"]
    param_dict = {
        'ha': 'DEF',
        'volpo': [probe_rad, probe_rad, isotparam['zeopp_volpo_samples']],
        'block': [probe_rad, isotparam['zeopp_block_samples']],
    }
    return ZeoppParameters(dict=param_dict)


@calcfunction
def choose_pressure_points(inp_param, geom, raspa_widom_out):
    """If 'presure_list' is not provide, model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm, in a List.
    """
    if inp_param["pressure_list"]:
        pressure_points = inp_param["pressure_list"]
    else:
        khenry = list(raspa_widom_out["framework_1"]["components"].values())[0]['henry_coefficient_average']  #mol/kg/Pa
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
def get_geometric_output(zeopp_out, molecule):
    """Return the geometric_output Dict from Zeopp results, including Qsat and is_porous"""
    geometric_output = zeopp_out.get_dict()
    geometric_output.update({
        'Estimated_saturation_loading': zeopp_out['POAV_cm^3/g'] * molecule['molsatdens'],
        'Estimated_saturation_loading_unit': 'mol/kg',
        'is_porous': geometric_output["POAV_A^3"] > 0.000
    })
    return Dict(dict=geometric_output)


@calcfunction
def get_isotherm_output(parameters, widom_out, pressures, **gcmc_out_dict):
    """ Extract Widom and GCMC results to isotherm Dict """
    widom_out_mol = list(widom_out["framework_1"]["components"].values())[0]

    isotherm_output = {
        'temperature': parameters['temperature'],
        'temperature_unit': 'K',
        'is_kh_enough': widom_out_mol['henry_coefficient_average'] > parameters['raspa_minKh']
    }

    widom_labels = [
        'henry_coefficient_average',
        'henry_coefficient_dev',
        'henry_coefficient_unit',
        'adsorption_energy_widom_average',
        'adsorption_energy_widom_dev',
        'adsorption_energy_widom_unit',
    ]

    for label in widom_labels:
        isotherm_output.update({label: widom_out_mol[label]})

    if isotherm_output['is_kh_enough']:

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
            gcmc_out = gcmc_out_dict['RaspaGCMC_{}'.format(i + 1)]["framework_1"]
            gcmc_out_mol = list(gcmc_out["components"].values())[0]
            conv_load = gcmc_out_mol["conversion_factor_molec_uc_to_mol_kg"]

            for label in ['loading_absolute_average', 'loading_absolute_dev']:
                isotherm[label].append(conv_load * gcmc_out_mol[label])

            for label in ['enthalpy_of_adsorption_average', 'enthalpy_of_adsorption_dev']:
                if gcmc_out['general'][label]:
                    isotherm[label].append(conv_ener * gcmc_out['general'][label])
                else:  # when there are no particles and Raspa return Null enthalpy
                    isotherm[label].append(None)

        isotherm_output.update({
            "isotherm": isotherm,
            'conversion_factor_molec_uc_to_cm3stp_cm3': gcmc_out_mol['conversion_factor_molec_uc_to_cm3stp_cm3'],
            'conversion_factor_molec_uc_to_gr_gr': gcmc_out_mol['conversion_factor_molec_uc_to_gr_gr'],
            'conversion_factor_molec_uc_to_mol_kg': gcmc_out_mol['conversion_factor_molec_uc_to_mol_kg'],
        })

    return Dict(dict=isotherm_output)


# Deafault parameters
ISOTHERMPARAMETERS_DEFAULT = Dict(
    dict={  #TODO: create IsothermParameters instead of Dict # pylint: disable=fixme
        "forcefield": "UFF",  # str, Forcefield of the structure
        "ff_tailcorr": True,  # bool, Apply tail corrections
        "ff_shift": False,  # bool, Shift or truncate at cutoff
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "temperature": 300,  # float, Temperature of the simulation
        "temperature_list": None,  # list, to be used by IsothermMultiTempWorkChain
        "zeopp_volpo_samples": int(1e5),  # int, Number of samples for VOLPO calculation (per UC volume)
        "zeopp_block_samples": int(100),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_minKh": 1e-10,  # float, If Henry coefiicient < raspa_minKh do not run the isotherm (mol/kg/Pa)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_widom_cycles": int(1e5),  # int, Number of widom cycles
        "raspa_gcmc_init_cycles": int(1e3),  # int, Number of GCMC initialization cycles
        "raspa_gcmc_prod_cycles": int(1e4),  # int, Number of GCMC production cycles
        "pressure_list": None,  # list, Pressure list for the isotherm (bar): if given it will skip  guess
        "pressure_precision": 0.1,  # float, Precision in the sampling of the isotherm: 0.1 ok, 0.05 better for high res
        "pressure_maxstep": 5,  # float, Max distance between pressure points (bar)
        "pressure_min": 0.001,  # float, Lower pressure to sample (bar)
        "pressure_max": 10  # float, upper pressure to sample (bar)
    })


class IsothermWorkChain(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient.
    """

    @classmethod
    def define(cls, spec):
        super(IsothermWorkChain, cls).define(spec)

        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])

        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])

        spec.input('structure', valid_type=CifData, help='Adsorbent framework CIF')

        spec.input("molecule",
                   valid_type=(Str, Dict),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')

        spec.input("parameters",
                   valid_type=Dict,
                   help='Parameters for the Isotherm workchain. See IsothermParameters_defaults for types and default')

        spec.input("geometric", valid_type=Dict, required=False, help='Already computed Geometric properties')

        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            if_(cls.should_run_widom)(  # run Widom only if porous
                cls.run_raspa_widom,  # run raspa widom calculation
                if_(cls.should_run_gcmc)(  # Kh is high enough
                    cls.init_raspa_gcmc,  # initializate setting for GCMC
                    while_(cls.should_run_another_gcmc)(  # new pressure
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation
                    ),
                    cls.return_isotherm,
                ),
            ),
        )

        spec.expose_outputs(ZeoppCalculation, include=['block'])  #only if porous

        spec.output(
            'geometric_output',
            valid_type=Dict,
            required=False,  # only if not skip_zeopp
            help='Results of the Zeo++ calculation (density, pore volume, etc.) plus some extra results (Qsat)')

        spec.output(
            'isotherm_output',
            valid_type=Dict,
            required=False,  # only if is_porous
            help='Results of the widom calculation and (if is_kh_enough) isotherm at a single temperature')

    def setup(self):
        """Initialize the parameters"""

        # Get the molecule Dict from the yaml or directly as an input
        if isinstance(self.inputs.molecule, Str):
            self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
        elif isinstance(self.inputs.molecule, Dict):
            self.ctx.molecule = self.inputs.molecule

        # Get the parameters Dict, merging defaults with user settings
        self.ctx.parameters = aiida_dict_merge(ISOTHERMPARAMETERS_DEFAULT, self.inputs.parameters)

        # Get integer temperature in context for easy reports
        self.ctx.temperature = int(round(self.ctx.parameters['temperature']))

    def run_zeopp(self):
        """Perform Zeo++ block and VOLPO calculations."""

        # Skip zeopp calculation if the geometric properties are already provided (IsothermMultiTemp)
        self.ctx.skip_zeopp = "geometric" in self.inputs
        if self.ctx.skip_zeopp:
            return None

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')

        # Set inputs for zeopp
        inputs.update({
            'metadata': {
                'label': "ZeoppVolpoBlock",
                'call_link_label': 'run_zeopp_block_and_volpo',
            },
            'structure': self.inputs.structure,
            'atomic_radii': get_atomic_radii(self.ctx.parameters),
            'parameters': get_zeopp_parameters(self.ctx.molecule, self.ctx.parameters)
        })

        running = self.submit(ZeoppCalculation, **inputs)
        self.report("Running zeo++ block and volpo Calculation<{}>".format(running.id))
        return ToContext(zeopp=running)

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume,
        also check the number of blocking spheres and estimate the saturation loading.
        Also, stop if called by IsothermMultiTemp for geometric results only."""

        # Use geometric properties if provided in the input (from IsothermMultiTemp)
        if self.ctx.skip_zeopp:
            self.ctx.geom = self.inputs.geometric
            return True

        self.out("geometric_output", get_geometric_output(self.ctx.zeopp.outputs.output_parameters, self.ctx.molecule))
        self.ctx.geom = self.outputs['geometric_output']

        if self.ctx.geom['is_porous']:
            self.report("Found accessible pore volume: continue")
            self.report("Found {} blocking spheres".format(self.ctx.geom['Number_of_blocking_spheres']))
            # Return block file only if blocking spheres are present
            if self.ctx.geom['Number_of_blocking_spheres'] > 0:
                self.out_many(self.exposed_outputs(self.ctx.zeopp, ZeoppCalculation))
        else:
            self.report("No accessible pore volume: stop")

        return self.ctx.geom['is_porous'] and not self.ctx.parameters['temperature_list']

    def _get_widom_param(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        param = {
            "GeneralSettings": {
                "SimulationType":
                    "MonteCarlo",
                "NumberOfInitializationCycles":
                    0,
                "NumberOfCycles":
                    self.ctx.parameters['raspa_widom_cycles'],
                "PrintPropertiesEvery":
                    self.ctx.parameters['raspa_widom_cycles'] / self.ctx.parameters['raspa_verbosity'],
                "PrintEvery":
                    int(1e10),
                "RemoveAtomNumberCodeFromLabel":
                    True,  # be careful!
                "Forcefield":
                    "{}_{}_{}_{}".format(self.ctx.parameters['forcefield'], self.ctx.molecule["forcefield"],
                                         ["notc", "tc"][self.ctx.parameters['ff_tailcorr']],
                                         ["trunc", "shift"][self.ctx.parameters['ff_shift']]),
                "UseChargesFromCIFFile":
                    "yes",
                "CutOff":
                    self.ctx.parameters['ff_cutoff'],
            },
            "System": {
                "framework_1": {
                    "type": "Framework",
                    "HeliumVoidFraction": self.ctx.geom["POAV_Volume_fraction"],
                    "ExternalTemperature": self.ctx.parameters['temperature'],
                }
            },
            "Component": {
                self.ctx.molecule['name']: {
                    "MoleculeDefinition": self.ctx.molecule["forcefield"],
                    "WidomProbability": 1.0,
                },
            },
        }

        # Check particular conditions and settings
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
        param["System"]["framework_1"]["UnitCells"] = "{} {} {}".format(mult[0], mult[1], mult[2])

        if self.ctx.geom['Number_of_blocking_spheres'] > 0:
            param["Component"][self.ctx.molecule['name']]["BlockPocketsFileName"] = "block_file"

        if self.ctx.molecule['charged']:
            param["GeneralSettings"].update({"ChargeMethod": "Ewald", "EwaldPrecision": 1e-6})
        return param

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa."""

        # Initialize the input for raspa_base, which later will need only minor updates for GCMC
        self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.inp['metadata']['label'] = "RaspaWidom"
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_widom"

        self.ctx.inp['raspa']['framework'] = {"framework_1": self.inputs.structure}
        if self.ctx.geom['Number_of_blocking_spheres'] > 0:
            self.ctx.inp["raspa"]["block_pocket"] = {"block_file": self.ctx.zeopp.outputs.block}

        self.ctx.raspa_param = self._get_widom_param()
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param).store()

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa Widom @ {}K for the Henry coefficient".format(self.ctx.temperature))

        return ToContext(raspa_widom=running)

    def should_run_gcmc(self):
        """Output the widom results and decide to compute the isotherm if kH > kHmin, as defined by the user"""

        self.ctx.is_kh_enough = list(self.ctx.raspa_widom.outputs['output_parameters']["framework_1"]["components"].
                                     values())[0]['henry_coefficient_average'] > self.ctx.parameters['raspa_minKh']

        if self.ctx.is_kh_enough:
            self.report("kH larger than the threshold: continue")
            return True

        self.report("kHh lower than the threshold: stop")
        return False

    def _update_param_for_gcmc(self):
        """Update Raspa input parameter, from Widom to GCMC"""

        param = self.ctx.raspa_param
        param["GeneralSettings"].update({
            "NumberOfInitializationCycles": self.ctx.parameters['raspa_gcmc_init_cycles'],
            "NumberOfCycles": self.ctx.parameters['raspa_gcmc_prod_cycles'],
            "PrintPropertiesEvery": int(1e6),
            "PrintEvery": self.ctx.parameters['raspa_gcmc_prod_cycles'] / self.ctx.parameters['raspa_verbosity']
        })
        param["Component"][self.ctx.molecule['name']].update({
            "WidomProbability": 0.0,
            "TranslationProbability": 1.0,
            "ReinsertionProbability": 1.0,
            "SwapProbability": 2.0,
        })
        # Check particular conditions
        if not self.ctx.molecule['singlebead']:
            param["Component"][self.ctx.molecule['name']].update({"RotationProbability": 1.0})
        return param

    def init_raspa_gcmc(self):
        """Choose the pressures we want to sample, report some details, and update settings for GCMC"""

        self.ctx.current_p_index = 0
        self.ctx.pressures = choose_pressure_points(self.ctx.parameters, self.ctx.geom,
                                                    self.ctx.raspa_widom.outputs.output_parameters)

        self.report("Computed Kh(mol/kg/Pa)={:.2e} POAV(cm3/g)={:.3f} Qsat(mol/kg)={:.2f}".format(
            list(self.ctx.raspa_widom.outputs['output_parameters']["framework_1"]["components"].values())[0]
            ['henry_coefficient_average'], self.ctx.geom['POAV_cm^3/g'], self.ctx.geom['Estimated_saturation_loading']))
        self.report("Now evaluating the isotherm @ {}K for {} pressure points".format(
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
        self.ctx.inp['metadata']['label'] = "RaspaGCMC_{}".format(self.ctx.current_p_index + 1)
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_gcmc_{}".format(self.ctx.current_p_index + 1)

        # Update pressure (NOTE: need to convert from bar to Pa)
        self.ctx.raspa_param["System"]["framework_1"]['ExternalPressure'] = \
            self.ctx.pressures[self.ctx.current_p_index] * 1e5

        # Update parameters Dict
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param).store()

        # Update restart (if present, i.e., if current_p_index>0)
        if self.ctx.current_p_index > 0:
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_gcmc[self.ctx.current_p_index -
                                                                                   1].outputs.retrieved

        # Create the calculation process, launch it and update pressure index
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa GCMC @ {}K/{:.3f}bar (pressure {} of {})".format(
            self.ctx.temperature, self.ctx.pressures[self.ctx.current_p_index], self.ctx.current_p_index + 1,
            len(self.ctx.pressures)))
        self.ctx.current_p_index += 1
        return ToContext(raspa_gcmc=append_(running))

    def return_isotherm(self):
        """If is_porous and is_kh_enough create the isotherm_output Dict and report the pks"""

        gcmc_out_dict = {}
        for calc in self.ctx.raspa_gcmc:
            gcmc_out_dict[calc.label] = calc.outputs.output_parameters

        self.out(
            "isotherm_output",
            get_isotherm_output(self.ctx.parameters, self.ctx.raspa_widom.outputs.output_parameters, self.ctx.pressures,
                                **gcmc_out_dict))

        # Report the pk of the results (for geometric results show input if provided from IsothermMultiTemp)
        self.report("Isotherm @ {}K computed: geom Dict<{}>, isotherm Dict<{}>".format(
            self.ctx.temperature, self.ctx.geom.pk, self.outputs['isotherm_output'].pk))
