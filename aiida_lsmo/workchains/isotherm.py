# -*- coding: utf-8 -*-
"""Isothem workchain."""
from __future__ import absolute_import

import os

from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Float, Int, Str, List, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, while_, if_
from aiida_lsmo.utils import check_resize_unit_cell

# sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# calculation objects
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name

# data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  # pylint: disable=invalid-name


def choose_pressure_points(kh, qsat, dpa, dpmax, prange):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm in a list.

    :float kh: Henry coefficient (mol/kg/Pa)
    :float qsat: saturations loading (mol/kg)
    :float dpa: precision of the sampling at low pressure (0.1 is a good one)
    :float dpmax: maximum distance between two pressure points (bar)
    :list prange: [min,max] pressures to sample (bar)
    """
    pmin, pmax = prange
    b_value = kh / qsat * 1e5  #(1/bar)

    pressure_points = [pmin]
    while True:
        pold = pressure_points[-1]
        delta_p = min(dpmax, dpa * (b_value * pold**2 + 2 * pold + 1 / b_value))
        pnew = pold + delta_p
        if pnew <= pmax:
            pressure_points.append(pnew)
        else:
            pressure_points.append(pmax)
            return pressure_points


@calcfunction
def get_atomic_radii(filename):
    """Get {structure_radii}.rad as SinglefileData"""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    fullfilename = filename.value + ".rad"
    return SinglefileData(file=os.path.join(thisdir, "isotherm_data", fullfilename))


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
def get_zeopp_parameters(molecule_dict, volpo_samples, block_samples):
    """Get the ZeoppParameters from the inputs of the workchain"""
    probe_rad = molecule_dict["proberad"]
    param_dict = {
        'ha': 'DEF',
        'volpo': [probe_rad, probe_rad, volpo_samples.value],
        'block': [probe_rad, block_samples.value],
    }
    return ZeoppParameters(dict=param_dict)


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

        spec.input("forcefield", valid_type=Str, help='Forcefield used in Raspa')

        spec.input('structure_radii',
                   valid_type=(Str, SinglefileData),
                   help='Zeopp .rad definition of the structure, coherent with FF.' +
                   'If Str is used it will take {Str}.rad as Singlefile.')

        spec.input('temperature',
                   valid_type=Float,
                   default=Float(300),
                   required=False,
                   help='Temperature of the simulation')

        spec.input("zeopp_volpo_samples",
                   valid_type=Int,
                   default=Int(1e5),
                   required=False,
                   help='Number of samples for VOLPO calculation (per UC volume)')

        spec.input("zeopp_block_samples",
                   valid_type=Int,
                   default=Int(100),
                   required=False,
                   help='Number of samples for BLOCK calculation (per A^3)')

        spec.input("raspa_minKh",
                   valid_type=Float,
                   default=Float(1e-10),
                   required=False,
                   help='If Henry coefiicient < raspa_minKh do not run the isotherm (mol/kg/Pa)')

        spec.input("raspa_verbosity",
                   valid_type=Int,
                   default=Int(10),
                   required=False,
                   help='Print stats every: number of cycles / raspa_verbosity')

        spec.input("raspa_widom_cycles",
                   valid_type=Int,
                   default=Int(1e5),
                   required=False,
                   help='Number of widom cycles')

        spec.input("raspa_gcmc_init_cycles",
                   valid_type=Int,
                   default=Int(1e3),
                   required=False,
                   help='Number of GCMC initialization cycles')

        spec.input("raspa_gcmc_prod_cycles",
                   valid_type=Int,
                   default=Int(1e4),
                   required=False,
                   help='Number of GCMC production cycles')

        spec.input("raspa_cutoff",
                   valid_type=Float,
                   default=Float(12.0),
                   required=False,
                   help='CutOff truncation for the VdW interactions (Angstrom)')

        #TODO: allow to modify the probability list

        spec.input("pressure_list",
                   valid_type=List,
                   required=False,
                   help='Pressure list for the isotherm (bar): if given it will' +
                   'use this list instead of guessing the pressure points.')

        spec.input("pressure_precision",
                   valid_type=Float,
                   default=Float(0.1),
                   required=False,
                   help='Precision in the sampling of the isotherm:' +
                   '0.1 ok for full isotherm, 0.05 better for lowP range')

        spec.input("pressure_maxstep",
                   valid_type=Float,
                   default=Float(5),
                   required=False,
                   help='Max distance between pressure points (bar)')
        #
        spec.input("pressure_range",
                   valid_type=List,
                   default=List(list=[0.001, 30]),
                   required=False,
                   help='Range of min and max pressure (bar)')

        # workflow
        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            if_(cls.should_run_widom)(  # run Widom only if porous
                cls.run_raspa_widom,  # run raspa widom calculation
                if_(cls.should_run_gcmc)(  # Kh is high enough
                    cls.init_raspa_gcmc,  # initializate setting for GCMC
                    while_(cls.should_run_another_gcmc)(  # new pressure
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation
                        cls.parse_raspa_gcmc,  # parse the result @ T,P
                    ),
                ),
            ),
            cls.return_results,
        )

        spec.outputs.dynamic = True  # any outputs are accepted

    def setup(self):
        """Initialize the parameters"""

        # Get the molecule Dict
        # if self.inputs.molecule == Str # TODO
        self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
        # else: use Dict directly # TODO

    def run_zeopp(self):
        """Perform Zeo++ block and VOLPO calculations."""

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')

        # Set specific metadata
        inputs['metadata']['call_link_label'] = 'run_zeopp_block_and_volpo'
        inputs['metadata']['label'] = "ZeoppVolpoBlock"

        # Set input structure
        inputs['structure'] = self.inputs.structure

        # Set atomic_radii: use the SinglefileData or load it as {Str}.rad
        # if self.inputs.structure_radii == Str: # TODO:
        inputs['atomic_radii'] = get_atomic_radii(self.inputs.structure_radii)
        # else: simply load it TODO

        # Set parameters
        inputs['parameters'] = get_zeopp_parameters(self.ctx.molecule, self.inputs.zeopp_volpo_samples,
                                                    self.inputs.zeopp_block_samples)

        running = self.submit(ZeoppCalculation, **inputs)
        self.report("Running zeo++ block and volpo Calculation<{}>".format(running.id))
        return ToContext(zeopp=running)

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume,
        also check the number of blocking spheres and estimate the saturation loading"""

        self.ctx.is_porous = self.ctx.zeopp.outputs.output_parameters['POAV_Volume_fraction'] > 0.001

        if self.ctx.is_porous:
            self.report("Found accessible pore volume: continue")

            self.ctx.n_block_spheres = int(self.ctx.zeopp.outputs.block.get_content().splitlines()[0].strip())
            if self.ctx.n_block_spheres > 0:
                self.report("Found {} blocking spheres".format(n_block_spheres))
            else:
                self.report("No blocking spheres found")
        else:
            self.ctx.n_block_spheres = None
            self.report("No accessible pore volume: stop")

        # Estimate the total loading qsat (mol/kg) and choose the pressures
        # Note: cm3/g = l/kg and (mol/kg) = l/kg * mol/l
        self.ctx.estimated_qsat = \
            self.ctx.zeopp.outputs.output_parameters['POAV_cm^3/g'] * self.ctx.molecule['molsatdens']

        return self.ctx.is_porous

    def _get_widom_param(self):
        """Write Raspa input parameters from scratch, for a Widom calculation"""

        vf = self.ctx.zeopp.outputs.output_parameters["POAV_Volume_fraction"]
        param = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfInitializationCycles": 0,
                "NumberOfCycles": self.inputs.raspa_widom_cycles,  #toCheck: not sure I can avoid to use .value
                "PrintPropertiesEvery": self.inputs.raspa_widom_cycles / self.inputs.raspa_verbosity,
                "PrintEvery": int(1e6),
                "RemoveAtomNumberCodeFromLabel": True,  # be careful!
                "Forcefield": self.inputs.forcefield,
                "UseChargesFromCIFFile": "yes",
                "CutOff": self.inputs.raspa_cutoff,
            },
            "System": {
                "framework_1": {
                    "type": "Framework",
                    "HeliumVoidFraction": vf,
                    "ExternalTemperature": self.inputs.temperature,
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
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.inputs.raspa_cutoff)
        param["System"]["framework_1"]["UnitCells"] = "{} {} {}".format(mult[0], mult[1], mult[2])

        if self.ctx.n_block_spheres > 0:
            param["Component"][self.ctx.molecule['name']].update({"BlockPocketsFileName": "block_pocket"})

        if self.ctx.molecule['charged']:
            param["GeneralSettings"].update({"ChargeMethod": "Ewald", "EwaldPrecision": 1e-6})
        return param

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa."""

        # Initialize the input for raspa_base, that later will need only
        # minor updates
        self.ctx.inp = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.inp['metadata']['label'] = "RaspaWidom"
        self.ctx.inp['metadata']['call_link_label'] = "run_raspa_widom"

        self.ctx.inp['raspa']['framework'] = {"framework_1": self.inputs.structure}
        if self.ctx.n_block_spheres > 0:
            self.ctx.inp["raspa"]["block_pocket"] = self.ctx.zeopp.outputs.block

        self.ctx.raspa_param = self._get_widom_param()
        self.ctx.inp['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param).store()

        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa Widom for the Henry coefficient")

        return ToContext(raspa_widom=running)

    def should_run_gcmc(self):
        """Compute the isotherm only if the material meets the user defined criteria."""

        self.ctx.kh = self.ctx.raspa_widom.outputs.output_parameters["framework_1"]["components"][
            self.ctx.molecule['name']]['henry_coefficient_average']
        self.ctx.is_kh_ehough = self.ctx.kh > self.inputs.raspa_minKh
        if self.ctx.is_kh_ehough:
            self.report("kH larger than the threshold: continue")
        else:
            self.report("kHh lower than the threshold: stop")
        return self.ctx.is_kh_ehough

    def _update_param_for_gcmc(self):
        """Update Raspa input parameter, from Widom to GCMC"""

        param = self.ctx.raspa_param
        param["GeneralSettings"].update({
            "NumberOfInitializationCycles": self.inputs.raspa_gcmc_init_cycles,
            "NumberOfCycles": self.inputs.raspa_gcmc_prod_cycles,
            "PrintPropertiesEvery": int(1e6),
            "PrintEvery": self.inputs.raspa_gcmc_prod_cycles / self.inputs.raspa_verbosity,
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
        """Initialize variables and the pressures we want to compute."""

        # Get the self.ctx.pressures as a list
        if "pressure_list" in self.inputs:
            self.ctx.pressures = self.inputs.pressure_list.get_list()
        else:
            self.ctx.pressures = choose_pressure_points(
                self.ctx.kh,  # (mol/kg/Pa)
                self.ctx.estimated_qsat,  # (mol/kg_frame)
                self.inputs.pressure_precision.value,  # (kg*Pa/mol)
                self.inputs.pressure_maxstep.value,  # (bar)
                self.inputs.pressure_range.get_list()  # (bar)
            )

        # Initializate counter, and isothem_output dict
        self.ctx.current_p_index = 0
        self.ctx.isotherm_output = { # TODO: use for the the to_context approach
            'Pressure_(bar)' : self.ctx.pressures,
            'Loading_average_(mol/kg)': [],
            'Loading_deviation_(mol/kg)': [],
            'Enthalpy_of_adsorption_average_(kJ/mol)': [],
            'Enthalpy_of_adsorption_deviation_(kJ/mol)': [],
        }

        self.report("Computed Kh(mol/kg/Pa)={:.2e} POAV(cm3/g)={:.3f} Qsat(mol/kg)={:.2f}".format(
            self.ctx.kh,
            self.ctx.zeopp.outputs.output_parameters['POAV_cm^3/g'],
            self.ctx.estimated_qsat,
        ))
        self.report("Now evaluating the isotherm for {} pressure points".format(len(self.ctx.pressures)))

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
            self.ctx.inp['raspa']['retrieved_parent_folder'] = self.ctx.raspa_gcmc.outputs.retrieved

        # Create the calculation process and launch it
        running = self.submit(RaspaBaseWorkChain, **self.ctx.inp)
        self.report("Running Raspa GCMC at p(bar)={:.3f} ({} of {})".format(
            self.ctx.pressures[self.ctx.current_p_index], self.ctx.current_p_index + 1, len(self.ctx.pressures)))
        return ToContext(raspa_gcmc=running)

    def parse_raspa_gcmc(self):
        """Extract the pressure and loading average of the last completed Raspa calculation"""

        pressure = self.ctx.pressures[self.ctx.current_p_index]
        output_param = self.ctx.raspa_gcmc.outputs.output_parameters["framework_1"]

        conv_load = output_param["components"][self.ctx.molecule['name']]["conversion_factor_molec_uc_to_mol_kg"]
        self.ctx.isotherm_output['Loading_average_(mol/kg)'].append(
            conv_load * output_param["components"][self.ctx.molecule['name']]['loading_absolute_average'])
        self.ctx.isotherm_output['Loading_deviation_(mol/kg)'].append(
            conv_load * output_param["components"][self.ctx.molecule['name']]['loading_absolute_dev'])

        conv_ener = 1.0 / 120.273  # K to kJ/mol
        if output_param['general']['enthalpy_of_adsorption_average']:
            self.ctx.isotherm_output['Enthalpy_of_adsorption_average_(kJ/mol)'].append(
                    conv_ener * output_param['general']['enthalpy_of_adsorption_average'])
            self.ctx.isotherm_output['Enthalpy_of_adsorption_deviation_(kJ/mol)'].append(
                    conv_ener * output_param['general']['enthalpy_of_adsorption_dev'])
        else:
            self.ctx.isotherm_output['Enthalpy_of_adsorption_average_(kJ/mol)'].append(None)
            self.ctx.isotherm_output['Enthalpy_of_adsorption_deviation_(kJ/mol)'].append(None)
        # Update counter
        self.ctx.current_p_index += 1

    def return_results(self):
        """Collect the results"""

        # Zeopp section
        zeopp_out = self.ctx.zeopp.outputs.output_parameters
        self.out("geometric_output", Dict(dict={
            'Density': zeopp_out['Density'],
            'Density_unit': "g/cm^3",
            'POAV_Volume_fraction': zeopp_out['POAV_Volume_fraction'],
            'PONAV_Volume_fraction': zeopp_out['PONAV_Volume_fraction'],
            'POAV_cm^3/g': zeopp_out['POAV_cm^3/g'],
            'number_blocking_spheres': self.ctx.n_block_spheres,
            'estimated_saturation_loading_mol/kg': self.ctx.estimated_qsat,
            'is_porous': self.ctx.is_porous
        }).store()) #TODO: move to a calcfunction?

        # Raspa Widom section
        if self.ctx.is_porous:
            widom_out = self.ctx.raspa_widom.outputs.output_parameters["framework_1"]["components"][
                self.ctx.molecule["name"]]

            self.out("widom_output", Dict(dict={
                "Temperature_(K)": [int(round(self.inputs.temperature.value))],
                "{} K".format(int(round(self.inputs.temperature.value))): {
                    'Henry_coefficient_average_(mol/kg/Pa)': widom_out['henry_coefficient_average'],
                    'Henry_coefficient_deviation_(mol/kg/Pa)': widom_out['henry_coefficient_dev'],
                    'Adsorption_energy_average_(kJ/mol)': widom_out['adsorption_energy_widom_average'],
                    'Adsorption_energy_deviation_(kJ/mol)': widom_out['adsorption_energy_widom_dev'],
                    'is_kh_ehough': self.ctx.is_kh_ehough,
                }
            }).store()) #TODO: move to a calcfunction?

        # Raspa GCMC section
        if self.ctx.is_porous and self.ctx.is_kh_ehough:
            gcmc_out = self.ctx.raspa_gcmc.outputs.output_parameters["framework_1"]["components"][
                self.ctx.molecule['name']]

            self.out("isotherm_output", Dict(dict={
                "Temperature_(K)": [int(round(self.inputs.temperature.value))],
                "{} K".format(int(round(self.inputs.temperature.value))) : self.ctx.isotherm_output,
                'conversion_factor_molec_uc_to_cm3stp_cm3': gcmc_out['conversion_factor_molec_uc_to_cm3stp_cm3'],
                'conversion_factor_molec_uc_to_gr_gr': gcmc_out['conversion_factor_molec_uc_to_gr_gr'],
                'conversion_factor_molec_uc_to_mol_kg': gcmc_out['conversion_factor_molec_uc_to_mol_kg']
            }).store()) #TODO: move to a calcfunction?

        self.out("blocking_spheres", self.ctx.zeopp.outputs.block)
        self.report("Workchain completed: geom Dict<{}>, widom Dict<{}>, isotherm Dict<{}>".format(
            self.outputs['geometric_output'].pk,
            self.outputs['widom_output'].pk,
            self.outputs['isotherm_output'].pk))
