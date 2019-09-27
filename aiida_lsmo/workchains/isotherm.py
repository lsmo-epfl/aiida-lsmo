# -*- coding: utf-8 -*-
"""Isothem workchain."""
from __future__ import absolute_import
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, Float, Int
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, while_, if_
from aiida_lsmo.utils import aiida_dict_merge, check_resize_unit_cell

# sub-workchains
RaspaBaseWorkChain = WorkflowFactory('raspa.base')  # pylint: disable=invalid-name

# calculation objects
ZeoppCalculation = CalculationFactory('zeopp.network')  # pylint: disable=invalid-name

# data objects
CifData = DataFactory('cif')  # pylint: disable=invalid-name


def choose_pressure_points(k_henry, qsat, dpa, dpmax, pmax):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm. Returns a list of pressures.

    :param k_henry: Henry coefficient (mol/kg/Pa)
    :param qsat: saturations loading (mol/kg)
    :param dpa: precision of the sampling at low pressure (0.1 is a good one)
    :param dpmax: maximum distance between two pressure points (Pa)
    :param pmax: max pressure to sample (Pa)
    """
    b_value = k_henry / qsat  #(1/Pa)
    pmin = 0.001e5  #(Pa)

    pressure_points = [pmin]
    while True:
        pold = pressure_points[-1]
        delta_p = min(dpmax,
                      dpa * (b_value * pold**2 + 2 * pold + 1 / b_value))
        pnew = pold + delta_p
        if pnew <= pmax:
            pressure_points.append(pnew)
        else:
            pressure_points.append(pmax)
            return pressure_points


@calcfunction
def multiply_unit_cell(params, structure):
    """Return parameter Dict object that is multiplied according to the cutoff value."""
    cutoff = params["GeneralSettings"]["CutOff"]
    new_params = params.get_dict()
    mult = check_resize_unit_cell(structure, 2 * cutoff)
    framework_name = list(params['System'].keys())[0]
    new_params["System"][framework_name]["UnitCells"] = "{} {} {}".format(
        mult[0], mult[1], mult[2])
    return Dict(dict=new_params)


@calcfunction
def update_aiida_dict(input_dict, updated_parameters):
    return aiida_dict_merge(input_dict, updated_parameters)


class IsothermWorkChain(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient.
    """
    ACC_VOL_FRAC_CRITERIA = 1e-5

    @classmethod
    def define(cls, spec):
        super(IsothermWorkChain, cls).define(spec)
        # structure
        spec.input('structure', valid_type=CifData)

        # zeopp main
        spec.expose_inputs(ZeoppCalculation,
                           namespace='zeopp',
                           exclude=['structure'])

        # raspa main
        spec.expose_inputs(RaspaBaseWorkChain,
                           namespace='raspa_base',
                           exclude=['raspa.framework'])

        #
        spec.input(
            "raspa_minKh",
            valid_type=Float,
            default=Float(1e-10),
            required=False
        )  # Henry coefiicient limit to identify whether to run isotherm or not

        # advanced settings (TODO: add and discuss)
        spec.input(
            "raspa_verbosity", valid_type=Int, default=Int(10), required=False
        )  # How rare do we want to print the output properties: number of cycles / raspa_verbosity
        spec.input("raspa_widom_cycle_mult",
                   valid_type=Int,
                   default=Int(10),
                   required=False
                   )  # Multiplicator for the number of widom insertion cylces

        # Parameters required to set the pressure points where to run the simulations
        spec.input("raspa_molsatdens", valid_type=Float
                   )  # Density of the liquid phase of the molecule in mol/m3
        spec.input(
            "raspa_gcmc_press_precision",
            valid_type=Float,
            default=Float(0.1),
            required=False
        )  # Precision in the sampling of the isotherm: 0.1 ok for full isotherm, 0.01 better for lowP range
        spec.input("raspa_gcmc_press_maxstep",
                   valid_type=Float,
                   default=Float(5e5),
                   required=False)  # Max distance between pressure points (Pa)
        spec.input("raspa_gcmc_press_max",
                   valid_type=Float,
                   default=Float(30e5),
                   required=False)  # Max P of the isotherm (Pa)

        # workflow
        spec.outline(
            cls.setup,
            cls.run_zeopp,  # computes volpo and blocks
            if_(cls.should_run_widom)(  # run Widom only if porous
                cls.run_raspa_widom,  # run raspa widom calculation
                if_(cls.should_run_gcmc)
                (  # make decision (e.g., Kh is high enough)
                    cls.init_raspa_gcmc,  # initializate setting for GCMC
                    while_(cls.should_run_another_gcmc)
                    (  # new pressure to compute
                        cls.run_raspa_gcmc,  # run raspa GCMC calculation
                        cls.parse_raspa_gcmc,  # parse the result @ T,P
                    ),
                ),
            ),
            cls.return_results,
        )

        spec.outputs.dynamic = True  # any outputs are accepted

    def setup(self):
        """Check that everything is ready to run the simulations."""
        # Check frameworks
        if len(self.inputs.raspa_base.raspa.parameters['Component'].keys()
               ) != 1:
            raise ValueError("We accept one framework only")
        self.ctx.framework_name = list(
            self.inputs.raspa_base.raspa.parameters['System'].keys())[0]

        # Check components
        if len(self.inputs.raspa_base.raspa.parameters['Component'].keys()
               ) != 1:
            raise ValueError("We accept one component only")
        self.ctx.component_name = list(
            self.inputs.raspa_base.raspa.parameters['Component'].keys())[0]

        # Check zeopp inputs
        if 'block' not in self.inputs.zeopp.parameters.get_dict():
            raise ValueError(
                "Block pocket calculation must be requested. Please provide 'block' parameter in the Zeo++ input."
            )

        if 'volpo' not in self.inputs.zeopp.parameters.get_dict():
            raise ValueError(
                "Accessible volume calculation calculation must be requested. "
                "Please provide 'volpo' parameter in the Zeo++ input.")

        if 'ha' not in self.inputs.zeopp.parameters.get_dict():
            raise ValueError(
                "Please enable high accuracty. Put 'ha': 'DEF' in the Zeo++ input"
            )

        self.ctx.print_every = int(
            self.inputs.raspa_base.raspa.parameters["GeneralSettings"]
            ["NumberOfCycles"] / self.inputs.raspa_verbosity.value)

        self.ctx.raspa_input_parameters = multiply_unit_cell(
            self.inputs.raspa_base.raspa.parameters, self.inputs.structure)

        # Initializate counter and set restart to None. Used only for the isotherm calculations
        self.ctx.current_p_index = 0
        self.ctx.restart_raspa_calc = None

    def run_zeopp(self):
        """Main function that performs zeo++ block and VOLPO calculations."""
        try:
            inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')
        except AttributeError:
            raise AttributeError(
                'No calculation input dictionary was defined in `self.inputs.zeopp`'
            )

        # Set input structure
        inputs['structure'] = self.inputs.structure

        # Set the `CALL` link label and calculation label
        inputs['metadata']['call_link_label'] = 'run_zeopp_block_and_volpo'

        # If the label is not provided - use the default value
        if "label" not in inputs['metadata']:
            inputs['metadata']['label'] = "ZeoppVolpoBlock"

        running = self.submit(ZeoppCalculation, **inputs)
        self.report(
            "pk: {} | Running zeo++ block and volpo calculations".format(
                running.id))
        return ToContext(zeopp=running)

    def _prepare_raspa(self,
                       label,
                       call_link_label,
                       updated_parameters,
                       block,
                       silent=True):
        """Prepare inputs for a raspa calculation."""
        try:
            inputs = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        except AttributeError:
            raise AttributeError(
                'No calculation input dictionary was defined in `self.inputs.raspa_base`'
            )

        # Set input framework
        inputs['raspa']['framework'] = {
            self.ctx.framework_name: self.inputs.structure
        }

        if "label" not in inputs['metadata']:
            inputs['metadata']['label'] = label

        # Set the `CALL` link label and calculation label
        inputs['metadata']['call_link_label'] = call_link_label

        # Block pockets
        self.ctx.number_blocking_spheres = int(
            block.get_content().splitlines()[0].strip())
        if self.ctx.number_blocking_spheres > 0:
            updated_parameters[self.ctx.component_name] = {
                "BlockPocketsFileName": "block_pocket"
            }
            inputs["raspa"]["block_pocket"] = self.ctx.zeopp.outputs.block
            if not silent:
                self.report(
                    "{} blocking spheres are present and used for Raspa".
                    format(self.ctx.number_blocking_spheres))
        else:
            if not silent:
                self.report("No blocking spheres found")

        inputs['raspa']['parameters'] = update_aiida_dict(
            self.ctx.raspa_input_parameters,
            Dict(dict=updated_parameters).store())

        return inputs

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume."""
        if self.ctx.zeopp.outputs.output_parameters[
                'POAV_Volume_fraction'] > self.ACC_VOL_FRAC_CRITERIA:
            self.report(
                "Found accessible pore volume: continue with widom calculation"
            )
            return True
        self.report(
            "Accessible pore volume wasn't found: stopping the work chain")
        return False

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa."""

        # Updated parametes to be used in raspa widom calculations
        updated_parameters = {
            "GeneralSettings": {
                "SimulationType":
                "MonteCarlo",
                "NumberOfInitializationCycles":
                0,
                "NumberOfCycles":
                self.inputs.raspa_widom_cycle_mult.value *
                self.inputs.raspa_base.raspa.parameters["GeneralSettings"]
                ["NumberOfCycles"],
                "PrintPropertiesEvery":
                self.ctx.print_every,
                "PrintEvery":
                int(1e6),
            },
            "Component": {
                self.ctx.component_name: {
                    "WidomProbability": 1.0,
                    "TranslationProbability": 0.0,
                    "ReinsertionProbability": 0.0,
                    "SwapProbability": 0.0,
                },
            },
        }
        inputs = self._prepare_raspa(label="RaspaWidom",
                                     call_link_label="run_raspa_widom",
                                     updated_parameters=updated_parameters,
                                     block=self.ctx.zeopp.outputs.block,
                                     silent=False)
        running = self.submit(RaspaBaseWorkChain, **inputs)
        self.report(
            "pk: {} | Running Raspa Widom for the Henry coefficient".format(
                running.id))

        return ToContext(raspa_widom=running)

    def should_run_gcmc(self):
        """Compute the isotherm only if the material meets the user defined criteria."""

        k_h = self.ctx.raspa_widom.outputs.output_parameters[
            self.ctx.framework_name]["components"][
                self.ctx.component_name]['henry_coefficient_average']
        if k_h > self.inputs.raspa_minKh.value:
            self.report("Kh larger than the threshold: compute isotherm")
            return True
        self.report("Kh lower than the threshold: don't compute isotherm")
        return False

    def init_raspa_gcmc(self):
        """Initialize variables and the pressures we want to compute."""

        self.ctx.isotherm_loading = []
        self.ctx.isotherm_enthalpy = []

        # Estimate the total loading qsat and choose the pressure points
        pore_vol = self.ctx.zeopp.outputs.output_parameters[
            'POAV_cm^3/g']  # (cm3/g = l/kg)
        self.ctx.estimated_qsat = pore_vol * self.inputs.raspa_molsatdens.value  # (mol/kg) = l/kg * mol/l

        self.ctx.pressures = choose_pressure_points(
            self.ctx.raspa_widom.outputs.output_parameters[
                self.ctx.framework_name]["components"][self.ctx.component_name]
            ['henry_coefficient_average'],  # (mol/kg/Pa)
            self.ctx.estimated_qsat,  # (mol/kg_frame)
            self.inputs.raspa_gcmc_press_precision.value,  # (kg*Pa/mol)
            self.inputs.raspa_gcmc_press_maxstep.value,  # (Pa)
            self.inputs.raspa_gcmc_press_max.value,  # (Pa)
        )

        self.report(
            "Computed Kh(mol/kg/Pa)={:.2e} POAV(cm3/g)={:.3f} Qsat(mol/kg)={:.2f}"
            .format(
                self.ctx.raspa_widom.outputs.output_parameters[
                    self.ctx.framework_name]["components"][
                        self.ctx.component_name]['henry_coefficient_average'],
                pore_vol,
                self.ctx.estimated_qsat,
            ))
        self.report(
            "Now evaluating the isotherm for {} pressure points".format(
                len(self.ctx.pressures)))

    def should_run_another_gcmc(self):
        """We run another raspa calculation only if the current iteration is
        smaller than the total number of pressures we want to compute.
        """
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_raspa_gcmc(self):
        """Run a GCMC calculation in Raspa @ T,P. """

        # Parameters to perform GCMC
        updated_parameters = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "PrintPropertiesEvery": int(1e6),
                "PrintEvery": self.ctx.print_every,
                "WidomProbability": 1.0,
            },
            "System": {
                self.ctx.framework_name: {
                    "ExternalPressure":
                    self.ctx.pressures[self.ctx.current_p_index],
                }
            },
            "Component": {
                self.ctx.component_name: {
                    "WidomProbability": 0.0,
                    "TranslationProbability": 1.0,
                    "ReinsertionProbability": 1.0,
                    "SwapProbability": 2.0,
                },
            },
        }

        inputs = self._prepare_raspa(label="RaspaGCMC",
                                     call_link_label="run_raspa_gcmc",
                                     updated_parameters=updated_parameters,
                                     block=self.ctx.zeopp.outputs.block)

        # Check if there is a previous calculation (lower p) to restart from
        if self.ctx.restart_raspa_calc is not None:
            inputs['raspa'][
                'retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        running = self.submit(RaspaBaseWorkChain, **inputs)
        self.report(
            "pk: {} | Running Raspa GCMC at p(bar)={:.3f} ({} of {})".format(
                running.id, self.ctx.pressures[self.ctx.current_p_index] / 1e5,
                self.ctx.current_p_index + 1, len(self.ctx.pressures)))
        return ToContext(raspa_gcmc=running)

    def parse_raspa_gcmc(self):
        """Extract the pressure and loading average of the last completed raspa calculation."""
        pressure = self.ctx.pressures[self.ctx.current_p_index] / 1e5
        output_params = self.ctx.raspa_gcmc.outputs.output_parameters[
            self.ctx.framework_name]
        conv1 = output_params["components"][
            self.ctx.component_name]["conversion_factor_molec_uc_to_mol_kg"]
        try:
            loading_average = conv1 * output_params["components"][
                self.ctx.component_name]['loading_absolute_average']
            loading_dev = conv1 * [
                "components"
            ][self.ctx.component_name]['loading_absolute_dev']
        except TypeError:
            loading_average = None
            loading_dev = None

        conv2 = 1.0 / 120.273  # K to kJ/mol
        try:
            enthalpy_of_adsorption = conv2 * output_params['general'][
                'enthalpy_of_adsorption_average']
            enthalpy_of_adsorption_dev = conv2 * output_params['general'][
                'enthalpy_of_adsorption_dev']
        except TypeError:
            enthalpy_of_adsorption = None
            enthalpy_of_adsorption_dev = None

        self.ctx.isotherm_loading.append(
            (pressure, loading_average, loading_dev))
        self.ctx.isotherm_enthalpy.append(
            (pressure, enthalpy_of_adsorption, enthalpy_of_adsorption_dev))

        # Update counter and parent folder for restart
        self.ctx.current_p_index += 1
        self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc.outputs.retrieved

    def return_results(self):
        """Attach the results to the output."""

        zeopp_output = self.ctx.zeopp.outputs.output_parameters
        raspa_widom_output = self.ctx.raspa_widom.outputs.output_parameters[
            self.ctx.framework_name]["components"][self.ctx.component_name]
        raspa_gcmc_output = self.ctx.raspa_gcmc.outputs.output_parameters[
            self.ctx.framework_name]["components"][self.ctx.component_name]
        result_dict = {}

        # Zeopp section
        result_dict.update({
            'Density':
            zeopp_output['Density'],
            'Density_unit':
            "g/cm^3",
            'POAV_Volume_fraction':
            zeopp_output['POAV_Volume_fraction'],
            'PONAV_Volume_fraction':
            zeopp_output['PONAV_Volume_fraction'],
            'POAV_cm^3/g':
            zeopp_output['POAV_cm^3/g'],
        })

        # Block pockets
        try:
            result_dict[
                'number_blocking_spheres'] = self.ctx.number_blocking_spheres
        except AttributeError:
            pass

        # Raspa Widom section
        try:
            result_dict.update({
                'temperature':
                self.inputs.raspa_base.raspa.parameters["System"][
                    self.ctx.framework_name]["ExternalTemperature"],
                'temperature_unit':
                "K",
                'henry_coefficient_average':
                raspa_widom_output['henry_coefficient_average'],  #(mol/kg/Pa)
                'henry_coefficient_dev':
                raspa_widom_output['henry_coefficient_dev'],
                'henry_coefficient_units':
                raspa_widom_output['henry_coefficient_units'],
                'adsorption_energy_average':
                raspa_widom_output[
                    'adsorption_energy_widom_average'],  #(kJ/mol)
                'adsorption_energy_dev':
                raspa_widom_output['adsorption_energy_widom_dev'],
                'adsorption_energy_units':
                raspa_widom_output['adsorption_energy_widom_units'],
            })
        except AttributeError:
            pass

        # Raspa GCMC section
        try:
            result_dict.update({
                'estimated_saturation_loading':
                self.ctx.estimated_qsat,
                'estimated_saturation_loading_units':
                "mol/kg",
                'isotherm_loading_header': [
                    'Pressure(bar)', 'Loading_average(molec/UC)',
                    'Loading_deviation(molec/UC)'
                ],
                'isotherm_loading':
                self.ctx.isotherm_loading,
                'isotherm_enthalpy_header': [
                    'Pressure(bar)', 'Enthalpy_of_adsorption_average(kJ/mol)',
                    'Enthalpy_of_adsorption_deviation(kJ/mol)'
                ],
                'isotherm_enthalpy':
                self.ctx.isotherm_enthalpy,
                'conversion_factor_molec_uc_to_cm3stp_cm3':
                raspa_gcmc_output['conversion_factor_molec_uc_to_cm3stp_cm3'],
                'conversion_factor_molec_uc_to_gr_gr':
                raspa_gcmc_output['conversion_factor_molec_uc_to_gr_gr'],
                'conversion_factor_molec_uc_to_mol_kg':
                raspa_gcmc_output['conversion_factor_molec_uc_to_mol_kg']
            })
        except AttributeError:
            pass

        self.out("results", Dict(dict=result_dict).store())
        self.out('blocking_spheres', self.ctx.zeopp.outputs.block)
        self.report("Workchain completed successfully")


# EOF
