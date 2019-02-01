from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.orm.data.base import Float, Int
from aiida.work import workfunction as wf
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, Outputs, while_, if_
from copy import deepcopy

# subworkflows
from aiida_raspa.workflows import RaspaConvergeWorkChain

# calculation objects
ZeoppCalculation = CalculationFactory('zeopp.network')

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')


def choose_pressure_points(Kh, qsat, dpa, dpmax, pmax):
    """Model the isotherm as single-site langmuir and return the most important
    pressure points to evaluate for an isotherm. Returns a list of pressures.

    :param Kh: Henry coefficient (mol/kg/Pa)
    :param qsat: saturations loading (mol/kg)
    :param dpa: precision of the sampling at low pressure (0.1 is a good one)
    :param dpmax: maximum distance between two pressure points (Pa)
    :param pmax: max pressure to sample (Pa)
    """
    R = 8.314/1000 #(kJ/mol/K)
    b = Kh/qsat #(1/Pa)
    pmin = 0.001e5 #(Pa)
    p = [pmin]
    while True:
        pold = p[-1]
        dp = min(dpmax,dpa*(b*pold**2+2*pold+1/b))
        pnew = pold+dp
        if pnew <= pmax:
            p.append(pnew)
        else:
            p.append(pmax)
            return p

class VolpoKhIsothermWorkChain(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient.
    """
    @classmethod
    def define(cls, spec):
        super(VolpoKhIsothermWorkChain, cls).define(spec)

        # structure
        spec.input('structure', valid_type=CifData)

        # zeopp main
        spec.input('zeopp_code', valid_type=Code)
        spec.input("_zeopp_options", valid_type=dict, default=None, required=False)
        spec.input("zeopp_probe_radius", valid_type=Float)
        spec.input("zeopp_atomic_radii", valid_type=SinglefileData, default=None, required=False)

        # raspa main
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData)
        spec.input("_raspa_options", valid_type=dict, default=None, required=False)
        spec.input("_raspa_usecharges", valid_type=bool, default=False, required=False)
        spec.input("raspa_minKh", valid_type=Float, default=Float(1e-10), required=False)
        spec.input("raspa_molsatdens", valid_type=Float) #density of the liquid phase of the molecule in mol/m3

        # advanced settings (TODO: add and discuss)
        spec.input("zeopp_block_samples_A3", valid_type=Int, default=Int(100), required=False) #100 samples / Ang^3: accurate for all the structures
        spec.input("zeopp_volpo_samples_UC", valid_type=Int, default=Int(100000), required=False) #100k samples, may need more for structures bigger than 30x30x30
        spec.input("raspa_verbosity", valid_type=Int, default=Int(10), required=False)
        spec.input("raspa_widom_cycle_mult", valid_type=Int, default=Int(10), required=False)
        spec.input("raspa_gcmc_press_precision", valid_type=Float, default=Float(0.1), required=False)
        spec.input("raspa_gcmc_press_maxstep", valid_type=Float, default=Float(5e5), required=False)
        spec.input("raspa_gcmc_press_max", valid_type=Float, default=Float(30e5), required=False)

        # workflow
        spec.outline(
            cls.run_zeopp,                #computes volpo and blocks
            if_(cls.should_run_widom)(    #run Widom only if porous
                cls.init_raspa_widom,     #initialize raspa widom calculation
                cls.run_raspa_widom,      #run raspa widom calculation
                if_(cls.should_run_gcmc)(     #make decision (e.g., Kh high enough)
                    cls.init_raspa_gcmc,      #initializate setting for GCMC
                    while_(cls.should_run_another_gcmc)( #new pressure to compute
                        cls.run_raspa_gcmc,   #run raspa GCMC calculation
                        cls.parse_raspa_gcmc, #parse the result @ T,P
                    ),
                ),
            ),
            cls.return_results,
        )

        spec.dynamic_output()


    def run_zeopp(self):
        """Main function that performs zeo++ block and VOLPO calculations."""
        params = {
                'ha': True,
                'block': [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_block_samples_A3.value],
                'volpo': [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_volpo_samples_UC.value,]
        }

        inputs = {
            'code'         : self.inputs.zeopp_code,
            'structure'    : self.inputs.structure,
            'parameters'   : NetworkParameters(dict=params).store(),
            '_options'     : self.inputs._zeopp_options,
            '_label'       : "ZeoppVolpoBlock",
        }

        # Use default zeopp atomic radii only if a .rad file is not specified
        try:
            inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Create the calculation process and launch it
        running = submit(ZeoppCalculation.process(), **inputs)
        self.report("pk: {} | Running zeo++ block and volpo calculations".format(running.pid))
        return ToContext(zeopp=Outputs(running))


    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume."""
        if self.ctx.zeopp['output_parameters'].get_dict()['POAV_Volume_fraction'] > 1e-5:
            self.report("Found accessible pore volume: continue")
            return True
        else:
            self.report("NOT Found any accessible pore volume: stop")
            return False


    def init_raspa_widom(self):
        """Write the ParameterData for a Raspa and use the blocking spheres.
        REMEMBER that we want only one component for widom!
        """

        # Create a deepcopy of the user parameters, to modify before submission
        self.ctx.raspa_parameters = deepcopy(self.inputs.raspa_parameters.get_dict())

        # Turn on charges if requested
        if self.inputs._raspa_usecharges:
            self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "Ewald"
            self.ctx.raspa_parameters['GeneralSettings']['EwaldPrecision'] = 1e-6
            self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
        else:
            self.ctx.raspa_parameters['GeneralSettings']['ChargeMethod'] = "None"

        # CORRECT the settings to have only Widom insertion
        self.ctx.raspa_parameters["GeneralSettings"]["SimulationType"] = "MonteCarlo"
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = 0
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_widom_cycle_mult.value * self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] / self.inputs.raspa_verbosity.value)
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(1e6) #never
        self.ctx.raspa_parameters["Component"][0]["WidomProbability"] = 1.0
        return

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa."""

        # Create the inputs dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.inputs.structure,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            '_options'   : self.inputs._raspa_options,
            '_label'     : "RaspaWidom",
        }

        # Check if there are blocking spheres (reading the header of the file) and use them for Raspa
        with open(self.ctx.zeopp['block'].get_abs_path() + '/path/' + \
                  self.ctx.zeopp['block'].get_folder_list()[0]) as f:
            self.ctx.number_blocking_spheres = int(f.readline().strip())
        if self.ctx.number_blocking_spheres > 0:
            inputs['block_component_0'] = self.ctx.zeopp['block']
            self.report("Blocking spheres ({}) are present and used for Raspa".format(self.ctx.number_blocking_spheres))
        else:
            self.report("No blocking spheres found")

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running Raspa Widom for the Henry coefficient".format(running.pid))

        return ToContext(raspa_widom=Outputs(running))


    def should_run_gcmc(self):
        """Compute the isotherm only if the material meets the user defined
        criteria.
        """

        Kh = self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_average']
        if Kh > self.inputs.raspa_minKh.value:
            self.report("Kh larger than the threshold: compute isotherm")
            return True
        else:
            self.report("Kh lower than the threshold: don't compute isotherm")
            return False


    def init_raspa_gcmc(self):
        """Initialize variables and the pressures we want to compute."""

        # Initializate counter and set restart to None
        self.ctx.current_p_index = 0
        self.ctx.restart_raspa_calc = None

        # Estimate the total loading qsat and choose the pressure points
        satDens = self.inputs.raspa_molsatdens.value #(mol/l)
        poreVol = self.ctx.zeopp['output_parameters'].get_dict()['POAV_cm^3/g'] #(cm3/g = l/kg)
        self.ctx.estimated_qsat = satDens * poreVol
        self.ctx.pressures = choose_pressure_points(
            self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_average'], #(mol/kg/Pa)
            self.ctx.estimated_qsat, #(mol/kg_frame)
            self.inputs.raspa_gcmc_press_precision.value, #(kg*Pa/mol)
            self.inputs.raspa_gcmc_press_maxstep.value, #(Pa)
            self.inputs.raspa_gcmc_press_max.value, #(Pa)
            )
        self.report("Computed Kh(mol/kg/Pa)={:.2e} POAV(cm3/g)={:.3f} Qsat(mol/kg)={:.2f}".format(
                    self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_average'],
                    poreVol,
                    self.ctx.estimated_qsat,
                    ))
        self.report("Now evaluating the isotherm for {} pressure points".format(len(self.ctx.pressures)))

        # CORRECT the parameters to perform GCMC
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfInitializationCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfInitializationCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"] = self.inputs.raspa_parameters.get_dict()["GeneralSettings"]["NumberOfCycles"]
        self.ctx.raspa_parameters["GeneralSettings"]["PrintPropertiesEvery"] = int(1e6) #never
        self.ctx.raspa_parameters["GeneralSettings"]["PrintEvery"] = int(self.ctx.raspa_parameters["GeneralSettings"]["NumberOfCycles"]/self.inputs.raspa_verbosity.value)
        self.ctx.raspa_parameters["Component"][0]["WidomProbability"] = 0.0
        self.ctx.raspa_parameters["Component"][0]["TranslationProbability"] = 1.0
        self.ctx.raspa_parameters["Component"][0]["RotationProbability"] = 1.0
        self.ctx.raspa_parameters["Component"][0]["ReinsertionProbability"] = 1.0
        self.ctx.raspa_parameters["Component"][0]["SwapProbability"] = 2.0
        return


    def should_run_another_gcmc(self):
        """We run another raspa calculation only if the current iteration is
        smaller than the total number of pressures we want to compute.
        """
        return self.ctx.current_p_index < len(self.ctx.pressures)


    def run_raspa_gcmc(self):
        """Run a GCMC calculation in Raspa @ T,P. """
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure'] = pressure

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.inputs.structure,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters).store(),
            '_options'   : self.inputs._raspa_options,
            '_label'     : "RaspaGCMC",
        }
        # Check if there are poket blocks to be loaded
        if self.ctx.number_blocking_spheres > 0:
            inputs['block_component_0'] = self.ctx.zeopp['block']

        # Check if there is a previous calculation (lower p) to restart from
        if self.ctx.restart_raspa_calc is not None:
            inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running Raspa GCMC at p(bar)={:.3f} ({} of {})".format(running.pid, pressure/1e5, self.ctx.current_p_index+1, len(self.ctx.pressures)))
        return ToContext(raspa_gcmc=Outputs(running))


    def parse_raspa_gcmc(self):
        """Extract the pressure and loading average of the last completed raspa
        calculation.
        """

        # Initializate variables
        if self.ctx.current_p_index == 0:
            self.ctx.isotherm_loading = []
            self.ctx.isotherm_enthalpy = []

        # Store results
        pressure = self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure']/1e5
        conv1 = self.ctx.raspa_gcmc["component_0"].get_dict()["conversion_factor_molec_uc_to_mol_kg"]
        loading_average = conv1 * self.ctx.raspa_gcmc["component_0"].get_dict()['loading_absolute_average']
        loading_dev = conv1 * self.ctx.raspa_gcmc["component_0"].get_dict()['loading_absolute_dev']
        conv2 = 1/120.273 # K to kJ/mol
        enthalpy_of_adsorption = conv2 * self.ctx.raspa_gcmc["output_parameters"].get_dict()['enthalpy_of_adsorption_average']
        enthalpy_of_adsorption_dev = conv2 * self.ctx.raspa_gcmc["output_parameters"].get_dict()['enthalpy_of_adsorption_dev']
        self.ctx.isotherm_loading.append((pressure, loading_average,loading_dev))
        self.ctx.isotherm_enthalpy.append((pressure, enthalpy_of_adsorption, enthalpy_of_adsorption_dev))

        # Update counter and parent folder for restart
        self.ctx.current_p_index += 1
        self.ctx.restart_raspa_calc = self.ctx.raspa_gcmc['retrieved_parent_folder']
        return

    def return_results(self):
        """Attach the results to the output."""

        result_dict = {}

        # Zeopp section
        result_dict['Density'] = self.ctx.zeopp['output_parameters'].get_dict()['Density']
        result_dict['Density_unit'] = "g/cm^3"
        result_dict['POAV_Volume_fraction'] = self.ctx.zeopp['output_parameters'].get_dict()['POAV_Volume_fraction']
        result_dict['PONAV_Volume_fraction'] = self.ctx.zeopp['output_parameters'].get_dict()['PONAV_Volume_fraction']
        result_dict['POAV_cm^3/g'] = self.ctx.zeopp['output_parameters'].get_dict()['POAV_cm^3/g']
        try:
            result_dict['number_blocking_spheres'] = self.ctx.number_blocking_spheres
        except AttributeError:
            pass

        # Raspa Widom section
        try:
            result_dict['temperature'] = self.ctx.raspa_parameters["GeneralSettings"]["ExternalTemperature"]
            result_dict['temperature_unit'] = "K"
            result_dict['henry_coefficient_average'] = self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_average'] #(mol/kg/Pa)
            result_dict['henry_coefficient_dev'] = self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_dev']
            result_dict['henry_coefficient_units'] = self.ctx.raspa_widom["component_0"].get_dict()['henry_coefficient_units']
            result_dict['adsorption_energy_average'] = self.ctx.raspa_widom["component_0"].get_dict()['adsorption_energy_widom_average'] #(kJ/mol)
            result_dict['adsorption_energy_dev'] = self.ctx.raspa_widom["component_0"].get_dict()['adsorption_energy_widom_dev']
            result_dict['adsorption_energy_units'] = self.ctx.raspa_widom["component_0"].get_dict()['adsorption_energy_widom_units']
        except AttributeError:
            pass

        # Raspa GCMC section
        try:
            result_dict['estimated_saturation_loading'] = self.ctx.estimated_qsat
            result_dict['estimated_saturation_loading_units'] = "mol/kg"
            result_dict['isotherm_loading_header'] = ['Pressure(bar)', 'Loading_average(molec/UC)', 'Loading_deviation(molec/UC)']
            result_dict['isotherm_loading'] = self.ctx.isotherm_loading
            result_dict['isotherm_enthalpy_header'] = ['Pressure(bar)', 'Enthalpy_of_adsorption_average(kJ/mol)', 'Enthalpy_of_adsorption_deviation(kJ/mol)']
            result_dict['isotherm_enthalpy'] = self.ctx.isotherm_enthalpy
            result_dict['conversion_factor_molec_uc_to_cm3stp_cm3'] = self.ctx.raspa_gcmc["component_0"].get_dict()['conversion_factor_molec_uc_to_cm3stp_cm3']
            result_dict['conversion_factor_molec_uc_to_gr_gr'] = self.ctx.raspa_gcmc["component_0"].get_dict()['conversion_factor_molec_uc_to_gr_gr']
            result_dict['conversion_factor_molec_uc_to_mol_kg'] = self.ctx.raspa_gcmc["component_0"].get_dict()['conversion_factor_molec_uc_to_mol_kg']
        except AttributeError:
            pass

        self.out("results", ParameterData(dict=result_dict).store())
        self.out('blocking_spheres', self.ctx.zeopp['block'])
        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return
# EOF
