"""A work chain."""
import os
import ruamel.yaml as yaml

# AiiDA modules
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Dict, SinglefileData
from aiida.engine import calcfunction
from aiida.engine import WorkChain, if_
from aiida_lsmo.utils import aiida_dict_merge, check_resize_unit_cell

RaspaBaseWorkChain = WorkflowFactory('raspa.base')  #pylint: disable=invalid-name

# Defining DataFactory and CalculationFactory
CifData = DataFactory("cif")  #pylint: disable=invalid-name
ZeoppParameters = DataFactory("zeopp.parameters")  #pylint: disable=invalid-name

ZeoppCalculation = CalculationFactory("zeopp.network")  #pylint: disable=invalid-name
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name

# Default parameters
PARAMETERS_DEFAULT = Dict(
    dict={
        "ff_framework": "UFF",  # str, Forcefield of the structure (used also as a definition of ff.rad for zeopp)
        "ff_shifted": False,  # bool, Shift or truncate at cutoff
        "ff_tail_corrections": True,  # bool, Apply tail corrections
        "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
        "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "zeopp_block_scaling": 1.0,  # float, scaling probe's diameter: use 0.0 for skipping block calc
        "zeopp_block_samples": int(1000),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_gcmc_init_cycles": int(1e3),  # int, Number of GCMC initialization cycles
        "raspa_gcmc_prod_cycles": int(1e4),  # int, Number of GCMC production cycles
    })


@calcfunction
def get_components_dict(conditions, parameters):
    """Construct components dict, like:
    {
        'xenon': {
            'name': 'Xe',
            'molfraction': xxx,
            'proberad': xxx,
            'zeopp': {...},
            ...
        },
        ...
    }
    """

    components_dict = {}
    thisdir = os.path.dirname(os.path.abspath(__file__))
    yamlfile = os.path.join(thisdir, "isotherm_data", "isotherm_molecules.yaml")
    with open(yamlfile, 'r') as stream:
        yaml_dict = yaml.safe_load(stream)
    for key, value in conditions.get_dict()['molfraction'].items():
        components_dict[key] = yaml_dict[key]
        components_dict[key]['molfraction'] = value
        probe_rad = components_dict[key]['proberad']
        components_dict[key]['zeopp'] = {
            'ha': 'DEF',
            'block': [probe_rad * parameters['zeopp_block_scaling'], parameters['zeopp_block_samples']],
        }
    return Dict(dict=components_dict)


@calcfunction
def get_ff_parameters(components, isotparams):
    """Get the parameters for ff_builder."""
    ff_params = {
        'ff_framework': isotparams['ff_framework'],
        'ff_molecules': {},
        'shifted': isotparams['ff_shifted'],
        'tail_corrections': isotparams['ff_tail_corrections'],
        'mixing_rule': isotparams['ff_mixing_rule'],
        'separate_interactions': isotparams['ff_separate_interactions']
    }
    for value in components.get_dict().values():
        ff = value['forcefield']  #pylint: disable=invalid-name
        ff_params['ff_molecules'][value['name']] = ff
    return Dict(dict=ff_params)


@calcfunction
def get_atomic_radii(isotparam):
    """Get {ff_framework}.rad as SinglefileData form workchain/isotherm_data. If not existing use DEFAULT.rad."""
    thisdir = os.path.dirname(os.path.abspath(__file__))
    filename = isotparam['ff_framework'] + ".rad"
    filepath = os.path.join(thisdir, "isotherm_data", filename)
    if not os.path.isfile(filepath):
        filepath = os.path.join(thisdir, "isotherm_data", "DEFAULT.rad")
    return SinglefileData(file=filepath)


@calcfunction
def get_geometric_output(zeopp_out):
    """Return the geometric_output Dict from Zeopp results, including is_porous"""
    geometric_output = zeopp_out.get_dict()
    geometric_output.update({'is_porous': geometric_output["POAV_A^3"] > 0.000})
    return Dict(dict=geometric_output)


#pylint: disable = too-many-branches
@calcfunction
def get_output_parameters(inp_conditions, components, **all_out_dicts):
    """Extract results to output_parameters Dict."""
    out_dict = {}

    # Add geometric results from Zeo++
    for comp_dict in components.get_dict().values():
        comp = comp_dict['name']
        zeopp_label = "Zeopp_{}".format(comp)
        if zeopp_label in all_out_dicts:  # skip if Zeopp not computed
            for key in ['Input_block', 'Number_of_blocking_spheres']:
                if key not in out_dict:
                    out_dict[key] = {}
                out_dict[key][comp] = all_out_dicts[zeopp_label][key]

    # Add GCMC results
    key_system = list(all_out_dicts[f'RaspaGCMC_Ads'].get_dict().keys())[0]
    conv_ener = 1.0 / 120.273  # K to kJ/mol

    out_dict.update({
        'temperatures': [inp_conditions[x]['temperature'] for x in ['adsorption', 'desorption']],
        'temperatures_unit': 'K',
        'pressures': [inp_conditions[x]['pressure'] for x in ['adsorption', 'desorption']],
        'pressures_unit': 'bar',
        'composition': {},
        'loading_absolute_unit': 'mol/kg' if key_system == "framework_1" else "(cm^3_STP/cm^3)",
        'enthalpy_of_adsorption_unit': 'kJ/mol',
    })

    for i in ["Ads", "Des"]:
        gcmc_out = all_out_dicts[f'RaspaGCMC_{i}'][key_system]  # can be framework_1 or box_1
        for comp_dict in components.get_dict().values():
            comp = comp_dict['name']
            out_dict['composition'][comp] = comp_dict['molfraction']
            key_conv_load = "conversion_factor_molec_uc_to_mol_kg" \
                            if key_system == "framework_1" else "conversion_factor_molec_uc_to_cm3stp_cm3"
            conv_load = gcmc_out['components'][comp][key_conv_load]
            for label in [
                    'loading_absolute_average', 'loading_absolute_dev', 'enthalpy_of_adsorption_average',
                    'enthalpy_of_adsorption_dev'
            ]:
                if label not in out_dict:
                    out_dict[label] = {}
                if comp not in out_dict[label]:
                    out_dict[label][comp] = []
                if label.startswith('loading'):
                    out_dict[label][comp].append(conv_load * gcmc_out['components'][comp][label])
                elif label.startswith('enthalpy'):
                    if gcmc_out['components'][comp][label] is None:
                        gcmc_out['components'][comp][label] = 0
                    out_dict[label][comp] = conv_ener * gcmc_out['components'][comp][label]

    return Dict(dict=out_dict)


class MulticompAdsDesWorkChain(WorkChain):
    """Compute Adsorption/Desorption in crystalline materials,
    for a mixture of componentes and at specific temperature/pressure conditions.
    """

    @classmethod
    def define(cls, spec):
        super(MulticompAdsDesWorkChain, cls).define(spec)

        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])
        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])
        spec.input('structure', valid_type=CifData, help='Adsorbent framework CIF.')
        spec.input("conditions",
                   valid_type=Dict,
                   help='A dictionary of components with their corresponding mol fractions in the mixture.')
        spec.input("parameters",
                   valid_type=Dict,
                   help='It provides the parameters which control the decision making behavior of workchain.')
        spec.outline(
            cls.setup,
            if_(cls.should_run_zeopp)(cls.run_zeopp, cls.inspect_zeopp_calc),
            cls.run_raspa_gcmc_ads,
            cls.run_raspa_gcmc_des,  #desorpion
            cls.return_output_parameters,
        )

        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='Results of the single temperature multi component workchain')

        spec.output_namespace('block_files',
                              valid_type=SinglefileData,
                              required=False,
                              dynamic=True,
                              help='Generated block pocket files.')

    def setup(self):
        """Initialize parameters"""
        self.ctx.sim_in_box = "structure" not in self.inputs.keys()
        self.ctx.parameters = aiida_dict_merge(PARAMETERS_DEFAULT, self.inputs.parameters)
        self.ctx.components = get_components_dict(self.inputs.conditions, self.ctx.parameters)
        self.ctx.ff_params = get_ff_parameters(self.ctx.components, self.ctx.parameters)

    def should_run_zeopp(self):
        """Return if it should run zeopp calculation."""
        return not self.ctx.sim_in_box and self.ctx.parameters['zeopp_block_scaling'] > 0.0

    def run_zeopp(self):
        """It performs the full zeopp calculation for all components."""
        zeopp_inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')
        zeopp_inputs.update({'structure': self.inputs.structure, 'atomic_radii': get_atomic_radii(self.ctx.parameters)})
        for key, value in self.ctx.components.get_dict().items():
            comp = value['name']
            zeopp_label = "Zeopp_{}".format(comp)
            zeopp_inputs['metadata']['label'] = zeopp_label
            zeopp_inputs['metadata']['call_link_label'] = "run_{}".format(zeopp_label)
            zeopp_inputs['parameters'] = ZeoppParameters(dict=self.ctx.components[key]['zeopp'])
            running = self.submit(ZeoppCalculation, **zeopp_inputs)
            self.report(f"Running zeo++ block calculation<{running.id}> for {comp}")
            self.to_context(**{zeopp_label: running})

    def inspect_zeopp_calc(self):
        """Asserts whether all widom calculations are finished ok. If so, manage zeopp results."""
        for value in self.ctx.components.get_dict().values():
            assert self.ctx["Zeopp_{}".format(value['name'])].is_finished_ok

        for value in self.ctx.components.get_dict().values():
            comp = value['name']
            zeopp_label = f"Zeopp_{comp}"
            if self.ctx[zeopp_label].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
                self.out(f"block_files.{comp}_block_file", self.ctx[zeopp_label].outputs.block)

    def _get_gcmc_inputs_adsorption(self):
        """Generate Raspa input parameters from scratch, for a multicomponent GCMC calculation."""
        self.ctx.raspa_inputs = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.raspa_inputs['raspa']['file'] = FFBuilder(self.ctx.ff_params)
        self.ctx.raspa_inputs["raspa"]["block_pocket"] = {}
        self.ctx.raspa_param = {
            "GeneralSettings": {
                "SimulationType": "MonteCarlo",
                "NumberOfInitializationCycles": self.ctx.parameters['raspa_gcmc_init_cycles'],
                "NumberOfCycles": self.ctx.parameters['raspa_gcmc_prod_cycles'],
                "PrintPropertiesEvery": int(1e10),
                "PrintEvery": self.ctx.parameters['raspa_gcmc_prod_cycles'] / self.ctx.parameters['raspa_verbosity'],
                "RemoveAtomNumberCodeFromLabel": True,  # github.com/aiidateam/aiida-core/issues/3304
                "Forcefield": "Local",
                "UseChargesFromCIFFile": "yes",
                "CutOff": self.ctx.parameters['ff_cutoff'],
                "ChargeMethod": "None",  #will be updated later if any molecule charged
            },
            "System": {},
            "Component": {}
        }

        if self.ctx.sim_in_box:
            boxlength = 3 * self.ctx.parameters['ff_cutoff']  #3x to avoid Coulomb self-interactions (sure?)
            self.ctx.raspa_param['System'] = {
                "box_1": {
                    "type": "Box",
                    "BoxLengths": f"{boxlength} {boxlength} {boxlength}"
                }
            }
        else:
            mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
            self.ctx.raspa_param['System'] = {
                "framework_1": {
                    "type": "Framework",
                    "UnitCells": f"{mult[0]} {mult[1]} {mult[2]}"
                }
            }
            self.ctx.raspa_inputs['raspa']['framework'] = {"framework_1": self.inputs.structure}

        system_key = list(self.ctx.raspa_param['System'].keys())[0]
        temp = self.inputs.conditions['adsorption']['temperature']  #K
        press = self.inputs.conditions['adsorption']['pressure'] * 1e5  #Ps
        self.ctx.raspa_param['System'][system_key]['ExternalTemperature'] = temp
        self.ctx.raspa_param['System'][system_key]['ExternalPressure'] = press

        for comp_dict in self.ctx.components.get_dict().values():
            comp = comp_dict['name']
            zeopp_label = "Zeopp_{}".format(comp)
            self.ctx.raspa_param["Component"][comp] = {
                "MoleculeDefinition": "Local",
                "MolFraction": comp_dict['molfraction'],
                "TranslationProbability": 1.0,
                "ReinsertionProbability": 1.0,
                "SwapProbability": 2.0,
                "IdentityChangeProbability": 2.0,
                "NumberOfIdentityChanges": len(list(self.ctx.components.get_dict())),  # todofuture: remove self-change
                "IdentityChangesList": [i for i in range(len(list(self.ctx.components.get_dict())))]
            }
            if comp_dict['charged']:  # will switch on if any molecule charged
                self.ctx.raspa_param["GeneralSettings"]['ChargeMethod'] = 'Ewald'
            if not comp_dict['singlebead']:
                self.ctx.raspa_param["Component"][comp]["RotationProbability"] = 1.0
            if zeopp_label in self.ctx and \
               self.ctx[zeopp_label].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
                self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"] = {"framework_1": comp + "_block_file"}
                self.ctx.raspa_inputs["raspa"]["block_pocket"][comp +
                                                               "_block_file"] = self.ctx[zeopp_label].outputs.block

    def _update_gcmc_inputs_desorption(self):
        """Update Raspa input parameters for desorption: Temperature, Pressure, Composition and Restart."""
        # Update temperature and pressure
        system_key = list(self.ctx.raspa_param['System'].keys())[0]
        temp = self.inputs.conditions['desorption']['temperature']  #K
        press = self.inputs.conditions['desorption']['pressure'] * 1e5  #Ps
        self.ctx.raspa_param['System'][system_key]['ExternalTemperature'] = temp
        self.ctx.raspa_param['System'][system_key]['ExternalPressure'] = press

        # Update composition reading the output from adsorption
        ads_results = self.ctx['RaspaGCMC_Ads'].outputs.output_parameters[system_key]['components']
        ads_molecxuc_tot = sum([x["loading_absolute_average"] for x in ads_results.values()])
        for key, val in ads_results.items():
            self.ctx.raspa_param["Component"][key]["MolFraction"] = val["loading_absolute_average"] / ads_molecxuc_tot

        # Link to restart file
        self.ctx.raspa_inputs["raspa"]['retrieved_parent_folder'] = self.ctx['RaspaGCMC_Ads'].outputs.retrieved

    def run_raspa_gcmc_ads(self):
        """Submit Raspa GCMC with adsorption T, P and composition."""
        self._get_gcmc_inputs_adsorption()
        gcmc_label = "RaspaGCMC_Ads"
        self.ctx.raspa_inputs['metadata']['label'] = gcmc_label
        self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_raspa_gcmc_ads"
        self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
        running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
        self.report("Running Raspa GCMC @ Adsorption conditions")
        self.to_context(**{gcmc_label: running})

    def run_raspa_gcmc_des(self):
        """Submit Raspa GCMC with adsorption T, P and composition."""
        self._update_gcmc_inputs_desorption()
        gcmc_label = "RaspaGCMC_Des"
        self.ctx.raspa_inputs['metadata']['label'] = gcmc_label
        self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_raspa_gcmc_des"
        self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
        running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
        self.report("Running Raspa GCMC @ Adsorption conditions")
        self.to_context(**{gcmc_label: running})

    def return_output_parameters(self):
        """Merge all the parameters into output_parameters, depending on is_porous and is_kh_ehough."""

        all_out_dicts = {}

        for key, val in self.ctx.items():
            if key.startswith('Zeopp_'):
                all_out_dicts[key] = val.outputs.output_parameters
            elif key.startswith('RaspaGCMC_'):
                all_out_dicts[key] = val.outputs.output_parameters
        self.out(
            "output_parameters",
            get_output_parameters(inp_conditions=self.inputs.conditions,
                                  components=self.ctx.components,
                                  **all_out_dicts))

        self.report("Workchain completed: output parameters Dict<{}>".format(self.outputs['output_parameters'].pk))
