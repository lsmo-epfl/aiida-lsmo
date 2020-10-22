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
ISOTHERMPARAMETERS_DEFAULT = Dict(
    dict={  #TODO: create IsothermParameters instead of Dict # pylint: disable=fixme
        "ff_framework": "UFF",  # str, Forcefield of the structure (used also as a definition of ff.rad for zeopp)
        "ff_shifted": False,  # bool, Shift or truncate at cutoff
        "ff_tail_corrections": True,  # bool, Apply tail corrections
        "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
        "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "zeopp_volpo_samples": int(1e5),  # int, Number of samples for VOLPO calculation (per UC volume)
        "zeopp_sa_samples": int(1e5),  # int, Number of samples for VOLPO calculation (per UC volume)
        "zeopp_block_samples": int(100),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_widom_cycles": int(1e5),  # int, Number of widom cycles
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
            'block': [probe_rad * 0.7, parameters['zeopp_block_samples']],  #Probe scaled by 0.7
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


#pylint: disable = too-many-branches
@calcfunction
def get_output_parameters(inp_conditions, components, **all_out_dicts):
    """ Extract Widom and GCMC results to isotherm Dict """
    print(all_out_dicts.keys())
    out_dict = {'components': [], 'geometric': {}, 'widom': {}, 'gcmc': {}}

    # Add geometric results from Zeo++
    for comp_dict in components.get_dict().values():
        comp = comp_dict['name']
        out_dict['components'].append(comp)
        zeopp_label = "Zeopp_{}".format(comp)
        for key, val in all_out_dicts[zeopp_label].get_dict().items():
            if key not in out_dict['geometric']:
                out_dict['geometric'][key] = {}
            out_dict['geometric'][key][comp] = val

    # Add Widom results
    out_dict['widom'] = {
        'temperatures': inp_conditions['t_widom'],
        'temperatures_unit': 'K',
        'henry_coefficient_unit': 'mol/kg/Pa',
        'adsorption_energy_widom_unit': 'kJ/mol',
    }

    widom_keys = [
        'henry_coefficient_average',
        'henry_coefficient_dev',
        'adsorption_energy_widom_average',
        'adsorption_energy_widom_dev',
    ]

    for comp_dict in components.get_dict().values():
        comp = comp_dict['name']
        for i, _ in enumerate(inp_conditions['t_widom']):
            widom_label = "RaspaWidom_{}_{}".format(comp, i)
            output_widom = all_out_dicts[widom_label].get_dict()
            for key in widom_keys:
                if key not in out_dict['widom']:
                    out_dict['widom'][key] = {}
                if comp not in out_dict['widom'][key]:
                    out_dict['widom'][key][comp] = []
                out_dict['widom'][key][comp].append(output_widom['framework_1']['components'][comp][key])

    # Add GCMC results
    conv_ener = 1.0 / 120.273  # K to kJ/mol

    out_dict['gcmc'] = {
        'temperatures': [tp[0] for tp in inp_conditions['tp_gcmc']],
        'temperatures_unit': 'K',
        'pressures': [tp[1] for tp in inp_conditions['tp_gcmc']],
        'pressures_unit': 'bar',
        'composition': {},
        'loading_absolute_unit': 'mol/kg',
        'enthalpy_of_adsorption_unit': 'kJ/mol',
    }

    for i, _ in enumerate(inp_conditions['tp_gcmc']):
        gcmc_out = all_out_dicts['RaspaGCMC_{}'.format(i)]["framework_1"]
        for comp_dict in components.get_dict().values():
            comp = comp_dict['name']
            out_dict['gcmc']['composition'][comp] = comp_dict['molfraction']
            conv_load = gcmc_out['components'][comp]["conversion_factor_molec_uc_to_mol_kg"]
            for label in [
                    'loading_absolute_average', 'loading_absolute_dev', 'enthalpy_of_adsorption_average',
                    'enthalpy_of_adsorption_dev'
            ]:
                if label not in out_dict['gcmc']:
                    out_dict['gcmc'][label] = {}
                if comp not in out_dict['gcmc'][label]:
                    out_dict['gcmc'][label][comp] = []
                if label.startswith('loading'):
                    out_dict['gcmc'][label][comp].append(conv_load * gcmc_out['components'][comp][label])
                elif label.startswith('enthalpy'):
                    out_dict['gcmc'][label][comp] = conv_ener * gcmc_out['components'][comp][label]

    return Dict(dict=out_dict)


class MulticompGridWorkChain(WorkChain):
    """Compute multicomponent GCMC in crystalline materials,
    for a mixture of componentes and at specific temperature/pressure conditions.
    """

    @classmethod
    def define(cls, spec):
        super(MulticompGridWorkChain, cls).define(spec)

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
            cls.run_zeopp,
            cls.inspect_zeopp_calc,
            cls.run_raspa_widom,
            cls.inspect_widom_calc,
            if_(cls.should_run_gcmc)(cls.run_raspa_gcmc,),
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
        self.ctx.parameters = aiida_dict_merge(ISOTHERMPARAMETERS_DEFAULT, self.inputs.parameters)
        self.ctx.components = get_components_dict(self.inputs.conditions, self.ctx.parameters)
        self.ctx.ff_params = get_ff_parameters(self.ctx.components, self.ctx.parameters)

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
            self.report("Running zeo++ volpo and block calculation<{}>".format(running.id))
            self.to_context(**{zeopp_label: running})

    def inspect_zeopp_calc(self):
        """Asserts whether all widom calculations are finished ok.
        If so, manage zeopp results.
        NOTE: in this workchain I'm NOT skipping any calculation if the POAV for a component is null.
        """
        for value in self.ctx.components.get_dict().values():
            assert self.ctx["Zeopp_{}".format(value['name'])].is_finished_ok

        for value in self.ctx.components.get_dict().values():
            comp = value['name']
            zeopp_label = "Zeopp_{}".format(comp)
            if self.ctx[zeopp_label].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
                self.out("block_files.{}_block_file".format(comp), self.ctx[zeopp_label].outputs.block)

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
                    True,  # BE CAREFULL: needed in AiiDA-1.0.0 because of github.com/aiidateam/aiida-core/issues/3304
                "Forcefield":
                    "Local",
                "UseChargesFromCIFFile":
                    "yes",
                "CutOff":
                    self.ctx.parameters['ff_cutoff'],
            },
            "System": {
                "framework_1": {
                    "type": "Framework",
                    # "HeliumVoidFraction": self.ctx.geom["POAV_Volume_fraction"], #not computed
                }
            },
            "Component": {},
        }
        mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
        param["System"]["framework_1"]["UnitCells"] = "{} {} {}".format(mult[0], mult[1], mult[2])
        return param

    def run_raspa_widom(self):
        """Run parallel Widom calculation in RASPA, at all temperature specified in the conditions setting."""
        self.ctx.raspa_inputs = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.raspa_inputs['raspa']['framework'] = {"framework_1": self.inputs.structure}
        self.ctx.raspa_inputs['raspa']['file'] = FFBuilder(self.ctx.ff_params)

        for value in self.ctx.components.get_dict().values():
            comp = value['name']
            for i, temp in enumerate(self.inputs.conditions['t_widom']):
                zeopp_label = "Zeopp_{}".format(comp)
                widom_label = "RaspaWidom_{}_{}".format(comp, i)
                self.ctx.raspa_inputs['metadata']['label'] = widom_label
                self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_{}".format(widom_label)
                self.ctx.raspa_param = self._get_widom_param()
                self.ctx.raspa_param["System"]["framework_1"]['ExternalTemperature'] = temp
                self.ctx.raspa_param["Component"][comp] = {}
                self.ctx.raspa_param["Component"][comp]["MoleculeDefinition"] = "Local"
                self.ctx.raspa_param["Component"][comp]["WidomProbability"] = 1.0
                self.ctx.raspa_inputs["raspa"]["block_pocket"] = {}
                if self.ctx[zeopp_label].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
                    self.ctx.raspa_inputs["raspa"]["block_pocket"] = {
                        comp + "_block_file": self.ctx[zeopp_label].outputs.block
                    }
                    self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"] = comp + "_block_file"
                if value['charged']:
                    self.ctx.raspa_param["GeneralSettings"].update({
                        "UseChargesFromCIFFile": "yes",
                        "ChargeMethod": "Ewald",
                        "EwaldPrecision": 1e-6
                    })
                self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)
                running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
                self.report("Running Raspa Widom @ {}K for the Henry coefficient of <{}>".format(temp, comp))
                self.to_context(**{widom_label: running})

    def inspect_widom_calc(self):
        """Asserts whether all widom calculations are finished ok."""
        for value in self.ctx.components.get_dict().values():
            for i, _ in enumerate(self.inputs.conditions['t_widom']):
                comp = value['name']
                widom_label = "RaspaWidom_{}_{}".format(comp, i)
                assert self.ctx[widom_label].is_finished_ok

    def should_run_gcmc(self):  #pylint: disable = no-self-use
        """To implement any decision whether to proceed with GCMC."""
        return True

    def _update_param_input_for_gcmc(self, temp, press):
        """Update Raspa input parameter, from Widom to GCMC"""
        self.ctx.raspa_param["GeneralSettings"].update({
            "NumberOfInitializationCycles": self.ctx.parameters['raspa_gcmc_init_cycles'],
            "NumberOfCycles": self.ctx.parameters['raspa_gcmc_prod_cycles'],
            "PrintPropertiesEvery": int(1e6),
            "PrintEvery": self.ctx.parameters['raspa_gcmc_prod_cycles'] / self.ctx.parameters['raspa_verbosity']
        })
        self.ctx.raspa_param["Component"] = {}
        self.ctx.raspa_inputs["raspa"]["block_pocket"] = {}

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
            if not comp_dict['singlebead']:
                self.ctx.raspa_param["Component"][comp]["RotationProbability"] = 1.0
            if self.ctx[zeopp_label].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
                self.ctx.raspa_param["Component"][comp]["BlockPocketsFileName"] = {"framework_1": comp + "_block_file"}
                self.ctx.raspa_inputs["raspa"]["block_pocket"][comp +
                                                               "_block_file"] = self.ctx[zeopp_label].outputs.block

        self.ctx.raspa_param["System"]["framework_1"]["ExternalPressure"] = press * 1e5  # Pa
        self.ctx.raspa_param["System"]["framework_1"]["ExternalTemperature"] = temp

    def run_raspa_gcmc(self):
        """Summits Raspa GCMC calculation for every condition (i.e, [temp, press] combination)."""

        for i, tp_gcmc in enumerate(self.inputs.conditions['tp_gcmc']):
            gcmc_label = "RaspaGCMC_{}".format(i)
            temp, press = tp_gcmc
            self._update_param_input_for_gcmc(temp, press)
            self.ctx.raspa_inputs['metadata']['label'] = gcmc_label
            self.ctx.raspa_inputs['metadata']['call_link_label'] = "run_{}".format(gcmc_label)

            self.ctx.raspa_inputs['raspa']['parameters'] = Dict(dict=self.ctx.raspa_param)

            running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
            self.report("Running Raspa GCMC #{} @ {} K, {} bar".format(i, temp, press))
            self.to_context(**{gcmc_label: running})

    def return_output_parameters(self):
        """Merge all the parameters into output_parameters, depending on is_porous and is_kh_ehough."""

        all_out_dicts = {}

        for key, val in self.ctx.items():
            if key.startswith('Zeopp_'):
                all_out_dicts[key] = val.outputs.output_parameters
            elif key.startswith('RaspaWidom_'):
                all_out_dicts[key] = val.outputs.output_parameters
            elif key.startswith('RaspaGCMC_'):
                all_out_dicts[key] = val.outputs.output_parameters
        self.out(
            "output_parameters",
            get_output_parameters(inp_conditions=self.inputs.conditions,
                                  components=self.ctx.components,
                                  **all_out_dicts))

        self.report("Workchain completed: output parameters Dict<{}>".format(self.outputs['output_parameters'].pk))
