# -*- coding: utf-8 -*-
"""A work chain."""
import functools

# AiiDA modules
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.orm import Str, Dict
from aiida.engine import calcfunction
from aiida.engine import WorkChain, if_
from aiida_lsmo.utils import check_resize_unit_cell

from .isotherm import get_molecule_dict, get_ff_parameters, get_atomic_radii, validate_dict
from .parameters_schemas import FF_PARAMETERS_VALIDATOR, NUMBER, Required

RaspaBaseWorkChain = WorkflowFactory('raspa.base')  #pylint: disable=invalid-name

# Defining DataFactory, CalculationFactory and default parameters
CifData = DataFactory('core.cif')  #pylint: disable=invalid-name
ZeoppParameters = DataFactory('zeopp.parameters')  #pylint: disable=invalid-name

ZeoppCalculation = CalculationFactory('zeopp.network')  #pylint: disable=invalid-name
FFBuilder = CalculationFactory('lsmo.ff_builder')  # pylint: disable=invalid-name


@calcfunction
def get_zeopp_parameters(molecule_dict, isotparam):
    """Get the ZeoppParameters from the inputs of the workchain"""
    probe_rad = molecule_dict['proberad'] * isotparam['zeopp_probe_scaling']
    param_dict = {
        'ha': 'DEF',
        'block': [probe_rad, isotparam['zeopp_block_samples']],
    }
    return ZeoppParameters(dict=param_dict)


@calcfunction
def get_output_parameters(inp_parameters, **all_out_dicts):
    """Extract results to output_parameters Dict."""

    out_dict = {
        'temperatures': inp_parameters['temperatures'],
        'temperatures_unit': 'K',
        'henry_coefficient_unit': 'mol/kg/Pa',
        'adsorption_energy_widom_unit': 'kJ/mol',
        'widom_rosenbluth_factor_unit': '-',
    }

    widom_keys = [
        'henry_coefficient_average', 'henry_coefficient_dev', 'adsorption_energy_widom_average',
        'adsorption_energy_widom_dev', 'widom_rosenbluth_factor_average', 'widom_rosenbluth_factor_dev'
    ]

    for temp in inp_parameters['temperatures']:
        widom_label = f'RaspaWidom_{int(temp)}'
        output_widom = all_out_dicts[widom_label].get_dict()
        for key in widom_keys:
            if key not in out_dict:
                out_dict[key] = []
            system_key = list(output_widom.keys())[0]
            comp_key = list(output_widom[system_key]['components'].keys())[0]
            out_dict[key].append(output_widom[system_key]['components'][comp_key][key])

    return Dict(out_dict)


class SinglecompWidomWorkChain(WorkChain):
    """Computes widom insertion for a framework/box at different temperatures."""

    parameters_schema = FF_PARAMETERS_VALIDATOR.extend({
        Required('zeopp_probe_scaling', default=1.0, description="scaling probe's diameter: molecular_rad * scaling"):
            NUMBER,
        Required('zeopp_block_samples',
                 default=int(100),
                 description='Number of samples for BLOCK calculation (per A^3).'):
            int,
        Required('raspa_verbosity', default=10, description='Print stats every: number of cycles / raspa_verbosity.'):
            int,
        Required('raspa_widom_cycles', default=int(1e5), description='Number of Widom cycles.'):
            int,
        Required('temperatures', default=[300, 400]): [NUMBER],
    })
    parameters_info = parameters_schema.schema  # shorthand for printing

    @classmethod
    def define(cls, spec):
        super(SinglecompWidomWorkChain, cls).define(spec)

        spec.expose_inputs(ZeoppCalculation, namespace='zeopp', include=['code', 'metadata'])
        spec.expose_inputs(RaspaBaseWorkChain, namespace='raspa_base', exclude=['raspa.structure', 'raspa.parameters'])
        spec.input('structure', valid_type=CifData, required=False, help='Adsorbent framework CIF or None for a box.')
        spec.input('molecule',
                   valid_type=(Str, Dict, CifData),
                   help='Adsorbate molecule: settings to be read from the yaml.' +
                   'Advanced: input a Dict for non-standard settings.')
        spec.input('parameters',
                   valid_type=Dict,
                   validator=functools.partial(validate_dict, schema=cls.parameters_schema),
                   help='Main parameters and settings for the calculations, to overwrite PARAMETERS_DEFAULT.')
        spec.outline(
            cls.setup,
            if_(cls.should_run_zeopp)(cls.run_zeopp, cls.inspect_zeopp_calc),
            cls.run_raspa_widom,
            cls.return_output_parameters,
        )

        spec.output('output_parameters', valid_type=Dict, required=True, help='Main results of the work chain.')

        spec.expose_outputs(ZeoppCalculation, include=['block'])

    def setup(self):
        """Initialize parameters"""
        self.ctx.sim_in_box = 'structure' not in self.inputs.keys()

        # Get the parameters Dict, merging defaults with user settings
        @calcfunction
        def get_valid_dict(dict_node):
            return Dict(self.parameters_schema(dict_node.get_dict()))

        self.ctx.parameters = get_valid_dict(self.inputs.parameters)
        if isinstance(self.inputs.molecule, Str):
            self.ctx.molecule = get_molecule_dict(self.inputs.molecule)
            self.ctx.cif_molecule = None
        elif isinstance(self.inputs.molecule, Dict):
            self.ctx.molecule = self.inputs.molecule
            self.ctx.cif_molecule = None
        elif isinstance(self.inputs.molecule, CifData):
            self.ctx.molecule = Dict(
                dict={
                    'name': 'MOL',
                    'forcefield': 'on-the-fly',
                    'molsatdens': 99.9,  # experm @Tb
                    'proberad': 9.9,  # sigma/2
                    'singlebead': False,
                    'charged': True,
                }).store()
            self.ctx.cif_molecule = self.inputs.molecule
        self.ctx.ff_params = get_ff_parameters(self.ctx.molecule, self.ctx.parameters)

    def should_run_zeopp(self):
        """Return if it should run zeopp calculation."""
        return not self.ctx.sim_in_box and self.ctx.parameters['zeopp_probe_scaling'] > 0.0

    def run_zeopp(self):
        """It performs the full zeopp calculation for all components."""
        zeopp_inputs = self.exposed_inputs(ZeoppCalculation, 'zeopp')
        zeopp_inputs.update({'structure': self.inputs.structure, 'atomic_radii': get_atomic_radii(self.ctx.parameters)})
        zeopp_inputs['metadata']['label'] = 'ZeoppBlock'
        zeopp_inputs['metadata']['call_link_label'] = 'run_zeopp_block'
        zeopp_inputs['parameters'] = get_zeopp_parameters(self.ctx.molecule, self.ctx.parameters)
        running = self.submit(ZeoppCalculation, **zeopp_inputs)
        self.report('Running zeo++ block calculation<{}>'.format(running.id))
        self.to_context(**{'ZeoppBlock': running})

    def inspect_zeopp_calc(self):
        """Asserts whether all widom calculations are finished ok and expose block file."""
        assert self.ctx['ZeoppBlock'].is_finished_ok

        if self.ctx['ZeoppBlock'].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
            self.out_many(self.exposed_outputs(self.ctx.zeopp, ZeoppCalculation))

    def _get_widom_inputs(self):
        """Generate Raspa input parameters from scratch, for a Widom calculation."""
        self.ctx.raspa_inputs = self.exposed_inputs(RaspaBaseWorkChain, 'raspa_base')
        self.ctx.raspa_inputs['raspa']['file'] = FFBuilder(self.ctx.ff_params, cif_molecule=self.ctx.cif_molecule)
        self.ctx.raspa_inputs['raspa']['block_pocket'] = {}
        verbosity = self.ctx.parameters['raspa_widom_cycles'] / self.ctx.parameters['raspa_verbosity']
        self.ctx.raspa_param = {
            'GeneralSettings': {
                'SimulationType': 'MonteCarlo',
                'NumberOfInitializationCycles': 0,
                'NumberOfCycles': self.ctx.parameters['raspa_widom_cycles'],
                'PrintPropertiesEvery': verbosity,
                'PrintEvery': int(1e10),
                'RemoveAtomNumberCodeFromLabel': True,  # github.com/aiidateam/aiida-core/issues/3304
                'Forcefield': 'Local',
                'UseChargesFromCIFFile': 'yes',
                'CutOff': self.ctx.parameters['ff_cutoff'],
                'ChargeMethod': 'Ewald' if self.ctx.molecule['charged'] else 'None'
            },
            'System': {},
            'Component': {
                self.ctx.molecule['name']: {
                    'MoleculeDefinition': 'Local',
                    'WidomProbability': 1.0
                }
            }
        }

        if self.ctx.sim_in_box:
            self.ctx.raspa_param['System'] = {
                'box_1': {
                    'type': 'Box',
                    'BoxLengths': '100 100 100'  #large to avoid Coulomb self interaction
                }
            }
        else:
            mult = check_resize_unit_cell(self.inputs.structure, 2 * self.ctx.parameters['ff_cutoff'])
            self.ctx.raspa_param['System'] = {
                'framework_1': {
                    'type': 'Framework',
                    'UnitCells': f'{mult[0]} {mult[1]} {mult[2]}'
                }
            }
            self.ctx.raspa_inputs['raspa']['framework'] = {'framework_1': self.inputs.structure}

        if 'ZeoppBlock' in self.ctx and \
           self.ctx['ZeoppBlock'].outputs.output_parameters['Number_of_blocking_spheres'] > 0:
            self.ctx.raspa_param['Component'][self.ctx.molecule['name']]['BlockPocketsFileName'] = 'block_file'
            self.ctx.inp['raspa']['block_pocket'] = {'block_file': self.ctx['ZeoppBlock'].outputs.block}

    def run_raspa_widom(self):
        """Run parallel Widom calculation in RASPA, at all temperature specified in the conditions setting."""

        self._get_widom_inputs()

        for temp in self.inputs.parameters['temperatures']:
            widom_label = f'RaspaWidom_{int(temp)}'
            self.ctx.raspa_inputs['metadata']['label'] = widom_label
            self.ctx.raspa_inputs['metadata']['call_link_label'] = f'run_{widom_label}'
            system_key = list(self.ctx.raspa_param['System'].keys())[0]
            self.ctx.raspa_param['System'][system_key]['ExternalTemperature'] = temp
            self.ctx.raspa_inputs['raspa']['parameters'] = Dict(self.ctx.raspa_param)
            running = self.submit(RaspaBaseWorkChain, **self.ctx.raspa_inputs)
            self.report(f"Running Raspa Widom @ {temp}K for the Henry coefficient of <{self.ctx.molecule['name']}>")
            self.to_context(**{widom_label: running})

    def return_output_parameters(self):
        """Merge all the parameters into output_parameters, depending on is_porous and is_kh_ehough."""

        all_out_dicts = {}

        for key, val in self.ctx.items():
            if key.startswith('RaspaWidom_'):
                all_out_dicts[key] = val.outputs.output_parameters
        self.out('output_parameters', get_output_parameters(inp_parameters=self.ctx.parameters, **all_out_dicts))

        self.report('Workchain completed: output parameters Dict<{}>'.format(self.outputs['output_parameters'].pk))
