# -*- coding: utf-8 -*-
"""Cp2kPhonopyWorkChain workchain"""

import io
import ase
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.units import CP2KToTHz

from aiida.orm import QueryBuilder, Node, load_node
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import while_, WorkChain
from aiida_lsmo.utils import aiida_dict_merge

# import sub-workchains
Cp2kBaseWorkChain = WorkflowFactory('cp2k.base')  # pylint: disable=invalid-name

# import aiida data
Str = DataFactory('str')  # pylint: disable=invalid-name
Int = DataFactory('int')  # pylint: disable=invalid-name
List = DataFactory('list')  # pylint: disable=invalid-name
Dict = DataFactory('dict')  # pylint: disable=invalid-name
CifData = DataFactory('cif')  # pylint: disable=invalid-name
StructureData = DataFactory('structure')  # pylint: disable=invalid-name
SinglefileData = DataFactory('singlefile')  # pylint: disable=invalid-name


class Cp2kPhonopyWorkChain(WorkChain):
    """A workchain to compute phonon frequencies using CP2K and Phonopy"""

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super().define(spec)

        spec.input(
            'structure',
            valid_type=StructureData,  # not a CifData because it needs to have atom tags for the oxidation states
            required=True,
            help='Input structure, output of some Cp2kCalculation.')
        spec.input(
            'cp2kcalc',
            valid_type=Str,
            required=False,
            help='Provide the UUID of a specific Cp2kCalc to be used. If not provided the WC search for an ancestor.')
        spec.input('mode',
                   valid_type=Str,
                   required=False,
                   default=lambda: Str('serial'),
                   help='Mode of the calculation: "serial" (default) or "parallel".')
        spec.input('max_displacements',
                   valid_type=Int,
                   required=False,
                   default=lambda: Int(10000),
                   help='Set a maximum number of displacements (or zero) for testing purpose.')
        spec.expose_inputs(Cp2kBaseWorkChain, namespace='cp2k_base',
                           exclude=['cp2k.structure', 'cp2k.parameters'
                                   ])  # Using the default parser: no need for lsmo.cp2k_advanced_parser

        spec.outline(
            cls.collect_cp2k_inputs,
            cls.generate_displacements,
            cls.run_cp2k_first,
            while_(cls.should_run_displacement)(cls.run_cp2k_displacement,),
            cls.results,
        )

        spec.output('initial_forces', valid_type=List, required=True, help='Forces computed on the input structure.')

        spec.output('phonopy_params',
                    valid_type=SinglefileData,
                    required=True,
                    help='File phonopy_params.yaml with displacements and forces, to be loaded by Phonopy.')

    def collect_cp2k_inputs(self):
        """Collect Cp2k inputs from the reference CP2K calculation."""
        # Use the input Cp2kCal if provided, or look for some ancestor calculation of the input structure
        if 'cp2kcalc' in self.inputs:
            ref_cp2k_calc = load_node(uuid=self.inputs.cp2kcalc.value)
        else:
            qb_cp2k_calc = QueryBuilder()
            qb_cp2k_calc.append(Node, filters={'id': self.inputs.structure.pk}, tag='cif_out')
            qb_cp2k_calc.append(Node,
                                filters={'attributes.process_label': 'Cp2kCalculation'},
                                with_descendants='cif_out',
                                tag='calc')
            ref_cp2k_calc = qb_cp2k_calc.distinct().all()[-1][0]

        self.report(f'Using inputs from Cp2kCalculation<{ref_cp2k_calc.pk}>, StructureData<{self.inputs.structure.pk}>')

        # Collect the other inputs for Cp2kBaseWC from the ref_cp2k_calc to self.ctx.base_inp
        self.ctx.base_inp = AttributeDict(self.exposed_inputs(Cp2kBaseWorkChain, 'cp2k_base'))
        self.ctx.base_inp['cp2k']['settings'] = Dict(dict={'additional_retrieve_list': ['aiida-forces-1_0.xyz']})
        self.ctx.base_inp['cp2k']['parent_calc_folder'] = ref_cp2k_calc.outputs.remote_folder
        self.ctx.base_inp['cp2k']['file'] = {}
        if 'file' in ref_cp2k_calc.inputs:
            for edge_label, file_node in ref_cp2k_calc.inputs.file.items():
                self.ctx.base_inp['cp2k']['file'][edge_label] = file_node

        param_modify = Dict(
            dict={
                'GLOBAL': {  # NOTE: I'm using the default parser that gets only the final energy: lowering verbosity
                    'PRINT_LEVEL': 'LOW',
                    'RUN_TYPE': 'ENERGY_FORCE'
                },
                'FORCE_EVAL': {
                    'STRESS_TENSOR': 'ANALYTICAL',
                    'PRINT': {
                        'FORCES': {
                            'FILENAME': 'forces'
                        }
                    },
                    'DFT': {
                        'WFN_RESTART_FILE_NAME': './parent_calc/aiida-RESTART.wfn',
                        'SCF': {
                            'SCF_GUESS':
                                'RESTART',  # if parent_calc was cleared, CP2K won't find the WFN and will set ATOMIC
                        },
                        'PRINT': {
                            'E_DENSITY_CUBE': {
                                '_': 'OFF'
                            },
                            'MO_CUBES': {
                                '_': 'OFF'
                            },
                            'MULLIKEN': {
                                '_': 'OFF'
                            }
                        }
                    }
                }
            }).store()

        self.ctx.base_inp['cp2k']['parameters'] = aiida_dict_merge(ref_cp2k_calc.inputs['parameters'], param_modify)

    def generate_displacements(self):
        """Generate displacements using Phonopy"""

        # CifData to PhonopyAtoms
        uc_opt_ase = self.inputs.structure.get_ase()
        uc_opt_pa = PhonopyAtoms(symbols=uc_opt_ase.get_chemical_symbols(),
                                 cell=uc_opt_ase.get_cell(),
                                 scaled_positions=uc_opt_ase.get_scaled_positions())

        # Generate displacements
        self.ctx.phonon = Phonopy(
            unitcell=uc_opt_pa,
            supercell_matrix=None,
            primitive_matrix=None,
            factor=CP2KToTHz,
        )
        self.ctx.phonon.generate_displacements(distance=0.01)
        uc_displ_pa_list = self.ctx.phonon.supercells_with_displacements

        # List of PhonopyAtoms to list of StructureData
        # The list will contain the reference optimized structure plus the 6N displacements.
        # For testing purpose one can use the input "max_displacements" to compute only some of the 6N displacements.
        self.ctx.sd_list = []
        for uc_pa in [uc_opt_pa] + uc_displ_pa_list[:self.inputs.max_displacements.value]:
            uc_ase = ase.Atoms(
                cell=uc_opt_pa.get_cell(),  # Same for all
                symbols=uc_opt_pa.get_chemical_symbols(),  # Same for all
                tags=uc_opt_ase.get_tags(),  # Important for atoms of the same element with different oxidation state
                positions=uc_pa.get_positions(),  # Different for each displacement
            )
            uc_sd = StructureData(ase=uc_ase)
            uc_sd.store()
            self.ctx.sd_list.append(uc_sd)
        self.ctx.ndisplacements = len(self.ctx.sd_list) - 1  # = s6N (unless inputs.max_displacements < 6N)

    def run_cp2k_first(self):
        """Run the first CP2K calculation from scratch for the original structure."""
        self.ctx.base_inp['metadata'].update({
            'label': 'cp2k_first',
            'call_link_label': 'run_cp2k_first',
        })
        self.ctx.base_inp['cp2k']['structure'] = self.ctx.sd_list[0]
        self.report('Submit first cp2k')
        running_base = self.submit(Cp2kBaseWorkChain, **self.ctx.base_inp)
        self.ctx.index = 0
        self.to_context(**{f'cp2kbase_{self.ctx.index}': running_base})

    def should_run_displacement(self):
        """Prepare the input for computing displacements, and check if all have been computed"""
        if self.ctx.index == 0:  # only cp2k_first computed
            self.ctx.base_inp['cp2k']['parent_calc_folder'] = self.ctx.cp2kbase_0.outputs.remote_folder
            self.ctx.index = 1
            return True
        if self.ctx.index < self.ctx.ndisplacements:  # still some to compute
            self.ctx.index += 1
            return True
        # if self.ctx.index == self.ctx.ndisplacements:  # all done: this will also be tweeked after parallel calcs!
        return False

    def run_cp2k_displacement(self):
        """Run the other CP2K calculations for the given structure with displacement."""
        if self.inputs.mode == 'serial':
            index_list = [self.ctx.index]
        elif self.inputs.mode == 'parallel':
            index_list = list(range(1, self.ctx.ndisplacements + 1))

        for index in index_list:
            self.ctx.index = index
            self.ctx.base_inp['metadata'].update({
                'label': f'cp2k_displ_{self.ctx.index}',
                'call_link_label': f'run_cp2k_displ_{self.ctx.index}',
            })
            self.ctx.base_inp['cp2k']['structure'] = self.ctx.sd_list[self.ctx.index]
            self.report(f'Submit displacement {self.ctx.index} of {len(self.ctx.sd_list)-1}')
            running_base = self.submit(Cp2kBaseWorkChain, **self.ctx.base_inp)
            self.to_context(**{f'cp2kbase_{self.ctx.index}': running_base})

    def results(self):
        """Parse forces and store them in a PhonopyYaml file."""

        sets_of_forces = []
        for i in range(0, self.ctx.ndisplacements + 1):
            cp2k_calc = self.ctx[f'cp2kbase_{i}'].called[-1]  # Get the Cp2kCalc from Cp2kBaseWC
            forces_path = cp2k_calc.get_retrieved_node()._repository._get_base_folder(  # pylint: disable=protected-access
            ).abspath + '/aiida-forces-1_0.xyz'
            with open(forces_path, 'r') as forces_file:
                lines_to_parse = forces_file.readlines()[4:-1]
            forces_parsed = [list(map(float, l.split()[3:])) for l in lines_to_parse]
            sets_of_forces.append(forces_parsed)

        # Output the forces of the initial structure as a List
        initial_forces = List(list=sets_of_forces[0])
        initial_forces.store()
        self.out('initial_forces', initial_forces)
        self.report(f'Forces computed on the input structure: List<{initial_forces.pk}>')

        # Output the forces and displacements as a YAML SinglefileData
        self.ctx.phonon.set_forces(sets_of_forces[1:])  # Exclude forces on the input structure
        phpy_yaml = PhonopyYaml()
        phpy_yaml.set_phonon_info(self.ctx.phonon)
        phpy_yaml_bytes = bytes(str(phpy_yaml), encoding='utf-8')
        phpy_yaml_sfd = SinglefileData(file=io.BytesIO(phpy_yaml_bytes), filename='phonopy_params.yaml')
        phpy_yaml_sfd.store()
        self.out('phonopy_params', phpy_yaml_sfd)
        self.report(f'Output phonopy_params.yaml: SinglefileData<{phpy_yaml_sfd.pk}>')
