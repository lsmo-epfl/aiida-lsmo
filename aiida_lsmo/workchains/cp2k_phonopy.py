# -*- coding: utf-8 -*-
"""Cp2kPhonopyWorkChain workchain"""

import io
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.units import CP2KToTHz

from aiida.orm import QueryBuilder, Node
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import append_, while_, WorkChain, ToContext
from aiida_lsmo.utils import aiida_dict_merge

# import sub-workchains
Cp2kBaseWorkChain = WorkflowFactory('cp2k.base')  # pylint: disable=invalid-name

# import aiida data
List = DataFactory('list')
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

        spec.input('structure', valid_type=(CifData, StructureData), required=True, help='Input structure')
        spec.expose_inputs(Cp2kBaseWorkChain, namespace='cp2k_base',
                           exclude=['cp2k.structure', 'cp2k.parameters'
                                   ])  # Using the default parser: no need for lsmo.cp2k_advanced_parser

        spec.outline(
            cls.generate_displacements,
            cls.collect_cp2k_inputs,
            cls.run_cp2k_first,
            while_(cls.should_run_displacement)(cls.run_cp2k_displacement,),
            cls.results,
        )

        spec.output('initial_forces', valid_type=List, required=True, help='Forces computed on the input structure.')

        spec.output('phonopy_params',
                    valid_type=SinglefileData,
                    required=True,
                    help='File phonopy_params.yaml with displacements and forces, to be loaded by Phonopy.')

    def generate_displacements(self):
        """Generate displacements using Phonopy"""

        # CifData to PhonopyAtoms
        ase = self.inputs.structure.get_ase()
        uc = PhonopyAtoms(symbols=ase.get_chemical_symbols(),
                          cell=ase.get_cell(),
                          scaled_positions=ase.get_scaled_positions())

        # Generate displacements
        self.ctx.phonon = Phonopy(
            unitcell=uc,
            supercell_matrix=None,
            primitive_matrix=None,
            factor=CP2KToTHz,
        )
        self.ctx.phonon.generate_displacements(distance=0.01)
        disp_ucs = self.ctx.phonon.supercells_with_displacements
        #phonon.save() to investigate

        # List of PhonopyAtoms to list of StructureData
        # TODO: create an ASE object first and use Structuredata(ase=uc_ase)
        self.ctx.sd_list = []
        for pa in [uc] + disp_ucs:  # 6N+1 in total, first is the original uc
            sd = StructureData(cell=pa.get_cell())
            for i in range(len(pa)):
                symbol = pa.get_chemical_symbols()[i]
                position = pa.get_positions()[i]
                sd.append_atom(symbols=symbol, position=position)
            self.ctx.sd_list.append(sd)

    def collect_cp2k_inputs(self):
        """Collect Cp2k inputs from the last CP2K calculation."""
        qb = QueryBuilder()
        qb.append(Node, filters={f'id': self.inputs.structure.pk}, tag='cif_out')
        qb.append(Node, filters={'attributes.process_label': 'Cp2kCalculation'}, with_descendants='cif_out', tag='calc')
        last_cp2k_calc = qb.distinct().all()[-1][0]
        self.report(f'Found last Cp2kCalculation<{last_cp2k_calc.pk}>: using its inputs.')

        self.ctx.base_inp = AttributeDict(self.exposed_inputs(Cp2kBaseWorkChain, 'cp2k_base'))
        self.ctx.base_inp['cp2k']['settings'] = Dict(dict={'additional_retrieve_list': ['aiida-forces-1_0.xyz']})
        self.ctx.base_inp['cp2k']['parent_calc_folder'] = last_cp2k_calc.outputs.remote_folder
        self.ctx.base_inp['cp2k']['file'] = {}

        for edge_label in last_cp2k_calc.inputs:
            if edge_label.startswith('file'):
                edge_label_file = edge_label.split('__')[-1]
                self.ctx.base_inp['cp2k']['file'][edge_label_file] = last_cp2k_calc.inputs[edge_label]

        param_modify = Dict(
            dict={
                'GLOBAL':
                    {  # NOTE: I'm using the default parser that gets only the final energy: I should try low verbosity and deactivate virtual orbitals'
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
                                'RESTART',  # if parent_calc was cleared, CP2K won't find the WFN and will switch to ATOMIC
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

        self.ctx.base_inp['cp2k']['parameters'] = aiida_dict_merge(last_cp2k_calc.inputs['parameters'], param_modify)

    def run_cp2k_first(self):
        """Run the first CP2K calculation from scratch for the original structure."""
        self.ctx.base_inp['metadata'].update({
            'label': 'cp2k_first',
            'call_link_label': 'run_cp2k_first',
        })
        self.ctx.base_inp['cp2k']['structure'] = self.ctx.sd_list[0]
        self.report('Submit first cp2k')
        running_base = self.submit(Cp2kBaseWorkChain, **self.ctx.base_inp)
        return ToContext(base_wcs=append_(running_base))

    def should_run_displacement(self):
        """Prepare the input for computing displacements, and check if all have been computed"""
        if len(self.ctx.base_wcs) == 1:  # only cp2k_first computed
            cp2k_first_calc = self.ctx.base_wcs[0]
            self.ctx.base_inp['cp2k']['parent_calc_folder'] = cp2k_first_calc.outputs.remote_folder
            self.ctx.index = 1
            return True
        # if len(self.ctx.base_wcs)==3: # for testing purpose
        #     return False
        if len(self.ctx.base_wcs) < len(self.ctx.sd_list):  # still some to compute
            self.ctx.index += 1
            return True
        if len(self.ctx.base_wcs) == len(self.ctx.sd_list):  # all done
            return False

    def run_cp2k_displacement(self):
        """Run the other CP2K calculations for the given structure with displacement."""
        self.ctx.base_inp['metadata'].update({
            'label': f'cp2k_displ_{self.ctx.index}',
            'call_link_label': f'run_cp2k_displ_{self.ctx.index}',
        })
        self.ctx.base_inp['cp2k']['structure'] = self.ctx.sd_list[self.ctx.index]
        self.report(f'Submit displacement {self.ctx.index} of {len(self.ctx.sd_list)-1}')
        running_base = self.submit(Cp2kBaseWorkChain, **self.ctx.base_inp)
        return ToContext(base_wcs=append_(running_base))

    def results(self):
        """Parse forces and store them in a PhonopyYaml file."""

        sets_of_forces = []
        for i, base_wc in enumerate(self.ctx.base_wcs):
            cp2k_calc = base_wc.called[-1]
            forces_path = cp2k_calc.get_retrieved_node()._repository._get_base_folder(
            ).abspath + '/aiida-forces-1_0.xyz'
            with open(forces_path, 'r') as f:
                lines_to_parse = f.readlines()[4:-1]
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

        # NOTE: if we decide to compute force constrants in this workchain, let's monitor the time if it is too expensive
