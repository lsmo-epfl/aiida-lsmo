# -*- coding: utf-8 -*-
"""Utilities related to CP2K."""
from aiida_lsmo.utils import HARTREE2EV


class Cp2kSubsys():
    """Class to represent CP2K subsystem.

     Uses AiiDA StructureData, possibly with additional metadata such as oxidation states.
     """

    def __init__(self, structure, protocol, oxidation_states=None):
        self.structure = structure
        self.oxidation_states = oxidation_states
        self.protocol = protocol

        self._atoms = None
        self._initial_magnetizations = None
        self._kinds = None

    @property
    def atoms(self):
        """Return ASE atoms instance"""
        if not self._atoms:
            self._atoms = self.structure.get_ase()
        return self._atoms

    @property
    def initial_magnetizations(self):
        """ Compute the initial magnetization.

        By default, the initial magnetization is determined based on the element.
        If oxidation states are provided, the guess may be improved based on those.
        """
        from aiida_lsmo.workchains.cp2k_multistage_protocols import INITIAL_MAGNETIZATION, get_default_magnetization
        if not self._initial_magnetizations:
            all_atoms = self.atoms.get_chemical_symbols()

            # Use magnetization table
            initial_magnetizations = [get_default_magnetization(element) for element in all_atoms]

            if self.oxidation_states is None:
                return initial_magnetizations

            # Update selected choices, if oxidation states provided
            for index, element, oxidation_state in zip(self.oxidation_states['metal_indices'],
                                                       self.oxidation_states['metal_symbols'],
                                                       self.oxidation_states['prediction']):
                initial_magnetizations[index] = INITIAL_MAGNETIZATION[element]['magnetization'][oxidation_state]

            self._initial_magnetizations = initial_magnetizations

        return self._initial_magnetizations

    def get_multiplicity_section(self):
        """ Compute the total multiplicity of the structure by summing the atomic magnetizations.

            multiplicity = 1 + sum_i ( natoms_i * magnetization_i ), for each atom_type i
                         = 1 + sum_i magnetization_j, for each atomic site j

        :returns: dict (for cp2k input)
        """
        multiplicity = 1 + sum(self.initial_magnetizations)
        multiplicity = int(round(multiplicity))

        multiplicity_dict = {'FORCE_EVAL': {'DFT': {'MULTIPLICITY': multiplicity}}}
        if multiplicity != 1:
            multiplicity_dict['FORCE_EVAL']['DFT']['UKS'] = True
        return multiplicity_dict

    def get_kinds_section(self):
        """ Write the &KIND sections given the structure and the settings_dict"""
        from aiida_lsmo.workchains.cp2k_multistage_protocols import get_default_magnetization

        if not self._kinds:
            kinds = []
            all_atoms = set(self.atoms.get_chemical_symbols())
            for atom in sorted(all_atoms):
                kinds.append({
                    '_': atom,
                    'BASIS_SET': self.protocol['basis_set'][atom],
                    'POTENTIAL': self.protocol['pseudopotential'][atom],
                    'MAGNETIZATION': get_default_magnetization(atom),
                })
        else:
            kinds = [{
                '_': k['kind'],
                'ELEMENT': k['element'],
                'BASIS_SET': self.protocol['basis_set'][k['element']],
                'POTENTIAL': self.protocol['pseudopotential'][k['element']],
                'MAGNETIZATION': k['magnetization'],
            } for k in self._kinds]

        return {'FORCE_EVAL': {'SUBSYS': {'KIND': kinds}}}

    def tag_atoms(self):
        """Tag atoms according to initial magnetization and store new kinds.

        For each metal atom species that occurs with more than one oxidation state, add tags (starting from 1) to
        differentiate them. Resulting kind names will be e.g. 'Cu1', 'Cu2'.

        Stores list of atom kinds in self._kinds (used by `get_kinds_section`).
        """
        all_atoms = self.atoms.get_chemical_symbols()

        n_sites = len(all_atoms)
        symbols = sorted(set(all_atoms))

        # Prepare dictionary [symbol][magnetization] => tag  (tag is 0, if only single oxidation state)
        element_magnetization_tag_map = {}
        kinds = []
        for symbol in symbols:
            magnetizations = list({self.initial_magnetizations[i] for i in range(n_sites) if all_atoms[i] == symbol})
            magnetization_tag_map = {}

            if len(magnetizations) == 1:
                kinds.append({'element': symbol, 'kind': symbol, 'magnetization': magnetizations[0]})
                magnetization_tag_map[magnetizations[0]] = None
            else:
                for index, magnetization in enumerate(magnetizations):
                    kinds.append({
                        'element': symbol,
                        'kind': f'{symbol}{str(index + 1)}',
                        'magnetization': magnetization
                    })
                    magnetization_tag_map[magnetization] = index + 1  # start from 1
            element_magnetization_tag_map[symbol] = magnetization_tag_map
        self._kinds = kinds

        # Set tags for metal atoms where needed
        for element, magnetization, atom in zip(all_atoms, self.initial_magnetizations, self.atoms):
            tag = element_magnetization_tag_map[element][magnetization]
            if tag is not None:  # only tag when needed
                atom.tag = tag

    def get_kinds_with_ghost_section(self):
        """Write the &KIND sections given the structure and the settings_dict, and add also GHOST atoms"""
        from aiida_lsmo.workchains.cp2k_multistage_protocols import get_default_magnetization
        kinds = []
        all_atoms = set(self.atoms.get_chemical_symbols())
        for atom in sorted(all_atoms):
            kinds.append({
                '_': atom,
                'BASIS_SET': self.protocol['basis_set'][atom],
                'POTENTIAL': self.protocol['pseudopotential'][atom],
                'MAGNETIZATION': get_default_magnetization(atom),
            })
            kinds.append({'_': atom + '_ghost', 'BASIS_SET': self.protocol['basis_set'][atom], 'GHOST': True})
        return {'FORCE_EVAL': {'SUBSYS': {'KIND': kinds}}}


def get_bsse_section(natoms_a, natoms_b, mult_a=1, mult_b=1, charge_a=0, charge_b=0):  # pylint: disable=too-many-arguments
    """Get the &FORCE_EVAL/&BSSE section."""
    bsse_section = {
        'FORCE_EVAL': {
            'BSSE' : {
                'FRAGMENT': [{
                    'LIST': '1..{}'.format(natoms_a)
                },
                {
                    'LIST': '{}..{}'.format(natoms_a + 1, natoms_a + natoms_b)
                }],
                'CONFIGURATION': [
                    { # A fragment with basis set A
                        'MULTIPLICITY': mult_a,
                        'CHARGE': charge_a,
                        'GLB_CONF': '1 0',
                        'SUB_CONF': '1 0',
                        },
                    { # B fragment with basis set B
                        'MULTIPLICITY': mult_b,
                        'CHARGE': charge_b,
                        'GLB_CONF': '0 1',
                        'SUB_CONF': '0 1',
                        },
                    { # A fragment with basis set A+B
                        'MULTIPLICITY': mult_a,
                        'CHARGE': charge_a,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '1 0',
                        },
                    { # B fragment with basis set A+B
                        'MULTIPLICITY': mult_b,
                        'CHARGE': charge_b,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '0 1',
                        },
                    { # A+B fragments with basis set A+B
                        'MULTIPLICITY': mult_a + mult_b - 1,
                        'CHARGE': charge_a + charge_b,
                        'GLB_CONF': '1 1',
                        'SUB_CONF': '1 1',
                        }
                ]
            }
        }
    }
    return bsse_section


# Functions to parse results


def ot_has_small_bandgap(cp2k_input, cp2k_output, bandgap_thr_ev):
    """ Returns True if the calculation used OT and had a smaller bandgap then the guess needed for the OT.
    (NOTE: It has been observed also negative bandgap with OT in CP2K!)
    cp2k_input: dict
    cp2k_output: dict
    bandgap_thr_ev: float [eV]
    """
    list_true = [True, 'T', 't', '.TRUE.', 'True', 'true']  #add more?
    try:
        ot_settings = cp2k_input['FORCE_EVAL']['DFT']['SCF']['OT']
        if '_' not in ot_settings.keys() or ot_settings['_'] in list_true:  #pylint: disable=simplifiable-if-statement
            using_ot = True
        else:
            using_ot = False
    except KeyError:
        using_ot = False
    min_bandgap_ev = min(cp2k_output['bandgap_spin1_au'], cp2k_output['bandgap_spin2_au']) * HARTREE2EV
    is_bandgap_small = (min_bandgap_ev < bandgap_thr_ev)
    return using_ot and is_bandgap_small
