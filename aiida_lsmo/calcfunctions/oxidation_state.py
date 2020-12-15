# -*- coding: utf-8 -*-
"""CalcFunction to compute the oxidation states of metals using oximachine"""
from aiida.engine import calcfunction
from aiida.orm import Dict, StructureData

try:
    import oximachinerunner
except ImportError as exc:
    raise ImportError("Please install the 'oximachine' extra for oxidation state prediction.") from exc

# Note: This loads the model into memory
OXIMACHINE_RUNNER = oximachinerunner.OximachineRunner(modelname='mof')


@calcfunction
def compute_oxidation_states(cif):
    """Compute the oxidation states of metals using oximachine

    :param cif:  AiiDA CifData instance
    :return: AiiDA Dict node
    """
    results_dict = OXIMACHINE_RUNNER.run_oximachine(cif.get_ase())
    results_dict['oximachine_version'] = str(OXIMACHINE_RUNNER)
    return Dict(dict=results_dict)


def tag_ase_atoms(atoms, oxi_results):
    """Tag metal atoms depending on oxidation state.

    For each metal atom species that occurs with different oxidation states, add a tag (starting from 1).

    :param atoms:  ASE Atoms instance
    :returns: ASE Atoms instance with tags and charges
    """
    n_sites = len(oxi_results['prediction'])
    symbols = set(oxi_results['metal_symbols'])

    # Prepare dictionary [symbol][state] => tag  (tag is 0, if only single oxidation state)
    symbol_state_tag_map = {}
    for symbol in symbols:
        oxidation_states = list(
            {oxi_results['prediction'][i] for i in range(n_sites) if oxi_results['metal_symbols'][i] == symbol})

        state_tag_dict = {}
        if len(oxidation_states) == 1:
            state_tag_dict[oxidation_states[0]] = None
        else:
            for index, state in enumerate(oxidation_states):
                state_tag_dict[state] = index + 1  # start from 1
        symbol_state_tag_map[symbol] = state_tag_dict

    # Set tags for metal atoms where needed
    for metal_index, atom_index in enumerate(oxi_results['metal_indices']):
        symbol = oxi_results['metal_symbols'][metal_index]
        state = oxi_results['prediction'][metal_index]
        atoms[atom_index].charge = -state  # charge is not currently used; perhaps useful later

        tag = symbol_state_tag_map[symbol][state]
        if tag is not None:
            atoms[atom_index].tag = tag

    return atoms


@calcfunction
def tag_structure_data(structure, oxi_results):
    """Tag StructureData depending on oxidation state.

    For use with aiida-cp2k.

    :param structure: AiiDA StructureData node
    :param oxi_results: AiiDA Dict node with results from oxidation state predictions
    """
    return StructureData(ase=tag_ase_atoms(structure.get_ase(), oxi_results.get_dict()))
