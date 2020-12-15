# -*- coding: utf-8 -*-
"""CalcFunction to compute the oxidation states of metals using oximachine"""
from aiida.engine import calcfunction
from aiida.orm import Dict

try:
    import oximachinerunner
except ImportError as exc:
    raise ImportError("Please install the 'oximachine' extra for oxidation state prediction.") from exc

# Note: This loads the model into memory
OXIMACHINE_RUNNER = oximachinerunner.OximachineRunner()

# def compute_oxidation_states(structure):
#     """compute the oxidation states of metals using oximachine"""
#     return OXIMACHINE_RUNNER.run_oximachine(structure)


@calcfunction
def compute_oxidation_states(cif):
    """Compute the oxidation states of metals using oximachine

    :param cif:  AiiDA CifData instance
    :return: AiiDA Dict node
    """
    results_dict = OXIMACHINE_RUNNER.run_oximachine(cif.get_ase())
    return Dict(dict=results_dict)
