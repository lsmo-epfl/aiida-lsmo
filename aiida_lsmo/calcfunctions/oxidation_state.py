# -*- coding: utf-8 -*-
"""CalcFunction to compute the oxidation states of metals using oximachine"""
from aiida.engine import calcfunction
from aiida.orm import Dict
import oximachinerunner

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
