"""Calcfunctions to compute gas-selectivity related applications."""

from aiida.engine import calcfunction
from aiida.orm import Dict


@calcfunction
def calc_selectivity(isot_dict_a, isot_dict_b):
    """Compute the selectivity of gas A on gas B as S = kH_a/kH_b.
    Note that if the material is not porous to one of the materials, the result is simply {'is_porous': False}.
    To maintain the comptaibility with v1, intead of checking 'is_porous', it checks for the henry_coefficient_average
    key in the Dict.
    """
    from math import sqrt

    out_dict = {}
    out_dict['is_porous'] = all([
        "henry_coefficient_average" in isot_dict_a.get_dict() and "henry_coefficient_average" in isot_dict_b.get_dict()
    ])

    if out_dict['is_porous']:
        out_dict[
            "selectivity_average"] = isot_dict_a["henry_coefficient_average"] / isot_dict_b["henry_coefficient_average"]
        out_dict["selectivity_dev"] = out_dict["selectivity_average"] * sqrt(
            (isot_dict_a["henry_coefficient_dev"] / isot_dict_a["henry_coefficient_average"]) +
            (isot_dict_b["henry_coefficient_dev"] / isot_dict_b["henry_coefficient_average"]))
    return Dict(dict=out_dict)
