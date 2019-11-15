# -*- coding: utf-8 -*-
"""Calcfunctions to compute working capacities for different gasses."""

from __future__ import absolute_import
from aiida.engine import calcfunction
from aiida.orm import Dict


@calcfunction
def calc_ch4_working_cap(isot_dict):  # Move this to aiida-lsmo/calcfunctions/calc_working_cap
    """Compute the CH4 working capacity from the output_parameters Dict of IsothermWorkChain.
    This must have run calculations at 5.8 and 65.0 bar (at 298K), which are the standard reference for the evaluation.
    The results can be compared with Simon2015 (10.1039/C4EE03515A).
    """
    from math import sqrt

    out_dict = {}
    out_dict['is_porous'] = isot_dict['is_porous']
    if out_dict['is_porous']:

        ip5 = isot_dict["isotherm"]["pressure"].index(5.8)
        ip65 = isot_dict["isotherm"]["pressure"].index(65.0)

        # in iso_dict uptakes are in mol/kg and enthalpies in kJ/mol
        conv1 = isot_dict["conversion_factor_molec_uc_to_cm3stp_cm3"] / \
                isot_dict["conversion_factor_molec_uc_to_mol_kg"]
        conv2 = isot_dict["conversion_factor_molec_uc_to_gr_gr"] / \
                isot_dict["conversion_factor_molec_uc_to_mol_kg"]

        wc_65bar_average = isot_dict["isotherm"]["loading_absolute_average"][ip65] - isot_dict["isotherm"][
            "loading_absolute_average"][ip5]
        wc_65bar_dev = sqrt(isot_dict["isotherm"]["loading_absolute_dev"][ip5]**2 +
                            isot_dict["isotherm"]["loading_absolute_dev"][ip65]**2)
        wc_65bar_fract = wc_65bar_average / isot_dict["isotherm"]["loading_absolute_average"][ip65]

        out_dict.update({
            "temperature": isot_dict["temperature"],
            "temperature_unit": isot_dict["temperature_unit"],
            "henry_coefficient_average": isot_dict["henry_coefficient_average"],
            "henry_coefficient_dev": isot_dict["henry_coefficient_dev"],
            "henry_coefficient_unit": isot_dict["henry_coefficient_unit"],
            "enthalpy_of_adsorption_5p8bar_average": isot_dict["isotherm"]["enthalpy_of_adsorption_average"][ip5],
            "enthalpy_of_adsorption_5p8bar_dev": isot_dict["isotherm"]["enthalpy_of_adsorption_dev"][ip5],
            "enthalpy_of_adsorption_5p8bar_unit": isot_dict["isotherm"]["enthalpy_of_adsorption_unit"],
            "enthalpy_of_adsorption_65bar_average": isot_dict["isotherm"]["enthalpy_of_adsorption_average"][ip65],
            "enthalpy_of_adsorption_65bar_dev": isot_dict["isotherm"]["enthalpy_of_adsorption_dev"][ip65],
            "enthalpy_of_adsorption_65bar_unit": isot_dict["isotherm"]["enthalpy_of_adsorption_unit"],
            "wc_65bar_cm3stp/cm3_average": wc_65bar_average * conv1,
            "wc_65bar_cm3stp/cm3_dev": wc_65bar_dev * conv1,
            "wc_65bar_cm3stp/cm3_unit": 'cm3 STP/cm3',
            "wc_65bar_g/g_average": wc_65bar_average * conv2,
            "wc_65bar_g/g_dev": wc_65bar_dev * conv2,
            "wc_65bar_g/g_unit": 'g/g',
            "wc_65bar_mol/kg_average": wc_65bar_average,
            "wc_65bar_mol/kg_dev": wc_65bar_dev,
            "wc_65bar_mol/kg_unit": 'mol/kg',
            "wc_65bar_fraction": wc_65bar_fract,
            "wc_65bar_fraction_unit": "-",
        })
    return Dict(dict=out_dict)
