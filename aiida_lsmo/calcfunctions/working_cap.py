# -*- coding: utf-8 -*-
"""Calcfunctions to compute working capacities for different gasses."""

from __future__ import absolute_import
from aiida.engine import calcfunction
from aiida.orm import Dict


def get_molec_uc_to_mg_g(isot_dict):
    """Fix the discrepancy coming from old Raspa calculations, having a typo in the conversion label."""
    if "conversion_factor_molec_uc_to_gr_gr" in isot_dict.get_dict():
        molec_uc_to_mg_g = isot_dict["conversion_factor_molec_uc_to_gr_gr"]
    elif "conversion_factor_molec_uc_to_mg_g" in isot_dict.get_dict():
        molec_uc_to_mg_g = isot_dict["conversion_factor_molec_uc_to_mg_g"]
    return molec_uc_to_mg_g


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

        # conversion factors form mol/kg to cm3STP/cm3 and wt%
        conv1 = isot_dict["conversion_factor_molec_uc_to_cm3stp_cm3"] / isot_dict["conversion_factor_molec_uc_to_mol_kg"]  # pylint: disable=line-too-long
        conv2 = get_molec_uc_to_mg_g(isot_dict) / isot_dict["conversion_factor_molec_uc_to_mol_kg"] / 10

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
            "wc_65bar_wt%_average": wc_65bar_average * conv2,
            "wc_65bar_wt%_dev": wc_65bar_dev * conv2,
            "wc_65bar_wt%_unit": 'g/g/100',
            "wc_65bar_mol/kg_average": wc_65bar_average,
            "wc_65bar_mol/kg_dev": wc_65bar_dev,
            "wc_65bar_mol/kg_unit": 'mol/kg',
            "wc_65bar_fraction": wc_65bar_fract,
            "wc_65bar_fraction_unit": "-",
        })
    return Dict(dict=out_dict)


@calcfunction
def calc_h2_working_cap(isot_dict):  # Move this to aiida-lsmo/calcfunctions/calc_working_cap
    """Compute the H2 working capacity from the output_parameters Dict of MultiTempIsothermWorkChain.
    This must have run calculations at 1, 5 and 100 bar at 77, 198, 298 K.
    The US DOE Target for the Onboard Storage of Hydrogen Vehicles set the bar to 4.5 wt% and 30 g/L (Kapelewski2018).
    Case-A: near-ambient-T adsorption, 100bar/198K to 5bar/298K (cf. Kapelewski2018, 10.1021/acs.chemmater.8b03276)
    ....... Ni2(m-dobdc), experimental: 23.0 g/L
    Case-B: low T adsorption, 100-5bar at 77K (cf. Ahmed2019, 10.1038/s41467-019-09365-w)
    ....... NU-100, best experimental: 35.5 g/L
    Case-C: low T adsorption at low discharge, 100-1bar at 77K (cf. Thornton2017, 10.1021/acs.chemmater.6b04933)
    ....... hypMOF-5059389, best simulated: 40.0 g/L
    """
    from math import sqrt

    out_dict = {}
    out_dict['is_porous'] = isot_dict['is_porous']

    if out_dict['is_porous']:
        press2index = {}
        temp2index = {}
        for press in 1, 5, 100:
            press2index[press] = isot_dict["isotherm"][0]["pressure"].index(press)
        for temp in 77, 198, 298:
            temp2index[temp] = isot_dict["temperature"].index(temp)

        case2pt = {"a": [[100, 198], [5, 298]], "b": [[100, 77], [5, 77]], "c": [[100, 77], [1, 77]]}

        unitconv = {
            "wt%":
                get_molec_uc_to_mg_g(isot_dict) / isot_dict["conversion_factor_molec_uc_to_mol_kg"] /
                10,  # mol/kg to wt%
            "g/L":
                isot_dict["conversion_factor_molec_uc_to_gr_gr"] / isot_dict["conversion_factor_molec_uc_to_mol_kg"] *
                isot_dict["Density"]  # mol/kg to g/L
        }

        out_dict = {}
        for case, presstemp in case2pt.items():
            for unit, conv in unitconv.items():
                load_average = isot_dict["isotherm"][temp2index[presstemp[0][1]]]["loading_absolute_average"][
                    press2index[presstemp[0][0]]]
                disc_average = isot_dict["isotherm"][temp2index[presstemp[1][1]]]["loading_absolute_average"][
                    press2index[presstemp[1][0]]]
                load_dev = isot_dict["isotherm"][temp2index[presstemp[0][1]]]["loading_absolute_dev"][press2index[
                    presstemp[0][0]]]
                disc_dev = isot_dict["isotherm"][temp2index[presstemp[1][1]]]["loading_absolute_dev"][press2index[
                    presstemp[1][0]]]
                out_dict.update({
                    "case-{}_{}_unit".format(case, unit): unit,
                    "case-{}_{}_average".format(case, unit): (load_average - disc_average) * conv,
                    "case-{}_{}_dev".format(case, unit): sqrt(load_dev**2 + disc_dev**2) * conv
                })

    return Dict(dict=out_dict)
