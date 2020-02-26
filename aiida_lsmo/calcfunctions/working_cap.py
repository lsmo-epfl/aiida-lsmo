"""Calcfunctions to compute working capacities for different gasses."""

from math import sqrt
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
def calc_ch4_working_cap(isot_dict):
    """Compute the CH4 working capacity from the output_parameters Dict of IsothermWorkChain.
    This must have run calculations at 5.8 and 65.0 bar (at 298K), which are the standard reference for the evaluation.

    The results can be compared with Simon2015 (10.1039/C4EE03515A).
    """

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
def calc_h2_working_cap(isotmt_dict):
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

    out_dict = {}
    out_dict['is_porous'] = isotmt_dict['is_porous']

    if out_dict['is_porous']:
        press2index = {}
        temp2index = {}
        for press in 1, 5, 100:
            press2index[press] = isotmt_dict["isotherm"][0]["pressure"].index(press)
        for temp in 77, 198, 298:
            temp2index[temp] = isotmt_dict["temperature"].index(temp)

        case2pt = {"a": [[100, 198], [5, 298]], "b": [[100, 77], [5, 77]], "c": [[100, 77], [1, 77]]}

        unitconv = {
            "wt%":  # convert mol/kg to wt%
                get_molec_uc_to_mg_g(isotmt_dict) / isotmt_dict["conversion_factor_molec_uc_to_mol_kg"] / 10,
            "g/L":  # convert mol/kg to g/L
                get_molec_uc_to_mg_g(isotmt_dict) / isotmt_dict["conversion_factor_molec_uc_to_mol_kg"] *
                isotmt_dict["Density"]
        }

        out_dict = {}
        for case, presstemp in case2pt.items():
            for unit, conv in unitconv.items():
                load_average = isotmt_dict["isotherm"][temp2index[presstemp[0][1]]]["loading_absolute_average"][
                    press2index[presstemp[0][0]]]
                disc_average = isotmt_dict["isotherm"][temp2index[presstemp[1][1]]]["loading_absolute_average"][
                    press2index[presstemp[1][0]]]
                load_dev = isotmt_dict["isotherm"][temp2index[presstemp[0][1]]]["loading_absolute_dev"][press2index[
                    presstemp[0][0]]]
                disc_dev = isotmt_dict["isotherm"][temp2index[presstemp[1][1]]]["loading_absolute_dev"][press2index[
                    presstemp[1][0]]]
                out_dict.update({
                    "case-{}_{}_unit".format(case, unit): unit,
                    "case-{}_{}_average".format(case, unit): (load_average - disc_average) * conv,
                    "case-{}_{}_dev".format(case, unit): sqrt(load_dev**2 + disc_dev**2) * conv
                })

    return Dict(dict=out_dict)


@calcfunction
def calc_o2_working_cap(isot_dict):
    """Compute the O2 working capacity from the output_parameters Dict of IsothermWorkChain.
    This must have run calculations at 5 and 140.0 bar (at 298K), to be consistent with the screening of Moghadam2018
    (10.1038/s41467-018-03892-8), for which the MOF ANUGIA (UMCM-152) was found to have a volumetric working capacity
    of 249 vSTP/v (simulations are nearly identical to experiments).
    Consider that, at the same conditions, an empty thank can only store 136 vSTP/v, and a comparable working capacity
    can only br obtained compressing till 300bar.
    """

    out_dict = {}
    out_dict['is_porous'] = isot_dict['is_porous']
    if out_dict['is_porous']:

        ip5 = isot_dict["isotherm"]["pressure"].index(5.0)
        ip140 = isot_dict["isotherm"]["pressure"].index(140.0)

        # conversion factors form mol/kg to cm3STP/cm3 and wt%
        conv1 = isot_dict["conversion_factor_molec_uc_to_cm3stp_cm3"] / isot_dict["conversion_factor_molec_uc_to_mol_kg"]  # pylint: disable=line-too-long
        conv2 = get_molec_uc_to_mg_g(isot_dict) / isot_dict["conversion_factor_molec_uc_to_mol_kg"] / 10

        wc_140bar_average = isot_dict["isotherm"]["loading_absolute_average"][ip140] - isot_dict["isotherm"][
            "loading_absolute_average"][ip5]
        wc_140bar_dev = sqrt(isot_dict["isotherm"]["loading_absolute_dev"][ip5]**2 +
                             isot_dict["isotherm"]["loading_absolute_dev"][ip140]**2)
        wc_140bar_fract = wc_140bar_average / isot_dict["isotherm"]["loading_absolute_average"][ip140]

        out_dict.update({
            "enthalpy_of_adsorption_5bar_average": isot_dict["isotherm"]["enthalpy_of_adsorption_average"][ip5],
            "enthalpy_of_adsorption_5bar_dev": isot_dict["isotherm"]["enthalpy_of_adsorption_dev"][ip5],
            "enthalpy_of_adsorption_5bar_unit": isot_dict["isotherm"]["enthalpy_of_adsorption_unit"],
            "enthalpy_of_adsorption_140bar_average": isot_dict["isotherm"]["enthalpy_of_adsorption_average"][ip140],
            "enthalpy_of_adsorption_140bar_dev": isot_dict["isotherm"]["enthalpy_of_adsorption_dev"][ip140],
            "enthalpy_of_adsorption_140bar_unit": isot_dict["isotherm"]["enthalpy_of_adsorption_unit"],
            "wc_140bar_cm3stp/cm3_average": wc_140bar_average * conv1,
            "wc_140bar_cm3stp/cm3_dev": wc_140bar_dev * conv1,
            "wc_140bar_cm3stp/cm3_unit": 'cm3 STP/cm3',
            "wc_140bar_wt%_average": wc_140bar_average * conv2,
            "wc_140bar_wt%_dev": wc_140bar_dev * conv2,
            "wc_140bar_wt%_unit": 'g/g/100',
            "wc_140bar_mol/kg_average": wc_140bar_average,
            "wc_140bar_mol/kg_dev": wc_140bar_dev,
            "wc_140bar_mol/kg_unit": 'mol/kg',
            "wc_140bar_fraction": wc_140bar_fract,
            "wc_140bar_fraction_unit": "-",
        })
    return Dict(dict=out_dict)
