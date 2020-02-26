"""Calculation functions that wrap some advanced script for process evaluation."""

from aiida.engine import calcfunction
from aiida.orm import Dict

PE_PARAMETERS_DEFAULT = { # Parameters used in 10.1021/acscentsci.9b00619
    'gasin': 'coal',
    'vf': 0.35,
    'process': 'TPSA',
    'cp': 985.0,
    'yd': 0.99,
    'eleff': 'carnot',
    'opt': 'PE',
}


@calcfunction
def calc_co2_parasitic_energy(isot_co2, isot_n2, pe_parameters):
    """ Submit calc_pe calculation using AiiDA, for the CO2 parasitic energy.
    :isot_co2: (Dict) CO2 IsothermWorkChainNode.outputs['output_parameters']
    :isot_n2: (Dict) N2 IsothermWorkChainNode.outputs['output_parameters']
    :pe_parameters: (Dict) See PE_PARAMETERS_DEFAULT
    """
    import pandas as pd
    from calc_pe import mainPE

    if not isot_co2['is_porous'] or not isot_n2['is_porous']:
        return Dict(dict={'is_porous': False})

    bar2pa = 1e5  # convert pressure from bar to Pa
    gcm2kgm = 1000  # convert density from g/cm3 to kg/m3

    t_iso = {}
    iso_df = {}
    for i in [0, 1]:
        gas = ['CO_2', 'N_2'][i]
        isot = [isot_co2, isot_n2][i].get_dict()
        t_iso[gas] = isot['temperature']
        iso_df[gas] = pd.DataFrame(columns=['pressure(Pa)', 'loading(mol/kg)', 'HoA(kJ/mol)'])
        iso_df[gas]['pressure(Pa)'] = [p * bar2pa for p in isot['isotherm']["pressure"]]
        iso_df[gas]['loading(mol/kg)'] = isot['isotherm']["loading_absolute_average"]
        iso_df[gas]['HoA(kJ/mol)'] = isot['isotherm']["enthalpy_of_adsorption_average"]
        # TRICK: use the enthalpy from widom (energy-RT) which is more accurate
        #        that the one at 0.001 bar (and which also is NaN for weakly interacting systems)
        iso_df[gas]['HoA(kJ/mol)'].loc[0] = isot['adsorption_energy_widom_average'] - isot['temperature'] / 120.027

    pe_dict = mainPE(
        rho=isot['Density'] * gcm2kgm,
        T_iso=t_iso,  # [T_iso_CO2, T_iso_N2]
        iso_df=iso_df,  # [iso_df_CO2, iso_df_N2] df containing pressure (Pa), loading (mol/kg), HoA (kJ/mol)
        **pe_parameters.get_dict()  # all other parameters that are specific of the modelling
    )
    pe_dict['is_porous'] = True

    return Dict(dict=pe_dict)
