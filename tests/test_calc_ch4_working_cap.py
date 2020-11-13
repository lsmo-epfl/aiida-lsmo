# -*- coding: utf-8 -*-
"""Test for calc_ch4_working_cap CalcFunction"""
# pylint: disable=invalid-name
import pytest
from aiida.orm import Dict
from aiida.plugins import CalculationFactory
from aiida import engine

# Calculation objects
calc_ch4_working_cap = CalculationFactory('lsmo.calc_ch4_working_cap')


def test_calc_ch4_working_cap():
    """Test extraction of working capacity data from IsothermWorkChain output dictionary."""
    OUTPUT_PARAMETERS = Dict(
        dict={
            'Density': 0.440527,
            'Density_unit': 'g/cm^3',
            'Estimated_saturation_loading': 41.985376,
            'Estimated_saturation_loading_unit': 'mol/kg',
            'Input_block': [1.865, 100],
            'Input_ha': 'DEF',
            'Input_structure_filename': 'tmp4a13iby3.cif',
            'Input_volpo': [1.865, 1.865, 100000],
            'Number_of_blocking_spheres': 0,
            'POAV_A^3': 8623.69,
            'POAV_A^3_unit': 'A^3',
            'POAV_Volume_fraction': 0.67999,
            'POAV_Volume_fraction_unit': None,
            'POAV_cm^3/g': 1.54358,
            'POAV_cm^3/g_unit': 'cm^3/g',
            'PONAV_A^3': 0.0,
            'PONAV_A^3_unit': 'A^3',
            'PONAV_Volume_fraction': 0.0,
            'PONAV_Volume_fraction_unit': None,
            'PONAV_cm^3/g': 0.0,
            'PONAV_cm^3/g_unit': 'cm^3/g',
            'Unitcell_volume': 12682.1,
            'Unitcell_volume_unit': 'A^3',
            'adsorption_energy_widom_average': -11.1626207486,
            'adsorption_energy_widom_dev': 0.02083606,
            'adsorption_energy_widom_unit': 'kJ/mol',
            'conversion_factor_molec_uc_to_cm3stp_cm3': 2.9347915768,
            'conversion_factor_molec_uc_to_mg_g': 4.7676018308,
            'conversion_factor_molec_uc_to_mol_kg': 0.2972320343,
            'henry_coefficient_average': 7.71003e-06,
            'henry_coefficient_dev': 1.65115e-08,
            'henry_coefficient_unit': 'mol/kg/Pa',
            'is_kh_enough': True,
            'is_porous': True,
            'isotherm': {
                'enthalpy_of_adsorption_average': [-13.510607783958, -10.787702310577],
                'enthalpy_of_adsorption_dev': [0.76886345231266, 1.0196832123586],
                'enthalpy_of_adsorption_unit': 'kJ/mol',
                'loading_absolute_average': [3.6279844874624, 16.11968088498],
                'loading_absolute_dev': [0.15865715470393, 0.075109385284932],
                'loading_absolute_unit': 'mol/kg',
                'pressure': [5.8, 65],
                'pressure_unit': 'bar'
            },
            'temperature': 298,
            'temperature_unit': 'K'
        })

    results = engine.run(calc_ch4_working_cap, OUTPUT_PARAMETERS)
    results_dict = results.get_dict()

    assert results_dict['wc_65bar_mol/kg_average'] == pytest.approx(12.5, abs=0.1)
