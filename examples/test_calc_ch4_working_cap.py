#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for ff_builder"""
from __future__ import absolute_import
from __future__ import print_function
import pprint
from aiida.orm import Dict
from aiida.plugins import CalculationFactory
from aiida.engine import run_get_node

# Calculation objects
calc_ch4_working_cap = CalculationFactory("lsmo.calc_ch4_working_cap")  # pylint: disable=invalid-name

# output_paramters from IsothermWorkChain
OUTPUT_PARAMETERS = Dict(
    dict={
        "Density": 0.440527,
        "Density_unit": "g/cm^3",
        "Estimated_saturation_loading": 41.985376,
        "Estimated_saturation_loading_unit": "mol/kg",
        "Input_block": [1.865, 100],
        "Input_ha": "DEF",
        "Input_structure_filename": "tmp4a13iby3.cif",
        "Input_volpo": [1.865, 1.865, 100000],
        "Number_of_blocking_spheres": 0,
        "POAV_A^3": 8623.69,
        "POAV_A^3_unit": "A^3",
        "POAV_Volume_fraction": 0.67999,
        "POAV_Volume_fraction_unit": None,
        "POAV_cm^3/g": 1.54358,
        "POAV_cm^3/g_unit": "cm^3/g",
        "PONAV_A^3": 0.0,
        "PONAV_A^3_unit": "A^3",
        "PONAV_Volume_fraction": 0.0,
        "PONAV_Volume_fraction_unit": None,
        "PONAV_cm^3/g": 0.0,
        "PONAV_cm^3/g_unit": "cm^3/g",
        "Unitcell_volume": 12682.1,
        "Unitcell_volume_unit": "A^3",
        "adsorption_energy_widom_average": -11.1626207486,
        "adsorption_energy_widom_dev": 0.02083606,
        "adsorption_energy_widom_unit": "kJ/mol",
        "conversion_factor_molec_uc_to_cm3stp_cm3": 2.9347915768,
        "conversion_factor_molec_uc_to_mg_g": 4.7676018308,
        "conversion_factor_molec_uc_to_mol_kg": 0.2972320343,
        "henry_coefficient_average": 7.71003e-06,
        "henry_coefficient_dev": 1.65115e-08,
        "henry_coefficient_unit": "mol/kg/Pa",
        "is_kh_enough": True,
        "is_porous": True,
        "isotherm": {
            "enthalpy_of_adsorption_average": [-13.510607783958, -10.787702310577],
            "enthalpy_of_adsorption_dev": [0.76886345231266, 1.0196832123586],
            "enthalpy_of_adsorption_unit": "kJ/mol",
            "loading_absolute_average": [3.6279844874624, 16.11968088498],
            "loading_absolute_dev": [0.15865715470393, 0.075109385284932],
            "loading_absolute_unit": "mol/kg",
            "pressure": [5.8, 65],
            "pressure_unit": "bar"
        },
        "temperature": 298,
        "temperature_unit": "K"
    })

results, node = run_get_node(calc_ch4_working_cap, OUTPUT_PARAMETERS)  # pylint: disable=invalid-name
print("Terminated calcfunction CH4 working capacity, pk:", node.pk)
print("Printing output Dict, pk:", results.pk)
pprint.pprint(results.get_dict())
