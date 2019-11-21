=====================================
LSMO calc functions and work chains
=====================================

In the following section all the calc functions and work chains of the `aiida-lsmo` plugin are listed and documented.

Force Field Builder
+++++++++++++++++++++++

The :py:func:`~aiida_lsmo.calcfunctions.ff_builder_module.ff_builder` calculation function allows to combine the force
field parameters (typically for a Lennard-Jones potential) for a
framework and the molecule(s), giving as an output the `.def` files required by `Raspa`.
To see the list of available parameterization for the frameworks and the available molecules, give a look to the file
`ff_data.yaml`.

What it can do:

#. Switch settings that are written in the `.def` files of Raspa, such as tail-corrections, truncation/shifting and
   mixing rules.
#. Decide to separate the interactions, so that framework-molecule interactions and molecule-molecule interactions are
   parametrized differently (e.g., TraPPE for molecule-mololecule and UFF, instead of UFF/TraPPE for framework-molecule).

What it currently can not do:

#. Deal with flexible molecules.
#. Take parameters from other files (e.g., YAML).
#. Generate `.def` files for a molecule, given just the geometry: it has to be included in the ff_data.yaml file.


Inputs details
--------------

* Parameters Dict::

    PARAMS_EXAMPLE = Dict( dict = {
       'ff_framework': 'UFF',              # See force fields available in ff_data.yaml as framework.keys()
       'ff_molecules': {                   # See molecules available in ff_data.yaml as ff_data.keys(
           'CO2': 'TraPPE',                    # See force fields available in ff_data.yaml as {molecule}.keys()
           'N2': 'TraPPE'
       },
       'shifted': True,                    # If True shift despersion interactions, if False simply truncate them
       'tail_corrections': False,          # If True apply tail corrections based on homogeneous-liquid assumption
       'mixing_rule': 'Lorentz-Berthelot', # Options: 'Lorentz-Berthelot' or 'Jorgensen'
       'separate_interactions': True       # If True use framework's force field for framework-molecule interactions
    })


Outputs details
---------------

* Dictionary containing the `.def` files as SinglefileData. This output dictionary is ready to be used as a `files` input
  of the `RaspaCalculation`: you can find and example of usage of this CalcFunction in the `IsothermWorkChain`, or a
  minimal test usage in the examples.

Usage
------

The CalcFunction can be imported with::

  from aiida.plugins import CalculationFactory

  FFBuilder = CalculationFactory('lsmo.ff_builder')


Working Capacity calculators
+++++++++++++++++++++++++++++

The module `calcfunctions/working_cap.py` contains a collections of calculation functions to compute the working capacities
for different compound (e.g., CH4, H2) at industrially reference/relevant conditions.
The working capacity is the usable amount of a stored adsorbed compound between the loading and discharging
temperature and pressure.
These are post-processing calculation from the `output_parameters` of Isotherm or IsothermMultiTemp work chains,
that needs to be run at specific conditions: see the header of the calc function to know them.
Their inner working is very simple but they are collected in this repository to be used as a reference in our group.
If you are investigating some different gas storage application, consider including a similar script here.

An example is :py:func:`~aiida_lsmo.calcfunctions.working_cap.calc_ch4_working_cap` for methane storage.


Isotherm work chain
+++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.isotherm.IsothermWorkChain` work function allows to compute a single-component
isotherm in a framework, from a few settings.

What it does, in order:

#. Run a geometry calculation (Zeo++) to assess the accessible probe-occubiable pore volume and the needed blocking spheres.
#. Stop if the structure is non-porous, i.e., not permeable to the molecule.
#. Get the parameters of the force field using the FFBuilder.
#. Get the number of unit cell replicas needed to have correct periodic boundary conditions at the given cutoff.
#. Compute the adsorption at zero loading (e.g., the Henry coefficient, kH) from a Widom insertion calculation using Raspa.
#. Stop if the kH is not more that a certain user-defined threshold: this can be used for screening purpose, or to
   intentionally compute only the kH using this work chain.
#. Given a min/max range, propose a list of pressures that sample the isotherm uniformly. However, the user can also
   specify a defined list of pressure and skip this automatic selection.
#. Compute the isotherm using Grand Canonical Mon* te Carlo (GCMC) sampling in series, and restarting each system from
   the previous one for a short and efficient equilibration.

What it can not do:

#. Compute isotherms at different temperatures (see IsothermMultiTemp work chain for this).
#. Compute multi-component isotherms, as it would complicate a lot the input, output and logic, and it is not trivial
   to assign the mixture composition of the bulk gas at different pressure, for studying a real case.
#. It is not currently possible to play too much with Monte Carlo probabilities and other advanced settings in Raspa.
#. Sample the isotherm uniformly in case of "type II" isotherms, i.e., like for water, having significant cooperative insertion.
#. Run the different pressures in parallel: this would be less efficient because you can not restart from the previous
   configuration, and not necessarily much faster considering that equilibrating the higher pressure calculation will be
   anyway the bottleneck.

.. aiida-workchain:: IsothermWorkChain
    :module: aiida_lsmo.workchains

Inputs details
--------------

* ``structure`` (``CifData``) is the framework with partial charges (provided as ``_atom_site_charge`` column in the CIF file)

* ``molecule`` can be provided both as a ``Str`` or ``Dict``. It contains information about the molecule force field and
  approximated spherical-probe radius for the geometry calculation. If provided as a string (e.g., ``co2``, ``n2``)
  the work chain looks up at the corresponding dictionary in ``isotherm_data/isotherm_molecules.yaml``.
  The input dictionary reads as, for example::

      co2:
        name: CO2          # Raspa's MoleculeName
        forcefield: TraPPE # Raspa's MoleculeDefinition
        molsatdens: 21.2   # Density of the liquid phase of the molecule in (mol/l). Typically I run a simulation at 300K/200bar
        proberad: 1.525    # radius used for computing VOLPO and Block (Angs). Typically FF's sigma/2
        singlebead: False  # if true: RotationProbability=0
        charged: True      # if true: ChargeMethod=Ewald

* ``parameters`` (``Dict``) goes to modify these default parameters::

  parameters = {
    "ff_framework": "UFF",  # (str) Forcefield of the structure.
    "ff_separate_interactions": False,  # (bool) Use "separate_interactions" in the FF builder.
    "ff_mixing_rule": "Lorentz-Berthelot",  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
    "ff_tail_corrections": True,  # (bool) Apply tail corrections.
    "ff_shifted": False,  # (bool) Shift or truncate the potential at cutoff.
    "ff_cutoff": 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
    "temperature": 300,  # (float) Temperature of the simulation.
    "temperature_list": None,  # (list) To be used by IsothermMultiTempWorkChain.
    "zeopp_volpo_samples": int(1e5),  # (int) Number of samples for VOLPO calculation (per UC volume).
    "zeopp_block_samples": int(100),  # (int) Number of samples for BLOCK calculation (per A^3).
    "raspa_minKh": 1e-10,  # (float) If Henry coefficient < raspa_minKh do not run the isotherm (mol/kg/Pa).
    "raspa_verbosity": 10,  # (int) Print stats every: number of cycles / raspa_verbosity.
    "raspa_widom_cycles": int(1e5),  # (int) Number of Widom cycles.
    "raspa_gcmc_init_cycles": int(1e3),  # (int) Number of GCMC initialization cycles.
    "raspa_gcmc_prod_cycles": int(1e4),  # (int) Number of GCMC production cycles.
    "pressure_list": None,  # (list) Pressure list for the isotherm (bar): if given it will skip to guess it.
    "pressure_precision": 0.1,  # (float) Precision in the sampling of the isotherm: 0.1 ok, 0.05 for high resolution.
    "pressure_maxstep": 5,  # (float) Max distance between pressure points (bar).
    "pressure_min": 0.001,  # (float) Lower pressure to sample (bar).
    "pressure_max": 10  # (float) Upper pressure to sample (bar).
  }

Note that if the ``pressure_list`` value is provided, the other pressure inputs are neglected and the automatic pressure
selection of the work chain is skipped.

* ``geometric`` is not meant to be used by the user, but by the IsothermMultiTemp work chains.

Outputs details
---------------

* ``output_parameters`` (``Dict``) whose length depends whether ``is_porous`` is ``True`` (if not, only geometric outputs are
  reported in the dictionary), and whether ``is_kh_enough`` (if ``False``, it prints only the output of the Widom calculation,
  otherwise it also reports the isotherm data).
  This is an example of a full isotherm with ``is_porous=True`` and ``is_kh_enough=True``, for 6 pressure points at 298K ::

      {
          "Density": 0.385817,
          "Density_unit": "g/cm^3",
          "Estimated_saturation_loading": 51.586704,
          "Estimated_saturation_loading_unit": "mol/kg",
          "Input_block": [
              1.865,
              100
          ],
          "Input_ha": "DEF",
          "Input_structure_filename": "19366N2.cif",
          "Input_volpo": [
              1.865,
              1.865,
              100000
          ],
          "Number_of_blocking_spheres": 0,
          "POAV_A^3": 8626.94,
          "POAV_A^3_unit": "A^3",
          "POAV_Volume_fraction": 0.73173,
          "POAV_Volume_fraction_unit": null,
          "POAV_cm^3/g": 1.89657,
          "POAV_cm^3/g_unit": "cm^3/g",
          "PONAV_A^3": 0.0,
          "PONAV_A^3_unit": "A^3",
          "PONAV_Volume_fraction": 0.0,
          "PONAV_Volume_fraction_unit": null,
          "PONAV_cm^3/g": 0.0,
          "PONAV_cm^3/g_unit": "cm^3/g",
          "Unitcell_volume": 11789.8,
          "Unitcell_volume_unit": "A^3",
          "adsorption_energy_widom_average": -9.7886451805,
          "adsorption_energy_widom_dev": 0.0204010566,
          "adsorption_energy_widom_unit": "kJ/mol",
          "conversion_factor_molec_uc_to_cm3stp_cm3": 3.1569089445,
          "conversion_factor_molec_uc_to_gr_gr": 5.8556741651,
          "conversion_factor_molec_uc_to_mol_kg": 0.3650669679,
          "henry_coefficient_average": 6.72787e-06,
          "henry_coefficient_dev": 3.94078e-08,
          "henry_coefficient_unit": "mol/kg/Pa",
          "is_kh_enough": true,
          "is_porous": true,
          "isotherm": {
              "enthalpy_of_adsorption_average": [
                  -12.309803364014,
                  -12.058276670574,
                  -10.952120841752,
                  -10.367181994296,
                  -9.8224431917388,
                  -9.6064899852835
              ],
              "enthalpy_of_adsorption_dev": [
                  0.34443269062882,
                  0.25307818878718,
                  0.53612978806548,
                  0.73138412611309,
                  0.50295849442518,
                  0.2598580313121
              ],
              "enthalpy_of_adsorption_unit": "kJ/mol",
              "loading_absolute_average": [
                  0.65880897694654,
                  3.2677144296729,
                  8.5184817556528,
                  12.108148744317,
                  14.891264154125,
                  17.302504097082
              ],
              "loading_absolute_dev": [
                  0.041847687204507,
                  0.064694332028498,
                  0.18617312433529,
                  0.15188945177344,
                  0.14186942260037,
                  0.14638828764266
              ],
              "loading_absolute_unit": "mol/kg",
              "pressure": [
                  1.0,
                  5.8,
                  20,
                  35,
                  50,
                  65
              ],
              "pressure_unit": "bar"
          },
          "temperature": 298,
          "temperature_unit": "K"
      }

* ``block`` (``SinglefileData``) file is outputted if blocking spheres are found and used for the isotherm. Therefore,
  this is ready to be used for a new, consistent, Raspa calculation.

Usage
-----

Import as::

  from aiida.plugins import WorkflowFactory

  IsothermWorkChain = WorkflowFactory('lsmo.isotherm')

and see the examples for many usage application.

IsothermMultiTemp work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.isotherm_multi_temp.IsothermMultiTempWorkChain` work chain can run in parallel the
Isotherm work chain at different temperatures. Since the
geometry initial calculation to get the pore volume and blocking spheres is not dependent on the temperature, this is
run only once. Inputs and outputs are very similar to the Isotherm work chain.

What it can do:

#. Compute the kH at every temperature and guess, for each temperature, the pressure points needed for an uniform
   sampling of the isotherm.

What it can not do:

#. Select specific pressure points (as ``pressure_list``) that are different at different temperatures.
#. Run an isobar curve (same pressure, different pressures) restarting each GCMC calculation from the previous system.

.. aiida-workchain:: IsothermMultiTempWorkChain
    :module: aiida_lsmo.workchains

Inputs details
--------------

* ``parameters`` (``Dict``), compared to the input of the Isotherm work chain, contains the key ``temperature_list``
  and neglects the key ``temperature``::

    "temperature_list": [278, 298.15, 318.0],

Outputs details
---------------

* ``output_parameters`` (``Dict``) contains the ``temperature`` and ``isotherm`` as lists. In this example 3 pressure
  points are computed at 77K, 198K and 298K::

    {
        "Density": 0.731022,
        "Density_unit": "g/cm^3",
        "Estimated_saturation_loading": 22.1095656,
        "Estimated_saturation_loading_unit": "mol/kg",
        "Input_block": [
            1.48,
            100
        ],
        "Input_ha": "DEF",
        "Input_structure_filename": "tmpQD_OdI.cif",
        "Input_volpo": [
            1.48,
            1.48,
            100000
        ],
        "Number_of_blocking_spheres": 0,
        "POAV_A^3": 1579.69,
        "POAV_A^3_unit": "A^3",
        "POAV_Volume_fraction": 0.45657,
        "POAV_Volume_fraction_unit": null,
        "POAV_cm^3/g": 0.624564,
        "POAV_cm^3/g_unit": "cm^3/g",
        "PONAV_A^3": 0.0,
        "PONAV_A^3_unit": "A^3",
        "PONAV_Volume_fraction": 0.0,
        "PONAV_Volume_fraction_unit": null,
        "PONAV_cm^3/g": 0.0,
        "PONAV_cm^3/g_unit": "cm^3/g",
        "Unitcell_volume": 3459.91,
        "Unitcell_volume_unit": "A^3",
        "adsorption_energy_widom_average": [
            -6.501026119,
            -3.7417828535,
            -2.9538187687
        ],
        "adsorption_energy_widom_dev": [
            0.0131402719,
            0.0109470973,
            0.009493264
        ],
        "adsorption_energy_widom_unit": "kJ/mol",
        "conversion_factor_molec_uc_to_cm3stp_cm3": 10.757306634,
        "conversion_factor_molec_uc_to_gr_gr": 1.3130795208,
        "conversion_factor_molec_uc_to_mol_kg": 0.6565397604,
        "henry_coefficient_average": [
            0.000590302,
            1.36478e-06,
            4.59353e-07
        ],
        "henry_coefficient_dev": [
            6.20272e-06,
            2.92729e-09,
            1.3813e-09
        ],
        "henry_coefficient_unit": "mol/kg/Pa",
        "is_kh_enough": [
            true,
            true,
            true
        ],
        "is_porous": true,
        "isotherm": [
            {
                "enthalpy_of_adsorption_average": [
                    -4.8763191239929,
                    -4.071414615084,
                    -3.8884980003825
                ],
                "enthalpy_of_adsorption_dev": [
                    0.27048724983995,
                    0.17838206413742,
                    0.30520201541493
                ],
                "enthalpy_of_adsorption_unit": "kJ/mol",
                "loading_absolute_average": [
                    8.8763231830174,
                    13.809017193987,
                    24.592736102413
                ],
                "loading_absolute_dev": [
                    0.10377880404968,
                    0.057485479697981,
                    0.1444399097573
                ],
                "loading_absolute_unit": "mol/kg",
                "pressure": [
                    1.0,
                    5.0,
                    100
                ],
                "pressure_unit": "bar"
            },
            {
                "enthalpy_of_adsorption_average": [
                    -5.3762452088166,
                    -5.304498349588,
                    -5.1469837785704
                ],
                "enthalpy_of_adsorption_dev": [
                    0.16413676386221,
                    0.23624406142692,
                    0.16877234291986
                ],
                "enthalpy_of_adsorption_unit": "kJ/mol",
                "loading_absolute_average": [
                    0.13688033329639,
                    0.64822632568393,
                    8.2218063857542
                ],
                "loading_absolute_dev": [
                    0.0022470007645714,
                    0.015908634630445,
                    0.063314699465606
                ],
                "loading_absolute_unit": "mol/kg",
                "pressure": [
                    1.0,
                    5.0,
                    100
                ],
                "pressure_unit": "bar"
            },
            {
                "enthalpy_of_adsorption_average": [
                    -5.3995609987279,
                    -5.5404431584811,
                    -5.410077906097
                ],
                "enthalpy_of_adsorption_dev": [
                    0.095159861315507,
                    0.081469905963932,
                    0.1393537452296
                ],
                "enthalpy_of_adsorption_unit": "kJ/mol",
                "loading_absolute_average": [
                    0.04589212925196,
                    0.22723251444794,
                    3.8118903657499
                ],
                "loading_absolute_dev": [
                    0.0018452227888317,
                    0.0031557689853122,
                    0.047824194130595
                ],
                "loading_absolute_unit": "mol/kg",
                "pressure": [
                    1.0,
                    5.0,
                    100
                ],
                "pressure_unit": "bar"
            }
        ],
        "temperature": [
            77,
            198,
            298
        ],
        "temperature_unit": "K"
    }

IosthermCalcPE work chain
++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.isotherm_calc_pe.IsothermCalcPEWorkChain` work chain takes as an input a structure
with partial charges, computes the isotherms for CO2 and N2 at
ambient temperature and models the process of carbon capture and compression for geological sequestration.
The final outcome informs about the performance of the adsorbent for this application, including the CO2 parasitic energy,
i.e., the energy that is required to separate and compress one kilogram of CO2, using that material.
Default input mixture is coal post-combustion flue gas, but also natural gas post-combustion and air mixtures are available.

.. aiida-workchain:: IsothermCalcPEWorkChain
    :module: aiida_lsmo.workchains

MultistageDdec work chain
++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.multistage_ddec.MultistageDdecWorkChain` work chain combines together the CP2K
multistage workchain and the DDEC calculation, with the scope of
optimizing the geometry of a structure and compute its partial charge using the DDEC protocol.

.. aiida-workchain:: MultistageDdecWorkChain
    :module: aiida_lsmo.workchains

ZeoppMultistageDdec work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.zeopp_multistage_ddec.ZeoppMultistageDdecWorkChain` work chain, is similar to MultistageDdec
but it runs a geometry characterization of the structure
using Zeo++ (NetworkCalculation) before and after, with the scope of assessing the structural changes due to the cell/geometry
optimization.

.. aiida-workchain:: ZeoppMultistageDdecWorkChain
    :module: aiida_lsmo.workchains
