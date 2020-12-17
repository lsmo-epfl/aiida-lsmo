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


**Inputs details**

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


**Outputs details**

* Dictionary containing the `.def` files as SinglefileData. This output dictionary is ready to be used as a `files` input
  of the `RaspaCalculation`: you can find and example of usage of this CalcFunction in the `IsothermWorkChain`, or a
  minimal test usage in the examples.

Selectivity calculators
+++++++++++++++++++++++

The :py:func:`~aiida_lsmo.calcfunctions.selectivity.calc_selectivity` calculation function
computes the selectivity of two gas in a material, as the ratio between their Henry coefficients.
In the future this module will host also different metrics to assess selectivity, for specific applications.

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
#. Compute the isotherm using Grand Canonical Monte Carlo (GCMC) sampling in series, and restarting each system from
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

**Inputs details**


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

* ``parameters`` (``Dict``) modifies the default parameters::

    parameters = {
      "ff_framework": "UFF",  # (str) Forcefield of the structure.
      "ff_separate_interactions": False,  # (bool) Use "separate_interactions" in the FF builder.
      "ff_mixing_rule": "Lorentz-Berthelot",  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
      "ff_tail_corrections": True,  # (bool) Apply tail corrections.
      "ff_shifted": False,  # (bool) Shift or truncate the potential at cutoff.
      "ff_cutoff": 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
      "temperature": 300,  # (float) Temperature of the simulation.
      "temperature_list": None,  # (list) To be used by IsothermMultiTempWorkChain.
      "zeopp_probe_scaling": 1.0, # (float), scaling probe's diameter: molecular_rad * scaling
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

Note that you can scale the probe radius to empirically account for some framework flexibility and avoid overblocking.
Setting ``zeopp_probe_scaling`` to zero (or a small value) basically corresponds to skipping the permeability check and
skips the calculation of blocking spheres.

* ``geometric`` is not meant to be used by the user, but by the IsothermMultiTemp work chains.

**Outputs details**

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
                  ...
                  -9.6064899852835
              ],
              "enthalpy_of_adsorption_dev": [
                  0.34443269062882,
                  ...
                  0.2598580313121
              ],
              "enthalpy_of_adsorption_unit": "kJ/mol",
              "loading_absolute_average": [
                  0.65880897694654,
                  ...
                  17.302504097082
              ],
              "loading_absolute_dev": [
                  0.041847687204507,
                  ...
                  0.14638828764266
              ],
              "loading_absolute_unit": "mol/kg",
              "pressure": [
                  1.0,
                  ...
                  65
              ],
              "pressure_unit": "bar"
          },
          "temperature": 298,
          "temperature_unit": "K"
      }

* ``block`` (``SinglefileData``) file is outputted if blocking spheres are found and used for the isotherm. Therefore,
  this is ready to be used for a new, consistent, Raspa calculation.

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

**Inputs details**

* ``parameters`` (``Dict``), compared to the input of the Isotherm work chain, contains the key ``temperature_list``
  and neglects the key ``temperature``::

    "temperature_list": [278, 298.15, 318.0],

**Outputs details**

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

IsothermCalcPE work chain
++++++++++++++++++++++++++

The :py:class:`~aiida_lsmo.workchains.isotherm_calc_pe.IsothermCalcPEWorkChain` work chain takes as an input a structure
with partial charges, computes the isotherms for CO2 and N2 at
ambient temperature and models the process of carbon capture and compression for geological sequestration.
The final outcome informs about the performance of the adsorbent for this application, including the CO2 parasitic energy,
i.e., the energy that is required to separate and compress one kilogram of CO2, using that material.
Default input mixture is coal post-combustion flue gas, but also natural gas post-combustion and air mixtures are available.

.. aiida-workchain:: IsothermCalcPEWorkChain
    :module: aiida_lsmo.workchains

CP2K multistage work chain
++++++++++++++++++++++++++++

The :py:class:`~aiida_lsmo.workchains.cp2k_multistage.Cp2kMultistageWorkChain` work chain in meant to automate
DFT optimizations in CP2K and guess some good parameters
for the simulation, but it is written in such a versatile fashion that it can be used for many other functions.

What it can do:

#. Given a protocol YAML with different settings, the work chains iterates until it converges the SCF calculation.
   The concept is to use general options for ``settings_0`` and more and more robust for the next ones.
#. The protocol YAML contains also a number of stages, i.e., different ``MOTION`` settings, that are executed one after
   the other, restarting from the previous calculation. During the first stage, ``stage_0``, different settings are
   tested until the SCF converges at the last step of ``stage_0``. If this dos not happening the work chain stops.
   Otherwise it continues running ``stage_1``, and all the other stages that are included in the protocol.
#. These stages can be used for running a robust cell optimization, i.e., combining first some MD steps to escape
   metastable geometries and later the final optimization, or ab-initio MD, first equilibrating the system with a shorter
   time constant for the thermostat, and then collecting statistics in the second stage.
#. Some default protocols are provided in ``workchains/multistage_protocols`` and they can be imported with simple tags
   such as ``test``, ``default``, ``robust_conv``. Otherwise, the user can take inspiration from these to write his
   own protocol and pass it to the work chain.
#. Compute the band gap.
#. You can restart from a previous calculation, e.g., from an already computed wavefunction.

What it can not do:

#. Run CP2K calculations with k-points.
#. Run CP2K advanced calculations, e.g., other than ``ENERGY``, ``GEO_OPT``, ``CELL_OPT`` and ``MD``.

.. aiida-workchain:: Cp2kMultistageWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

* ``structure`` (``StructureData``, NOTE this is not a ``CifData``) is the system to investigate. It can be also a molecule
  in a box and not necessarily a 2D/3D framework.

* ``protocol_tag`` (``Str``) calls a default protocol. Currently available:

+----------------+-----------------------------------------------------------------------------------------------------+
| ``default``    | Main choice, uses PBE-D3(BJ) with 600Ry/DZVP basis set and GTH pseudopotential.                     |
|                | First settings are with OT, and if not working it switches to diagonalization and smearing.         |
|                | As for the stages it runs a cell optimization, a short NPT MD and again cell optimization.          |
+----------------+-----------------------------------------------------------------------------------------------------+
| ``test``       | Quick protocol for testing purpose.                                                                 |
+----------------+-----------------------------------------------------------------------------------------------------+
|``robust_conv`` | Similar to ``default`` but using more robust and more expensive settings for the SCF convergence.   |
+----------------+-----------------------------------------------------------------------------------------------------+
|``singlepoint`` | Same settings as ``default`` but running only one stage for a single point calculation.             |
|                | Used to exploit the automation of this work chain for a simple energy calculation.                  |
+----------------+-----------------------------------------------------------------------------------------------------+

* ``protocol_yaml`` (``SinglefileData``) is used to specify a custom protocol through a YAML file. See the
  `default YAML file <https://github.com/aiidateam/aiida-cp2k/tree/master/aiida_cp2k/workchains/multistage_protocols/standard.yaml>`_
  as an example. Note that the dictionary need to contain the following keys:

+---------------------------+------------------------------------------------------------------------------------------+
| ``protocol_description``  | An user friendly description of the protocol.                                            |
+---------------------------+------------------------------------------------------------------------------------------+
| ``initial_magnetization`` | ``"element"`` for choice based on element, ``"oxidation_state"`` for choice based on     |
|                           | predicted oxidation state, or ``"zero"`` for no initial magnetization.                   |
|                           | To override default values, pass a dictionary element=>magnetization or a dictionary in  |
|                           | the form of the ``initial_magnetizations.yaml``.                                         |
+---------------------------+------------------------------------------------------------------------------------------+
| ``basis_set``             | Dictionary of ``KIND/BASIS_SET`` for each element.                                       |
+---------------------------+------------------------------------------------------------------------------------------+
| ``pseudopotential``       | Dictionary of ``KIND/POTENTIAL`` for each element.                                       |
+---------------------------+------------------------------------------------------------------------------------------+
| ``bandgap_thr_ev``        | Any ```stage_0`` using OT and evaluating a band gap below this threshold                 |
|                           | will be considered as a failure.                                                         |
+---------------------------+------------------------------------------------------------------------------------------+
| * ``settings_0``          | Settings updated in ``stage_0`` until the SCF converges.                                 |
| * ``settings_1``          |                                                                                          |
| * ...                     |                                                                                          |
+---------------------------+------------------------------------------------------------------------------------------+
| * ``stage_0``             | CP2K settings that are updated at every stage.                                           |
| * ``stage_1``             |                                                                                          |
| * ...                     |                                                                                          |
+---------------------------+------------------------------------------------------------------------------------------+

Other keys may be add in future to introduce new functionalities to the Multistage work chain.

* ``starting_settings_idx`` (``Int``) is used to start from a custom index of the settings. If for example you know that
  the material is conductive and needs for smearing, you can use ``Int(1)`` to update directly the settings to ``settings_1``
  that applies electron smearing: this is the case of ``default`` protocol.

* ``min_cell_size`` (``Float``) is used to extend the unit cell, so that the minimum perpendicular width of the cell is
  bigger than a certain specified value. This needed when a cell length is too narrow and the plane wave auxiliary basis
  set is not accurate enough at the Gamma point only. Also this may be needed for hybrid range-separated potentials that
  require a sufficient non-overlapping cutoff.

.. note:: Need to explain it further in Technicalities.

* ``parent_calc_folder`` (``RemoteData``) is used to restart from a previously computed wave function.

* ``cp2k_base.cp2k.parameters`` (``Dict``) can be used to specify some cp2k parameters that will be always overwritten
  just before submitting every calculation.

**Outputs details**

* ``output_structure`` (``StructureData``) is the final structure at the end of the last stage. It is not outputted in
  case of a single point calculation, since it does not update the geometry of the system.

* ``output_parameters`` (``Dict``), here it is an example for Aluminum, where the ``settings_0`` calculation is discarded
  because of a negative band gap, and therefore switched to ``settings_1`` which make the SCF converge and they are
  used for 2 stages::

    {
        "cell_resized": "1x1x1",
        "dft_type": "RKS",
        "final_bandgap_spin1_au": 6.1299999999931e-06,
        "final_bandgap_spin2_au": 6.1299999999931e-06,
        "last_tag": "stage_1_settings_1_valid",
        "natoms": 4,
        "nsettings_discarded": 1,
        "nstages_valid": 2,
        "stage_info": {
            "bandgap_spin1_au": [
                0.0,
                6.1299999999931e-06
            ],
            "bandgap_spin2_au": [
                0.0,
                6.1299999999931e-06
            ],
            "final_edens_rspace": [
                -3e-09,
                -3e-09
            ],
            "nsteps": [
                1,
                2
            ],
            "opt_converged": [
                true,
                false
            ]
        },
        "step_info": {
            "cell_a_angs": [
                4.05,
                4.05,
                4.05,
                4.05
            ],
            "cell_alp_deg": [
                90.0,
                90.0,
                90.0,
                90.0
            ],
            "cell_b_angs": [
                4.05,
                4.05,
                4.05,
                4.05
            ],
            "cell_bet_deg": [
                90.0,
                90.0,
                90.0,
                90.0
            ],
            "cell_c_angs": [
                4.05,
                4.05,
                4.05,
                4.05
            ],
            "cell_gam_deg": [
                90.0,
                90.0,
                90.0,
                90.0
            ],
            "cell_vol_angs3": [
                66.409,
                66.409,
                66.409,
                66.409
            ],
            "dispersion_energy_au": [
                -0.04894693184602,
                -0.04894693184602,
                -0.04894696543385,
                -0.04894705992872
            ],
            "energy_au": [
                -8.0811276714482,
                -8.0811276714483,
                -8.0811249649336,
                -8.0811173120933
            ],
            "max_grad_au": [
                null,
                0.0,
                null,
                null
            ],
            "max_step_au": [
                null,
                0.0,
                null,
                null
            ],
            "pressure_bar": [
                null,
                null,
                58260.2982324,
                58201.2710544
            ],
            "rms_grad_au": [
                null,
                0.0,
                null,
                null
            ],
            "rms_step_au": [
                null,
                0.0,
                null,
                null
            ],
            "scf_converged": [
                true,
                true,
                true,
                true
            ],
            "step": [
                0,
                1,
                1,
                2
            ]
        }
    }

* ``last_input_parameters`` (``Dict``) reports the inputs that were used for the last CP2K calculation. They are possibly
  the ones that make the SCF converge, so the user can inspect them and use them for other direct CP2K calculations in AiiDA.


**Usage**

See examples provided with the `plugin <https://github.com/aiidateam/aiida-cp2k/tree/master/examples/workchains>`_.
The report provides very useful insight on what happened during the run. Here it is the example of Aluminum::

  2019-11-22 16:54:52 [90962 | REPORT]: [266248|Cp2kMultistageWorkChain|setup_multistage]: Unit cell was NOT resized
  2019-11-22 16:54:52 [90963 | REPORT]: [266248|Cp2kMultistageWorkChain|run_stage]: submitted Cp2kBaseWorkChain for stage_0/settings_0
  2019-11-22 16:54:52 [90964 | REPORT]:   [266252|Cp2kBaseWorkChain|run_calculation]: launching Cp2kCalculation<266253> iteration #1
  2019-11-22 16:55:13 [90965 | REPORT]:   [266252|Cp2kBaseWorkChain|inspect_calculation]: Cp2kCalculation<266253> completed successfully
  2019-11-22 16:55:13 [90966 | REPORT]:   [266252|Cp2kBaseWorkChain|results]: work chain completed after 1 iterations
  2019-11-22 16:55:14 [90967 | REPORT]:   [266252|Cp2kBaseWorkChain|on_terminated]: remote folders will not be cleaned
  2019-11-22 16:55:14 [90968 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_settings_stage0]: Bandgaps spin1/spin2: -0.058 and -0.058 ev
  2019-11-22 16:55:14 [90969 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_settings_stage0]: BAD SETTINGS: band gap is < 0.100 eV
  2019-11-22 16:55:14 [90970 | REPORT]: [266248|Cp2kMultistageWorkChain|run_stage]: submitted Cp2kBaseWorkChain for stage_0/settings_1
  2019-11-22 16:55:15 [90971 | REPORT]:   [266259|Cp2kBaseWorkChain|run_calculation]: launching Cp2kCalculation<266260> iteration #1
  2019-11-22 16:55:34 [90972 | REPORT]:   [266259|Cp2kBaseWorkChain|inspect_calculation]: Cp2kCalculation<266260> completed successfully
  2019-11-22 16:55:34 [90973 | REPORT]:   [266259|Cp2kBaseWorkChain|results]: work chain completed after 1 iterations
  2019-11-22 16:55:34 [90974 | REPORT]:   [266259|Cp2kBaseWorkChain|on_terminated]: remote folders will not be cleaned
  2019-11-22 16:55:35 [90975 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_settings_stage0]: Bandgaps spin1/spin2: 0.000 and 0.000 ev
  2019-11-22 16:55:35 [90976 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_stage]: Structure updated for next stage
  2019-11-22 16:55:35 [90977 | REPORT]: [266248|Cp2kMultistageWorkChain|run_stage]: submitted Cp2kBaseWorkChain for stage_1/settings_1
  2019-11-22 16:55:35 [90978 | REPORT]:   [266266|Cp2kBaseWorkChain|run_calculation]: launching Cp2kCalculation<266267> iteration #1
  2019-11-22 16:55:53 [90979 | REPORT]:   [266266|Cp2kBaseWorkChain|inspect_calculation]: Cp2kCalculation<266267> completed successfully
  2019-11-22 16:55:53 [90980 | REPORT]:   [266266|Cp2kBaseWorkChain|results]: work chain completed after 1 iterations
  2019-11-22 16:55:54 [90981 | REPORT]:   [266266|Cp2kBaseWorkChain|on_terminated]: remote folders will not be cleaned
  2019-11-22 16:55:54 [90982 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_stage]: Structure updated for next stage
  2019-11-22 16:55:54 [90983 | REPORT]: [266248|Cp2kMultistageWorkChain|inspect_and_update_stage]: All stages computed, finishing...
  2019-11-22 16:55:55 [90984 | REPORT]: [266248|Cp2kMultistageWorkChain|results]: Outputs: Dict<266273> and StructureData<266271>

Cp2kMultistageDdec work chain
++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.cp2k_multistage_ddec.Cp2kMultistageDdecWorkChain` work chain combines together the CP2K
Multistage workchain and the DDEC calculation, with the scope of
optimizing the geometry of a structure and compute its partial charge using the DDEC protocol.

.. aiida-workchain:: Cp2kMultistageDdecWorkChain
    :module: aiida_lsmo.workchains

ZeoppMultistageDdec work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.zeopp_multistage_ddec.ZeoppMultistageDdecWorkChain` work chain, is similar to Cp2kMultistageDdec
but it runs a geometry characterization of the structure
using Zeo++ (NetworkCalculation) before and after, with the scope of assessing the structural changes due to the cell/geometry
optimization.

.. aiida-workchain:: ZeoppMultistageDdecWorkChain
    :module: aiida_lsmo.workchains

SimAnnealing work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.sim_annealing.SimAnnealingWorkChain` work chain
is designed to find the minimum energy configuration for a given number of gas molecules in the pores of a framework.
It runs several NVT Monte-Carlo simulations in RASPA at decreasing temperature in order to let the system move to its global minimum (simulated annealing),
and then performs a geometry optimization for the final fine tuning of the adsorbate position(s).

.. aiida-workchain:: SimAnnealingWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

* ``parameters`` (``Dict``) modifies the default parameters::

    PARAMETERS_DEFAULT = {
        "ff_framework": "UFF",  # (str) Forcefield of the structure.
        "ff_separate_interactions": False,  # (bool) Use "separate_interactions" in the FF builder.
        "ff_mixing_rule": "Lorentz-Berthelot",  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
        "ff_tail_corrections": True,  # (bool) Apply tail corrections.
        "ff_shifted": False,  # (bool) Shift or truncate the potential at cutoff.
        "ff_cutoff": 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
        "temperature_list": [300, 250, 200, 250, 100, 50],  # (list) List of decreasing temperatures for the annealing.
        "mc_steps": int(1e3),  # (int) Number of MC cycles.
        "number_of_molecules": 1  # (int) Number of molecules loaded in the framework.
    }


**Outputs details**

* ``output_parameters`` (``Dict``), example::

    {
        "description": [
            "NVT simulation at 300 K",
            "NVT simulation at 250 K",
            "NVT simulation at 200 K",
            "NVT simulation at 250 K",
            "NVT simulation at 100 K",
            "NVT simulation at 50 K",
            "Final energy minimization"
        ],
        "energy_ads/ads_coulomb_final": [
            -0.00095657162276787,
            ...
            3.5423777787399e-06
        ],
        "energy_ads/ads_tot_final": [
            -0.00095657162276787,
            ...
            3.5423777787399e-06
        ],
        "energy_ads/ads_vdw_final": [
            0.0,
            ...
            0.0
        ],
        "energy_host/ads_coulomb_final": [
            -12.696035310164,
            ...
            -15.592788991158
        ],
        "energy_host/ads_tot_final": [
            -30.545798720022,
            ...
            -36.132005060753
        ],
        "energy_host/ads_vdw_final": [
            -17.849763409859,
            ...
            -20.539216069678
        ],
        "energy_unit": "kJ/mol",
        "number_of_molecules": 1
    }


  In particular:

  * ``host/ads`` describes the interaction between host and (all) adsorbates, ``ads/ads`` the interaction between adsorbates
  * The binding energy is the final value of ``energy_host/ads_coulomb_final`` (you still need to divide by ``number_of_molecules``)
  * ``energy_ads/ads_coulomb_final`` and ``energy_ads/ads_tot_final`` may be non-zero even for a single molecule due to rounding errors in Ewald summation

Cp2kBindingEnergy work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.cp2k_binding_energy.Cp2kBindingEnergyWorkChain` work chain
takes as an input a CIF structure and the initial position of a molecule in its pore,
optimizes the molecule's geometry keeping the framework rigid and computes the BSSE corrected interactions energy.
The work chain is similar to CP2K's MulstistageWorkChain in reading the settings from YAML protocol,
and resubmitting the calculation with updated settings in case of failure,
but the only step is an hard-coded ``GEO_OPT`` simulation with 200 max steps.

NOTE:

#. It is better to start with the settings of a previous working MulstistageWorkChain, if already available.
   Otherwise, it may run for 200 steps before realizing that the settings are not good an switch them.
#. No restart is allowed, since the system is changing the number of atoms for the BSSE calculation: therefore, the
   wave function is recomputed 5 times from scratch. This needs to be fixed in the future.
#. If ``structure`` and ``molecule`` ``StructureData`` do not have the same size for the unit cell,
   the work chain will complain and stop.

.. aiida-workchain:: Cp2kBindingEnergyWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

Look at the inputs details of the Multistage work chain for more information about the choice of the protocol
(i.e., DFT settings).

**Outputs details**

* ``output_parameters`` (``Dict``), example::

    {
        "binding_energy_bsse": -1.7922110202537,
        "binding_energy_corr": -23.072114381515,
        "binding_energy_dispersion": -18.318476834858,
        "binding_energy_raw": -24.864325401768,
        "binding_energy_unit": "kJ/mol",
        "motion_opt_converged": false,
        "motion_step_info": {
            "dispersion_energy_au": [
                -0.1611999344803,
                ...
                -0.16105256797101
            ],
            "energy_au": [
                -829.9150365907,
                ...
                -829.91870835924
            ],
            "max_grad_au": [
                null,
                0.0082746554,
                ...
                0.0030823925
            ],
            "max_step_au": [
                null,
                0.0604411557,
                ...
                0.0215865148
            ],
            "rms_grad_au": [
                null,
                0.000915767,
                ...
                0.0003886735
            ],
            "rms_step_au": [
                null,
                0.0071240711,
                ...
                0.0026174255
            ],
            "scf_converged": [
                true,
                ...
                true
            ]
        }
    }

BindingSite work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.binding_site.BindingSiteWorkChain` work chain
simply combines :py:func:`~aiida_lsmo.workchains.sim_annealing.SimAnnealingWorkChain`
and :py:func:`~aiida_lsmo.workchains.cp2k_binding_energy.Cp2kBindingEnergyWorkChain`.
The outputs from the two workchain are collected under the ``ff`` and ``dft`` namespaces, respectively.

.. aiida-workchain:: BindingSiteWorkChain
    :module: aiida_lsmo.workchains


SinglecompWidom work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.singlecomp_widom.SinglecompWidomWorkChain` work chain
allows to compute the Henry's coefficient of a molecule via the Widom insertions algorithm.
The user can specify a list of temperatures to perform these calculations, and the results from the ``output_parameters``
Dict will be presented as lists as well, one for each temperature.

Blocking spheres are computed for the molecule before the calculation.

.. aiida-workchain:: SinglecompWidomWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**


* ``structure`` (``CifData``), if missing the calculation will be performed for an empty box, which is convenient to get the
  ``widom_rosenbluth_factor_average`` for flexible molecules.

* ``molecule`` (see IsothermWorkChain)

* ``parameters`` (``Dict``) modifies the default parameters::

        "ff_framework": "UFF",  # str, Forcefield of the structure (used also as a definition of ff.rad for zeopp)
        "ff_shifted": False,  # bool, Shift or truncate at cutoff
        "ff_tail_corrections": True,  # bool, Apply tail corrections
        "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
        "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "zeopp_probe_scaling": 1.0,  # float, scaling probe's diameter: use 0.0 for skipping block calc
        "zeopp_block_samples": int(1000),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_widom_cycles": int(1e5),  # int, Number of widom cycles
        "temperatures": [300, 400]

**Outputs details**

* ``output_parameters`` (``Dict``), example::

    {
        "adsorption_energy_widom_average": [
            -34.9698999639,
            -34.8262538296,
            -34.6772296828
        ],
        "adsorption_energy_widom_dev": [
            0.0166320673,
            0.0129078639,
            0.015868202
        ],
        "adsorption_energy_widom_unit": "kJ/mol",
        "henry_coefficient_average": [
            0.00783847,
            0.0025542,
            0.000964045
        ],
        "henry_coefficient_dev": [
            0.000100367,
            1.78042e-05,
            4.69145e-06
        ],
        "henry_coefficient_unit": "mol/kg/Pa",
        "temperatures": [
            273,
            293,
            313
        ],
        "temperatures_unit": "K",
        "widom_rosenbluth_factor_average": [
            21180.0,
            7407.21,
            2986.58
        ],
        "widom_rosenbluth_factor_dev": [
            271.198648,
            51.63243,
            14.533953
        ],
        "widom_rosenbluth_factor_unit": "-"
    }

MulticompGcmc work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.multicomp_gcmc.MulticompGcmcWorkChain` work chain
performs in parallel GCMC calcultions, at all the conditions of temperature and pressure specified,
for a given mixture of molecules.

Blocking spheres are computed for each molecule before the calculation.

.. aiida-workchain:: MulticompGcmcWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

* ``parameters`` (``Dict``) modifies the default parameters::

        "ff_framework": "UFF",  # str, Forcefield of the structure (used also as a definition of ff.rad for zeopp)
        "ff_shifted": False,  # bool, Shift or truncate at cutoff
        "ff_tail_corrections": True,  # bool, Apply tail corrections
        "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
        "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "zeopp_probe_scaling": 1.0,  # float, scaling probe's diameter: use 0.0 for skipping block calc
        "zeopp_block_samples": int(1000),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_gcmc_init_cycles": int(1e5),  # int, Number of GCMC initialization cycles
        "raspa_gcmc_prod_cycles": int(1e5),  # int, Number of GCMC production cycles

* ``conditions`` (``Dict``), example::

        'molfraction': {
            'co': 0.2,
            'ethene': 0.3,
            'ethane': 0.5,
        },
        'temp_press': [
            [200, 0.1],
            [300, 0.5],
            [400, 0.7],
        ]

**Outputs details**

* ``output_parameters`` (``Dict``), example::

    "Input_block": {
        "C2H4": [
            1.647,
            10
        ],
        "C2H6": [
            1.683,
            10
        ],
        "CO": [
            1.584,
            10
        ]
    },
    "Number_of_blocking_spheres": {
        "C2H4": 0,
        "C2H6": 0,
        "CO": 0
    },
    "composition": {
        "C2H4": 0.3,
        "C2H6": 0.5,
        "CO": 0.2
    },
    "enthalpy_of_adsorption_average": {
        "C2H4": -18.893613196644,
        "C2H6": -23.953846166638,
        "CO": -18.67727295403
    },
    "enthalpy_of_adsorption_dev": {
        "C2H4": 8.3425044773141,
        "C2H6": 7.6573330506431,
        "CO": 12.788154764577
    },
    "enthalpy_of_adsorption_unit": "kJ/mol",
    "loading_absolute_average": {
        "C2H4": [
            1.674941508006,
            0.4649745087,
            0.225254317548
        ],
        "C2H6": [
            8.558630790138,
            1.716272575446,
            0.634431885204
        ],
        "CO": [
            0.153958226214,
            0.044430897498,
            0.018598980348
        ]
    },
    "loading_absolute_dev": {
        "C2H4": [
            0.37769902489165,
            0.14568834858221,
            0.060107934935359
        ],
        "C2H6": [
            0.20733617827358,
            0.1483606861503,
            0.1634276608098
        ],
        "CO": [
            0.098979061794361,
            0.033704465508572,
            0.014016046900518
        ]
    },
    "loading_absolute_unit": "mol/kg",
    "pressures": [
        0.1,
        0.5,
        1.0
    ],
    "pressures_unit": "bar",
    "temperatures": [
        200,
        300,
        400
    ],
    "temperatures_unit": "K"

MulticompAdsDes work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.multicomp_ads_des.MulticompAdsDesWorkChain` work chain
is similar to MulticompGcmc, but it performs one simulation at given adsorption temperature, pressure and composition,
and a second one at given temperature and pressure for desorption. For the desorption mixure of the gas reservoir,
the workchains uses the composition previously obtained at adsorption conditions inside the framework.

Note that this is an approximation - in order to arrive at the appropriate mixture for the gas reservoir at desorption, one should iterate, taking as the next desorption condition trial the difference between the mixture inside the framework at adsorption and the mixture inside the framework at desorption.
The approximation may induce artifacts such as negative working capacity for certain components, which are in any case a warning sign that that the desorption (partial) pressure is not low enough to evacuate the component from the framework.

.. aiida-workchain:: MulticompAdsDesWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

* ``parameters`` (``Dict``) modifies the default parameters::

        "ff_framework": "UFF",  # str, Forcefield of the structure (used also as a definition of ff.rad for zeopp)
        "ff_shifted": False,  # bool, Shift or truncate at cutoff
        "ff_tail_corrections": True,  # bool, Apply tail corrections
        "ff_mixing_rule": 'Lorentz-Berthelot',  # str, Mixing rule for the forcefield
        "ff_separate_interactions": False,  # bool, if true use only ff_framework for framework-molecule interactions
        "ff_cutoff": 12.0,  # float, CutOff truncation for the VdW interactions (Angstrom)
        "zeopp_probe_scaling": 1.0,  # float, scaling probe's diameter: use 0.0 for skipping block calc
        "zeopp_block_samples": int(1000),  # int, Number of samples for BLOCK calculation (per A^3)
        "raspa_verbosity": 10,  # int, Print stats every: number of cycles / raspa_verbosity
        "raspa_gcmc_init_cycles": int(1e5),  # int, Number of GCMC initialization cycles
        "raspa_gcmc_prod_cycles": int(1e5),  # int, Number of GCMC production cycles

* ``conditions`` (``Dict``), example::

            'molfraction': {
                'xenon': 0.2,
                'krypton': 0.8
            },
            'adsorption': {
                'temperature': 298, #K
                'pressure': 1, #bar
            },
            'desorption': {
                'temperature': 308,
                'pressure': 0.1,
            },

**Outputs details**

* ``output_parameters`` (``Dict``), example::

    "Input_block": {
        "Kr": [
            1.647,
            10
        ],
        "Xe": [
            1.7865,
            10
        ]
    },
    "Number_of_blocking_spheres": {
        "Kr": 0,
        "Xe": 0
    },
    "composition": {
        "Kr": [
            0.8,
            0.37188099808061
        ],
        "Xe": [
            0.2,
            0.62811900191939
        ]
    },
    "loading_absolute_average": {
        "Kr": [
            0.80078943165,
            0.042364344126
        ],
        "Xe": [
            1.352559181974,
            0.684029166132
        ]
    },
    "loading_absolute_dev": {
        "Kr": [
            0.18747637335777,
            0.02392208478975
        ],
        "Xe": [
            0.20357386402562,
            0.20233235593516
        ]
    },
    "loading_absolute_unit": "mol/kg",
    "pressures": [
        1,
        0.1
    ],
    "pressures_unit": "bar",
    "temperatures": [
        298,
        308
    ],
    "temperatures_unit": "K",
    "working_capacity": {
        "Kr": 0.758425087524,
        "Xe": 0.668530015842
    },
    "working_capacity_unit": "mol/kg"

IsothermInflection work chain
++++++++++++++++++++++++++++++++++++++++++++++

The :py:func:`~aiida_lsmo.workchains.isotherm_inflection.IsothermInflectionWorkChain` work chain
is designed to compute those isotherms that may have hysteresis between adsorption and desorption.
The work chain computes in parallel the uptake via GCMC at all pressure points from both the bare and the saturated
framework. The saturated framework is obtaining by running a "quasi-NVT" simulation, initializated with a number
of molecules equal to ``90% * pore volume * fluid density``. "Quasi-NVT" is defined as a GCMC calculations where the
swap move is only rarely attempted.

Note that this work chain may run many calculations in parallel.

.. aiida-workchain:: IsothermInflectionWorkChain
    :module: aiida_lsmo.workchains

**Inputs details**

* ``parameters`` (``Dict``) modifies the default parameters::

        "ff_framework": "UFF",  # (str) Forcefield of the structure.
        "ff_separate_interactions": False,  # (bool) Use "separate_interactions" in the FF builder.
        "ff_mixing_rule": "Lorentz-Berthelot",  # (string) Choose 'Lorentz-Berthelot' or 'Jorgensen'.
        "ff_tail_corrections": True,  # (bool) Apply tail corrections.
        "ff_shifted": False,  # (bool) Shift or truncate the potential at cutoff.
        "ff_cutoff": 12.0,  # (float) CutOff truncation for the VdW interactions (Angstrom).
        "temperature": 300,  # (float) Temperature of the simulation.
        "zeopp_probe_scaling": 1.0,  # float, scaling probe's diameter: use 0.0 for skipping block calc
        "zeopp_volpo_samples": int(1e5),  # (int) Number of samples for VOLPO calculation (per UC volume).
        "zeopp_block_samples": int(100),  # (int) Number of samples for BLOCK calculation (per A^3).
        "raspa_verbosity": 10,  # (int) Print stats every: number of cycles / raspa_verbosity.
        "raspa_widom_cycles": int(1e5),  # (int) Number of Widom cycles.
        "raspa_gcmc_init_cycles": int(1e3),  # (int) Number of GCMC initialization cycles.
        "raspa_gcmc_prod_cycles": int(1e4),  # (int) Number of GCMC production cycles.
        "pressure_min": 0.001,  # (float) Min pressure in P/P0 TODO: MIN selected from the henry coefficient!
        "pressure_max": 1.0,  # (float) Max pressure in P/P0
        "pressure_num": 20,  # (int) Number of pressure points considered, eqispaced in a log plot
        "pressure_list": None,  # (list) Pressure list in P/P0. If 'None' pressure points are computed from min/max/num.

* ``molecule`` (``Dict``), example::

            'name': 'Ar',
            'forcefield': 'HIRSCHFELDER',
            "ff_cutoff": 8,
            'molsatdens': 35.4, # NOTE: very important to define the initial amount of molecules!
            'proberad': 1.7,
            'singlebead': True,
            'charged': False,
            'pressure_zero': 1, # Saturation pressure @ T (bar)

**Outputs details**

* ``output_parameters`` (``Dict``), example::

    "Density": 0.380639,
    "Density_unit": "g/cm^3",
    "Estimated_saturation_loading": 77.944428,
    "Estimated_saturation_loading_unit": "mol/kg",
    "Input_block": [
        1.7,
        100
    ],
    "Input_ha": "DEF",
    "Input_structure_filename": "Graphite_20A.cif",
    "Input_volpo": [
        1.7,
        1.7,
        10000
    ],
    "Number_of_blocking_spheres": 0,
    "POAV_A^3": 175.659,
    "POAV_A^3_unit": "A^3",
    "POAV_Volume_fraction": 0.8381,
    "POAV_Volume_fraction_unit": null,
    "POAV_cm^3/g": 2.20182,
    "POAV_cm^3/g_unit": "cm^3/g",
    "PONAV_A^3": 0.0,
    "PONAV_A^3_unit": "A^3",
    "PONAV_Volume_fraction": 0.0,
    "PONAV_Volume_fraction_unit": null,
    "PONAV_cm^3/g": 0.0,
    "PONAV_cm^3/g_unit": "cm^3/g",
    "Unitcell_volume": 209.592,
    "Unitcell_volume_unit": "A^3",
    "adsorption_energy_widom_average": -10.349783334,
    "adsorption_energy_widom_dev": 0.0203871821,
    "adsorption_energy_widom_unit": "kJ/mol",
    "henry_coefficient_average": 0.387019,
    "henry_coefficient_dev": 0.0244542,
    "henry_coefficient_unit": "mol/kg/Pa",
    "is_porous": true,
    "isotherm": {
        "conversion_factor_molec_uc_to_cm3stp_cm3": 177.5796535584,
        "conversion_factor_molec_uc_to_mg_g": 831.5477157974,
        "conversion_factor_molec_uc_to_mol_kg": 20.814711284,
        "enthalpy_of_adsorption_average_from_dil": [
            -10.813400929552,
            -7.4639574135508,
            -9.9392383993082,
            null
        ],
        "enthalpy_of_adsorption_average_from_sat": [
            null,
            null,
            -14.825644658402,
            null
        ],
        "enthalpy_of_adsorption_dev_from_dil": [
            4.8201465665611,
            3.0478392822994,
            5.5941478469815,
            null
        ],
        "enthalpy_of_adsorption_dev_from_sat": [
            null,
            null,
            7.2192916198981,
            null
        ],
        "enthalpy_of_adsorption_unit": "kJ/mol",
        "loading_absolute_average_from_dil": [
            27.31930856025,
            29.798489351576,
            62.091027143015,
            72.617323992055
        ],
        "loading_absolute_average_from_sat": [
            29.151746536085,
            30.534438070785,
            72.907243184345,
            77.211428125155
        ],
        "loading_absolute_dev_from_dil": [
            0.94383239273726,
            1.8944736272227,
            18.62478172839,
            2.3400142831813
        ],
        "loading_absolute_dev_from_sat": [
            0.32025900270015,
            0.26491967913899,
            1.1170544140921,
            0.40142657506021
        ],
        "loading_absolute_unit": "mol/kg",
        "pressure": [
            0.001,
            0.01,
            0.1,
            1.0
        ],
        "pressure_unit": "bar"
        "temperature": 87,
        "temperature_unit": "K"
    }
