cp2k_code=cp2k@localhost
raspa_code=raspa@localhost
ddec_code=ddec@localhost
zeopp_code=network@localhost

# verdi run run_BindingSiteWorkChain_MOF74_CO2.py \
#     --cp2k_code $cp2k_code \
#     --raspa_code $raspa_code
# DONE

# verdi run run_BindingSiteWorkChain_MOF74_O2.py \ 
#     --cp2k_code $cp2k_code \
#     --raspa_code $raspa_code
# DONE

# verdi run run_Cp2kMultistageWorkChain_2_H2O.py \
#     --cp2k-code $cp2k_code \
# DONE

# verdi run run_Cp2kMultistageWorkChain_Cu_HKUST-1.py \
#     --cp2k-code $cp2k_code \
# TAKES VERY LONG...

# verdi run run_Cp2kPhonopyWorkChain.py \
#     --cp2k-code $cp2k_code \
#     --structure-pk 73
# NOT SURE HOW THIS ONE WORKS

# verdi run run_IsothermAccurateWorkChain_Graphite_Ar.py $raspa_code $zeopp_code 
# GOT STUCK

# verdi run run_IsothermCalcPE_HKUST-1.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_IsothermMultiTempWorkChain_HKUST-1.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_IsothermWorkChain_COF-1.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_IsothermWorkChain_HKUST-1_onlyKh.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_MulticompGcmcWorkChain_HKUST-1_2comp.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_SimAnnealingWorkChain_HKUST-1_3xCO2.py \
#     --raspa_code $raspa_code \
# DONE

# verdi run run_SinglecompWidomWorkChain_Box.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_SinglecompWidomWorkChain_HKUST-1.py \
#     --raspa_code $raspa_code \
#     --zeopp_code $zeopp_code
# DONE

# verdi run run_ZeoppMultistageDdecWorkChain_H2O.py \
#   --zeopp_code $zeopp_code \
#   --cp2k_code $cp2k_code \
#   --ddec_code $ddec_code \
#   /home/lsmo/project/lsmo/git/aiida-lsmo-codes/data/chargemol/atomic_densities
# TAKES VERY LONG...

# verdi run test_Cp2kBindingEnergy_CO2_MOF74.py \
#   --cp2k-code $cp2k_code \
# TAKES VERY LONG...

# verdi run test_Cp2kMultistageDdecWorkChain_H2O.py \
#   --ddec-code $ddec_code \
#   --cp2k-code $cp2k_code 
# TAKES VERY LONG

# verdi run test_Cp2kMultistageWorkChain_Al.py \
#   --cp2k-code $cp2k_code
# TAKES VERY LONG

# verdi run test_IsothermInflectionWorkChain_graphite_Ar.py \
#     --raspa-code $raspa_code \
#     --zeopp-code $zeopp_code
# DONE

# verdi run test_IsothermWorkChain_Mg-MOF74.py \
#     --raspa-code $raspa_code \
#     --zeopp-code $zeopp_code
# DONE

# verdi run test_MulticompAdsDesWorkChain_HKUST-1.py \
#     --raspa-code $raspa_code \
#     --zeopp-code $zeopp_code
# THIS ONE FAILED ONCE but then worked...

# verdi run test_MulticompGcmcWorkChain_Box_3comp.py \
#     --raspa-code $raspa_code \
#     --zeopp-code $zeopp_code
# DONE

verdi run test_SimAnnealingWorkChain_MOF74_CO2.py \
    --raspa-code $raspa_code \
# DONE

# Worrying warning:
# /home/lsmo/project/lsmo/git/oximachinerunner/oximachinerunner/__init__.py:319: UserWarning: Oximachine can only predict oxidation states of metals.                     This structure contains no metals.
#   warnings.warn(
