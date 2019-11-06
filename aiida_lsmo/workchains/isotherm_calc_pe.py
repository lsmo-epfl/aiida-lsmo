# -*- coding: utf-8 -*-
"""IsothermCalcPE work chain."""

from __future__ import absolute_import

from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Dict, Str
from aiida.engine import calcfunction
from aiida.engine import WorkChain, if_

IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  #pylint: disable=invalid-name
CifData = DataFactory('cif')  #pylint: disable=invalid-name


@calcfunction
def get_pe(geom_co2, isot_co2, geom_n2, isot_n2, pe_parameters):
    """ Submit calc_pe calculation using AiiDA, for the CO2 parasitic energy."""
    import pandas as pd
    from calc_pe import mainPE
    bar2pa = 1e5  # convert pressure from bar to Pa
    gcm2kgm = 1000  # convert density from g/cm3 to kg/m3

    t_iso = {}
    iso_df = {}
    for i in [0, 1]:
        gas = ['CO_2', 'N_2'][i]
        geom = [geom_co2, geom_n2][i].get_dict()
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
        rho=geom['Density'] * gcm2kgm,
        T_iso=t_iso,  # [T_iso_CO2, T_iso_N2]
        iso_df=iso_df,  # [iso_df_CO2, iso_df_N2] df containing pressure (Pa), loading (mol/kg), HoA (kJ/mol)
        **pe_parameters.get_dict()  # all other parameters that are specific of the modelling
    )
    return Dict(dict=pe_dict)


PE_PARAMETERS_DEFAULT = Dict(dict={
    'gasin': 'coal',
    'vf': 0.35,
    'process': 'TPSA',
    'cp': 985.0,
    'yd': 0.99,
    'eleff': 'carnot',
    'opt': 'PE',
})

ISOTHERM_PARAMETERS_DEFAULT = Dict(
    dict={
        "forcefield": "UFF",  # valid_type=Str, help='Forcefield of the structure.'
        "ff_tailcorr": True,  # apply tail corrections
        "ff_shift": False,  # shift or truncate at cutoff
        "ff_cutoff": 12.0,  # valid_type=Float, help='CutOff truncation for the VdW interactions (Angstrom)'
        "temperature": 300,  # valid_type=Float, help='Temperature of the simulation'
        "zeopp_volpo_samples": 1e5,  # valid_type=Int,help='Number of samples for VOLPO calculation (per UC volume)'
        "zeopp_block_samples": 100,  # valid_type=Int, help='Number of samples for BLOCK calculation (per A^3)'
        "raspa_minKh": 1e-10,
        "raspa_verbosity": 10,  # valid_type=Int,help='Print stats every: number of cycles / raspa_verbosity'
        "raspa_widom_cycles": 1e5,  # valid_type=Int, help='Number of widom cycles'
        "raspa_gcmc_init_cycles": 1e3,  # valid_type=Int, help='Number of GCMC initialization cycles'
        "raspa_gcmc_prod_cycles": 1e4,  # valid_type=Int, help='Number of GCMC production cycles'
        "pressure_precision": 0.1,
        "pressure_maxstep": 5,  # valid_type=Float, help='Max distance between pressure points (bar)'
        "pressure_min": 0.001,  # valid_type=Float, help='Lower pressure to sample (bar)'
        "pressure_max": 30
    })


class IsothermCalcPEWorkChain(WorkChain):
    """ Compute CO2 parassitic energy (PE) after running IsothermWorkChain for CO2 and N2 at 300K."""

    @classmethod
    def define(cls, spec):
        super(IsothermCalcPEWorkChain, cls).define(spec)

        spec.expose_inputs(IsothermWorkChain, exclude=['molecule', 'parameters'])
        spec.input('parameters',
                   valid_type=Dict,
                   default=ISOTHERM_PARAMETERS_DEFAULT,
                   help='Parameters for Isotherm work chain')

        spec.input('pe_parameters',
                   valid_type=Dict,
                   default=PE_PARAMETERS_DEFAULT,
                   help='Parameters for PE process modelling')

        spec.outline(
            cls.run_isotherms,  # run IsothermWorkChain in parallel for co2 and n2
            if_(cls.should_continue)(  # if porous to both co2 and n2
                cls.run_calcpe,),
        )

        spec.expose_outputs(IsothermWorkChain, namespace='co2')
        spec.expose_outputs(IsothermWorkChain, namespace='n2')
        spec.output('calc_pe', valid_type=Dict, required=False, help='Output of CalcPe')

    def run_isotherms(self):
        """Run Isotherm work chain for CO2 and N2."""

        inputs = self.exposed_inputs(IsothermWorkChain)  ###needs to be AttributeDict(se.fexposed/..)?????????
        self.report("Run Isotherm work chain for CO2 and N2, in CifData<{}> (label: {} )".format(
            self.inputs.structure.pk, self.inputs.structure.label))

        for mol in ['co2', 'n2']:

            inputs.update({
                'metadata': {
                    'call_link_label': 'run_isotherm_for_{}'.format(mol),
                },
                'molecule': Str(mol),  #Shoud I store????????????
                'parameters': self.inputs.parameters
            })

            running = self.submit(IsothermWorkChain, **inputs)
            self.to_context(**{'isotherm_{}'.format(mol): running})

    def should_continue(self):
        """Expose outputs and continue if porous to both."""

        self.out_many(self.exposed_outputs(self.ctx.isotherm_co2, IsothermWorkChain, namespace='co2'))
        self.out_many(self.exposed_outputs(self.ctx.isotherm_n2, IsothermWorkChain, namespace='n2'))

        co2_porous = self.ctx.isotherm_co2.outputs['geometric_output'].get_dict()['is_porous']
        n2_porous = self.ctx.isotherm_n2.outputs['geometric_output'].get_dict()['is_porous']

        self.report("Porosity test for CO2/N2: {}/{}".format(co2_porous, n2_porous))
        return co2_porous and n2_porous

    def run_calcpe(self):
        """Prepare isotherms results for calc_pe, run it and return the output."""
        from calc_pe.utils import get_PE_results_string

        calcpe = get_pe(geom_co2=self.ctx.isotherm_co2.outputs['geometric_output'],
                        isot_co2=self.ctx.isotherm_co2.outputs['isotherm_output'],
                        geom_n2=self.ctx.isotherm_n2.outputs['geometric_output'],
                        isot_n2=self.ctx.isotherm_n2.outputs['isotherm_output'],
                        pe_parameters=self.inputs['pe_parameters'])

        self.out("calc_pe", calcpe)
        self.report('calc_pe Dict<{}>'.format(calcpe.uuid))
        self.report("{}".format(get_PE_results_string(adsorbent_name=self.inputs.structure.label,
                                                      res=calcpe.get_dict())))
