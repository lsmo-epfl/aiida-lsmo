"""IsothermMultiTemp workchain."""

from aiida.plugins import WorkflowFactory
from aiida.orm import Dict
from aiida.engine import calcfunction
from aiida.engine import WorkChain, ToContext, if_
from aiida_lsmo.utils import dict_merge

# import sub-workchains
IsothermWorkChain = WorkflowFactory('lsmo.isotherm')  # pylint: disable=invalid-name


def get_parameters_singletemp(i, parameters):
    parameters_singletemp = parameters.get_dict()
    parameters_singletemp['temperature'] = parameters_singletemp['temperature_list'][i]
    parameters_singletemp['temperature_list'] = None
    return Dict(dict=parameters_singletemp)


@calcfunction
def get_output_parameters(geom_dict, **isotherm_dict):
    """Gather together all the results, returning lists for the multi temperature values"""

    multi_temp_labels = [
        'temperature', 'henry_coefficient_average', 'henry_coefficient_dev', 'adsorption_energy_widom_average',
        'adsorption_energy_widom_dev', 'is_kh_enough', 'isotherm'
    ]

    single_temp_labels = [
        'temperature_unit', 'henry_coefficient_unit', 'adsorption_energy_widom_unit',
        'conversion_factor_molec_uc_to_cm3stp_cm3', 'conversion_factor_molec_uc_to_mg_g',
        'conversion_factor_molec_uc_to_mol_kg'
    ]

    out_dict = geom_dict.get_dict()
    if out_dict['is_porous']:
        for label in multi_temp_labels:
            out_dict[label] = []

        for i in range(len(isotherm_dict)):
            for label in multi_temp_labels:
                isotherm_out_i = isotherm_dict['isotherm_out_{}'.format(i)]
                out_dict[label].append(isotherm_out_i[label])

        # Same for all, take the last for convenience
        for label in single_temp_labels:
            out_dict[label] = isotherm_out_i[label]

    return Dict(dict=out_dict)


class IsothermMultiTempWorkChain(WorkChain):
    """ Run IsothermWorkChain for multiple temperatures: first compute geometric properties
    and then submit Widom+GCMC at different temperatures in parallel
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(IsothermWorkChain)

        spec.outline(
            cls.run_geometric,
            if_(cls.should_continue)(  # if porous
                cls.run_isotherms,  # run IsothermWorkChain in parallel at different temperatures
            ),
            cls.collect_isotherms)

        spec.expose_outputs(IsothermWorkChain, include=['block'])

        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='Results of isotherms run at different temperatures.')

    def run_geometric(self):
        """Perform Zeo++ block and VOLPO calculation with IsothermWC."""

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(IsothermWorkChain)

        # Set inputs for zeopp
        dict_merge(inputs, {
            'metadata': {
                'label': "IsothermGeometric",
                'call_link_label': 'run_geometric',
            },
        })

        running = self.submit(IsothermWorkChain, **inputs)
        self.report("Computing common gemetric properties")
        return ToContext(geom_only=running)

    def should_continue(self):
        """Continue if porous"""
        # Put the geometric_dict and in context for quick use
        self.ctx.geom = self.ctx.geom_only.outputs.output_parameters
        # Ouput block file
        if 'block' in self.ctx.geom_only.outputs:
            self.out_many(self.exposed_outputs(self.ctx.geom_only, IsothermWorkChain))
        return self.ctx.geom['is_porous']

    def run_isotherms(self):
        """Compute isotherms at different temperatures."""

        self.ctx.ntemp = len(self.inputs.parameters['temperature_list'])

        # create inputs: exposed are code and metadata
        inputs = self.exposed_inputs(IsothermWorkChain)
        inputs['geometric'] = self.ctx.geom
        if 'block' in self.ctx.geom_only.outputs:
            inputs['raspa_base']['raspa']["block_pocket"] = {"block_file": self.ctx.geom_only.outputs.block}

        # Update the parameters with only one temperature and submit
        for i in range(self.ctx.ntemp):
            self.ctx.parameters_singletemp = get_parameters_singletemp(i, self.inputs.parameters)

            dict_merge(
                inputs, {
                    'metadata': {
                        'label': "Isotherm_{}".format(i),
                        'call_link_label': 'run_isotherm_{}'.format(i),
                    },
                    'parameters': self.ctx.parameters_singletemp
                })

            running = self.submit(IsothermWorkChain, **inputs)
            self.to_context(**{'isotherm_{}'.format(i): running})

    def collect_isotherms(self):
        """ Collect all the results in one Dict """

        output_dict = {}
        if self.ctx.geom['is_porous']:
            for i in range(self.ctx.ntemp):
                output_dict['isotherm_out_{}'.format(i)] = self.ctx['isotherm_{}'.format(
                    i)].outputs['output_parameters']

        self.out("output_parameters", get_output_parameters(self.ctx.geom, **output_dict))

        self.report("All the isotherms computed: output Dict<{}>".format(self.outputs['output_parameters'].pk))
