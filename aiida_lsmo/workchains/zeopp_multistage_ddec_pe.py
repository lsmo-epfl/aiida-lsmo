# -*- coding: utf-8 -*-
"""ZeoppMultistageDdecPeWorkChain workchain"""

from __future__ import absolute_import

from aiida.orm import Group
from aiida.plugins import WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext

# import sub-workchains
ZeoppMultistageDdecWorkChain = WorkflowFactory('lsmo.zeoppmultistageddec')  # pylint: disable=invalid-name
IsothermCalcPEWorkChain = WorkflowFactory('lsmo.isotherm_calc_pe')  # pylint: disable=invalid-name


def include_node(tag, node, group):
    """Given an aiida-node and a (string) tag, add the node in the curated-cof_XXX_vX group,
    and set the tag as the extra of the node for the query.
    """
    node.set_extra('curated-cof_tag', tag)
    group.add_nodes(node)


class ZeoppMultistageDdecPeWorkChain(WorkChain):
    """A workchain that combines: ZeoppMultistageDdecWorkChain (wc1) + IsothermCalcPEWorkChain (wc2)"""

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super(ZeoppMultistageDdecPeWorkChain, cls).define(spec)

        spec.expose_inputs(ZeoppMultistageDdecWorkChain)
        spec.expose_inputs(IsothermCalcPEWorkChain, exclude=['structure'])

        # specify the chain of calculations to be performed
        spec.outline(cls.run_wc1, cls.run_wc2, cls.return_results)

        spec.expose_outputs(ZeoppMultistageDdecWorkChain)
        spec.expose_outputs(IsothermCalcPEWorkChain)

    def run_wc1(self):
        """Run work chain 1."""
        wc1_inp = AttributeDict(self.exposed_inputs(ZeoppMultistageDdecWorkChain))
        wc1_inp['metadata']['call_link_label'] = 'call_wc1'
        running = self.submit(ZeoppMultistageDdecWorkChain, **wc1_inp)
        self.report('Running work chain 1.')
        return ToContext(wc1=running)

    def run_wc2(self):
        """Run work chain 2."""
        wc2_inp = AttributeDict(self.exposed_inputs(IsothermCalcPEWorkChain))
        wc2_inp['structure'] = self.ctx.wc1.outputs.structure_ddec
        wc2_inp['metadata']['call_link_label'] = 'call_wc2'
        running = self.submit(IsothermCalcPEWorkChain, **wc2_inp)
        return ToContext(wc2=running)

    def return_results(self):
        """Return results and include them to the group.
        To be added also the work chain nodes:
        - dftopt_wc
        - ddec_wc
        - isot_co2_wc
        - isot_n2_wc
        - pe_calc
        """
        self.out_many(self.exposed_outputs(self.ctx.wc1, ZeoppMultistageDdecWorkChain))
        self.out_many(self.exposed_outputs(self.ctx.wc2, IsothermCalcPEWorkChain))
        # This part is now very specific to CURATED-COFs (to make more flexible later)
        orig_cif = self.exposed_inputs(ZeoppMultistageDdecWorkChain)['structure']
        group_label = "curated-cof_{}_v3".format(orig_cif.label)
        group = Group(
            label=group_label,
            description="Group collecting the results of CURATED-COFs: v3 is the update of Feb 2020, adding more COFs")
        include_node("orig_cif", orig_cif, group)
        include_node("orig_zeopp_out", self.ctx.wc1.zeopp_before_opt__output_parameters, group)
        include_node("dftopt_out", self.ctx.wc1.outputs.output_parameters, group)
        self.ctx.wc1.outputs.structure_ddec.label = "{}_DDEC".format(orig_cif.label)
        include_node("opt_cif_ddec", self.ctx.wc1.outputs.structure_ddec, group)
        include_node("opt_zeopp_out", self.ctx.wc1.zeopp_after_opt__output_parameters, group)
        include_node("isot_co2_out", self.ctx.wc2.outputs.co2__output_parameters, group)
        include_node("isot_n2_out", self.ctx.wc2.outputs.co2__output_parameters, group)
        include_node("pe_out", self.ctx.wc2.outputs.calc_pe, group)
