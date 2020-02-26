"""ZeoppMultistageDdecPeWorkChain workchain"""

from aiida.orm import Group
from aiida.plugins import WorkflowFactory
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext

# import sub-workchains
ZeoppMultistageDdecWorkChain = WorkflowFactory('lsmo.zeopp_multistage_ddec')  # pylint: disable=invalid-name
IsothermCalcPEWorkChain = WorkflowFactory('lsmo.isotherm_calc_pe')  # pylint: disable=invalid-name


def include_node(tag, node, group):
    """Given an aiida-node and a (string) tag, add the node in the curated-cof_XXX_vX group,
    and set the tag as the extra of the node for the query.
    """
    node.set_extra('curated-cof_tag', tag)
    group.add_nodes(node)


class NanoporousScreening1WorkChain(WorkChain):
    """A workchain that combines: ZeoppMultistageDdecWorkChain wc1 and IsothermCalcPEWorkChain wc2.
    In future I will use this to include more applications to run in parallel.
    """

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super().define(spec)

        spec.expose_inputs(ZeoppMultistageDdecWorkChain)
        spec.expose_inputs(IsothermCalcPEWorkChain, exclude=['structure'])

        # specify the chain of calculations to be performed
        spec.outline(cls.make_group, cls.run_wc1, cls.include_results_wc1, cls.run_wc2, cls.include_results_wc2)

        # do not check for output
        spec.outputs.dynamic = True

    def make_group(self):
        """Create curated-xxx_XXXX_vx group and put the orig_cif inside, and exit if it already exists."""
        self.ctx.orig_cif = self.exposed_inputs(ZeoppMultistageDdecWorkChain)['structure']
        self.ctx.group = Group(
            label="curated-{}_{}_v3".format(self.ctx.orig_cif.extras["class_material"], self.ctx.orig_cif.label),
            description=
            "Group collecting the results of CURATED-COFs/MOFs/ZEOs: v3 is consistent with the API of Feb 2020")
        self.ctx.group.store()  # REMEMBER: this will crash if a node with the same label exists!
        include_node("orig_cif", self.ctx.orig_cif, self.ctx.group)

    def run_wc1(self):
        """Run work chain 1."""
        wc1_inp = AttributeDict(self.exposed_inputs(ZeoppMultistageDdecWorkChain))
        wc1_inp['metadata']['call_link_label'] = 'call_wc1'
        running = self.submit(ZeoppMultistageDdecWorkChain, **wc1_inp)
        self.report('Running work chain 1.')
        return ToContext(wc1=running)

    def include_results_wc1(self):
        """Include results of work chain 1 in group."""
        include_node("orig_zeopp_out", self.ctx.wc1.outputs.zeopp_before_opt__output_parameters, self.ctx.group)
        include_node("dftopt_out", self.ctx.wc1.outputs.output_parameters, self.ctx.group)
        self.ctx.wc1.outputs.structure_ddec.label = "{}_DDEC".format(self.ctx.orig_cif.label)
        include_node("opt_cif_ddec", self.ctx.wc1.outputs.structure_ddec, self.ctx.group)
        include_node("opt_zeopp_out", self.ctx.wc1.outputs.zeopp_after_opt__output_parameters, self.ctx.group)
        include_node("dftopt_wc", self.ctx.wc1.called[-2], self.ctx.group)
        include_node("ddec_wc", self.ctx.wc1.called[-3], self.ctx.group)

    def run_wc2(self):
        """Run work chain 2."""
        wc2_inp = AttributeDict(self.exposed_inputs(IsothermCalcPEWorkChain))
        wc2_inp['structure'] = self.ctx.wc1.outputs.structure_ddec
        wc2_inp['metadata']['call_link_label'] = 'call_wc2'
        running = self.submit(IsothermCalcPEWorkChain, **wc2_inp)
        return ToContext(wc2=running)

    def include_results_wc2(self):
        """Include results of work chain 2 in group."""
        include_node("isot_co2_out", self.ctx.wc2.outputs.co2__output_parameters, self.ctx.group)
        include_node("isot_n2_out", self.ctx.wc2.outputs.n2__output_parameters, self.ctx.group)
        include_node("pe_out", self.ctx.wc2.outputs.output_parameters, self.ctx.group)
        include_node("isot_co2_wc", self.ctx.wc2.called[-1], self.ctx.group)
        include_node("isot_n2_wc", self.ctx.wc2.called[-2], self.ctx.group)
