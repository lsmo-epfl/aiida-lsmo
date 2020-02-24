"""BindingSite workchain."""

from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import WorkChain, ToContext
from aiida_lsmo.utils import get_structure_from_cif

# import sub-workchains
SimAnnealingWorkChain = WorkflowFactory('lsmo.sim_annealing')  # pylint: disable=invalid-name
Cp2kBindingEnergyWorkChain = WorkflowFactory('lsmo.cp2k_binding_energy')  # pylint: disable=invalid-name

# import aiida data
StructureData = DataFactory('structure')  # pylint: disable=invalid-name
CifData = DataFactory('cif')  # pylint: disable=invalid-name


class BindingSiteWorkChain(WorkChain):
    """A workchain that combines SimAnnealing & Cp2kBindingEnergy"""

    @classmethod
    def define(cls, spec):
        """Define workflow specification."""
        super().define(spec)

        spec.expose_inputs(SimAnnealingWorkChain)
        spec.expose_inputs(Cp2kBindingEnergyWorkChain, exclude=['structure', 'molecule'])

        spec.outline(cls.run_sim_annealing, cls.run_cp2k_binding_energy, cls.return_results)

        spec.expose_outputs(SimAnnealingWorkChain, namespace='ff')
        spec.expose_outputs(Cp2kBindingEnergyWorkChain, namespace='dft')

    def run_sim_annealing(self):
        """Run SimAnnealing"""
        sa_inputs = self.exposed_inputs(SimAnnealingWorkChain)
        sa_inputs['metadata']['call_link_label'] = 'call_sim_annealing'
        running = self.submit(SimAnnealingWorkChain, **sa_inputs)
        return ToContext(sa_wc=running)

    def run_cp2k_binding_energy(self):
        """Pass the ouptput molecule's geometry to Cp2kBindingEnergy."""
        be_inputs = self.exposed_inputs(Cp2kBindingEnergyWorkChain)
        be_inputs['structure'] = get_structure_from_cif(self.exposed_inputs(SimAnnealingWorkChain)['structure'])
        be_inputs['molecule'] = get_structure_from_cif(self.ctx.sa_wc.outputs['loaded_molecule'])
        be_inputs['metadata']['call_link_label'] = 'call_cp2k_binding_energy'

        running = self.submit(Cp2kBindingEnergyWorkChain, **be_inputs)
        return ToContext(be_wc=running)

    def return_results(self):
        """Return exposed outputs and info."""
        self.out_many(self.exposed_outputs(self.ctx.sa_wc, SimAnnealingWorkChain, namespace='ff'))
        self.out_many(self.exposed_outputs(self.ctx.be_wc, Cp2kBindingEnergyWorkChain, namespace='dft'))
        self.report("Completed! See 'ff' and 'dft' namespaces for results.")
