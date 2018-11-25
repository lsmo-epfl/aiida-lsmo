from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.orm.data.base import Float
from aiida.work import workfunction as wf
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, Outputs

# subworkflows
from aiida_raspa.workflows import RaspaConvergeWorkChain
from aiida_zeopp.workflows import ZeoppBlockPocketsWorkChain

from copy import deepcopy

def dict_merge(dct, merge_dct):
    """ Taken from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    import collections
    for k, v in merge_dct.iteritems():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

@wf
def merge_ParameterData(p1, p2):
    p1_dict = p1.get_dict()
    p2_dict = p2.get_dict()
    dict_merge(p1_dict, p2_dict)
    return ParameterData(dict=p1_dict).store()

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')

def multiply_unit_cell (cif, threshold):
    """Resurns the multiplication factors (tuple of 3 int) for the cell vectors
    that are needed to respect: min(perpendicular_width) > threshold
    """
    from math import cos, sin, sqrt, pi
    import numpy as np
    deg2rad=pi/180.
    threshold /= 2.0

    struct=cif.values.dictionary.itervalues().next()

    a = float(struct['_cell_length_a'])
    b = float(struct['_cell_length_b'])
    c = float(struct['_cell_length_c'])

    alpha = float(struct['_cell_angle_alpha'])*deg2rad
    beta  = float(struct['_cell_angle_beta'])*deg2rad
    gamma = float(struct['_cell_angle_gamma'])*deg2rad

    # first step is computing cell parameters according to  https://en.wikipedia.org/wiki/Fractional_coordinates
    # Note: this is the algorithm implemented in Raspa (framework.c/UnitCellBox). There also is a simpler one but it is less robust.
    v = sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
    cell=np.zeros((3,3))
    cell[0,:] = [a, 0, 0]
    cell[1,:] = [b*cos(gamma), b*sin(gamma),0]
    cell[2,:] = [c*cos(beta), c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma)),c*v/sin(gamma)]
    cell=np.array(cell)

    # diagonalizing the cell matrix: note that the diagonal elements are the perpendicolar widths because ay=az=bz=0
    diag = np.diag(cell)
    return tuple(int(i) for i in np.ceil(threshold/diag*2.))

class Isotherm(WorkChain):
    """Workchain that for a given matherial will compute an isotherm of a
    certain gaz adsorption."""
    @classmethod
    def define(cls, spec):
        super(Isotherm, cls).define(spec)

        # structure, adsorbant, pressures
        spec.input('structure', valid_type=CifData)
        spec.input("probe_molecule", valid_type=ParameterData)
        spec.input("pressures", valid_type=ArrayData)
        spec.input("min_cell_size", valid_type=Float)

        # zeopp
        spec.input('zeopp_code', valid_type=Code)
        spec.input("_zeopp_options", valid_type=dict, default=None, required=False)

        # raspa
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData)
        spec.input("_raspa_options", valid_type=dict, default=None, required=False)

        # settings
        spec.input("_interactive", valid_type=bool, default=False, required=False)
        spec.input("_usecharges", valid_type=bool, default=False, required=False)

        # workflow
        spec.outline(
            cls.init,                    #read pressures, switch on cif charges if _usecharges=True
            cls.run_block_pockets,       #computes sa, vol, povol, res, e chan. Block pore?
            cls.init_raspa_calc,         #assign HeliumVoidFraction=POAV and UnitCells
            cls.run_henry_raspa,         #NumberOfInitializationCycles=0, remove ExternalPressure, WidomProbability=1
            while_(cls.should_run_loading_raspa)( 
                cls.run_loading_raspa,   #for each P, recover the last snapshoot of the previous and run GCMC
                cls.parse_loading_raspa,
            ),
            cls.return_results,
        )

        # TODO: once the workflow is ready, explicitely specify the outputs
        spec.dynamic_output()

    def init(self):
        """Initialize variables and the pressures we want to compute"""
        self.ctx.structure = self.inputs.structure
        self.ctx.pressures = self.inputs.pressures.get_array("pressures")
        self.ctx.current_p_index = 0
        self.ctx.result = []

        self.ctx.cp2k_parameters = None

        self.ctx.raspa_parameters = self.inputs.raspa_parameters.get_dict()

        if self.inputs._usecharges:
            self.ctx.raspa_parameters['ChargeMethod'] = "Ewald"
            self.ctx.raspa_parameters['EwaldPrecision'] = 1e-6
            self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
        self.ctx.restart_raspa_calc = None


    def run_geom_zeopp(self):
        """Zeo++ calculation for geometric properties"""
        # network parameters
        sigma = self.inputs.probe_molecule.dict.sigma

        # Create the input dictionary
        inputs = {
            'probe_radius' : Float(sigma),
            'structure'    : self.ctx.structure,
            'zeopp_code'   : self.inputs.zeopp_code,
            '_options'     : self.inputs._zeopp_options,
            '_label'       : "ZeoppBlockPocketsWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(ZeoppBlockPocketsWorkChain, **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(running.pid))

        return ToContext(zeopp=Outputs(running))

    def init_raspa_calc(self):
        """Parse the output of Zeo++ and instruct the input for Raspa. """
        # Use probe-occupiable available void fraction as the helium void fraction (for excess uptake)
        self.ctx.raspa_parameters['GeneralSettings']['HeliumVoidFraction'] = \
        self.ctx.zeopp["output_parameters"].dict.POAV_Volume_fraction
        # Compute the UnitCells expansion considering the CutOff
        cutoff = self.ctx.raspa_parameters['GeneralSettings']['CutOff']                           
        ucs = multiply_unit_cell(self.ctx.structure,2*cutoff)
        self.ctx.raspa_parameters['GeneralSettings']['UnitCells'] = "{} {} {}".format(ucs[0], ucs[1], ucs[2])

    def run_henry_raspa(self):
        """Run a Widom insertion calculation in Raspa"""
        raspa_parameters_widom = deepcopy(self.ctx.raspa_parameters)

        # Remove the pressure and InitializationCycles (not needed for Widom insertion)
        raspa_parameters_widom['GeneralSettings'].pop('ExternalPressure')
        raspa_parameters_widom['GeneralSettings']['NumberOfInitializationCycles'] = 0

        # Switch the settings to Widom insertion
        for i, comp in enumerate(raspa_parameters_widom['Component']):
            name = comp['MoleculeName']
            raspa_parameters_widom['Component'][0] = {
                "MoleculeName"                     : name,
                "MoleculeDefinition"               : "TraPPE",
                "WidomProbability"                 : 1.0,
                "CreateNumberOfMolecules"          : 0,
            }

        parameters = ParameterData(dict=raspa_parameters_widom).store()

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : parameters,
            '_options'   : self.inputs._raspa_options,
            '_label'     : "RaspaConvergeWorkChain",
        }

        # Check if there are poket blocks to be loaded
        try:
            inputs['block_component_0'] = self.ctx.zeopp['block']
        except:
            pass

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the Henry coefficients".format(running.pid))

        return ToContext(raspa_henry=Outputs(running))

    def should_run_loading_raspa(self):
        """We run another raspa calculation only if the current iteration is smaller than
        the total number of pressures we want to compute."""
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_loading_raspa(self):
        """This function will run RaspaConvergeWorkChain for the current pressure"""
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure'] = pressure

        parameters = ParameterData(dict=self.ctx.raspa_parameters).store()
        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : parameters,
            '_options'    : self.inputs._raspa_options,
            '_label'     : "run_loading_raspa",
        }
        # Check if there are poket blocks to be loaded
        try:
            inputs['block_component_0'] = self.ctx.zeopp['block']
        except:
            pass

        if self.ctx.restart_raspa_calc is not None:
            inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the pressure {} [bar]".format(running.pid, pressure/1e5))
        self.ctx.current_p_index += 1
        return ToContext(raspa_loading=Outputs(running))

    def parse_loading_raspa(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""
        self.ctx.restart_raspa_calc = self.ctx.raspa_loading['retrieved_parent_folder']
        pressure = self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure']/1e5
        loading_average = self.ctx.raspa_loading["component_0"].dict.loading_absolute_average
        self.ctx.result.append((pressure, loading_average))

    def return_results(self):
        """Attach the results of the raspa calculation and the initial structure to the outputs."""
        self.out("result", ParameterData(dict={"isotherm": self.ctx.result}).store())
        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return

# EOF
