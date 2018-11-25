from aiida.orm import CalculationFactory, DataFactory
from aiida.orm import load_node
from aiida.orm.code import Code
from aiida.orm.data.base import Bool, Str, Float
from aiida.work import workfunction as wf
from aiida.work.run import run, submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, Outputs

# subworkflows
from aiida_cp2k.workflows import Cp2kRobustGeoOptWorkChain
from aiida_ddec.workflows import DdecCp2kChargesWorkChain
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

spin = {
        "H"  : 0.0,
        "Li" : 0.0,
        "Be" : 0.0,
        "B"  : 0.0,
        "C"  : 0.0,
        "N"  : 0.0,
        "O"  : 0.0,
        "F"  : 0.0,
        "Na" : 0.0,
        "Mg" : 0.0,
        "Al" : 0.0,
        "Si" : 0.0,
        "P"  : 0.0,
        "S"  : 0.0,
        "Cl" : 0.0,
        "K"  : 0.0,
        "Ca" : 0.0,
        "Sc" : 1.0 / 2.0,  # oxidation state +2, high spin
        "Ti" : 2.0 / 2.0,  # oxidation state +2, high spin
        "V"  : 3.0 / 2.0,  # oxidation state +2, high spin
        "Cr" : 4.0 / 2.0,  # oxidation state +2, high spin
        "Mn" : 5.0 / 2.0,  # oxidation state +2, high spin
        "Fe" : 4.0 / 2.0,  # oxidation state +2, high spin
        "Co" : 3.0 / 2.0,  # oxidation state +2, high spin
        "Ni" : 2.0 / 2.0,  # oxidation state +2, high spin
        "Cu" : 1.0 / 2.0,  # oxidation state +2, high spin
        "Zn" : 0.0,        # oxidation state +2, high spin
        "Zr" : 2.0 / 2.0,  # oxidation state +2, high spin
        }

@wf
def guess_multiplicity(structure):
    multiplicity = 1
    all_atoms = structure.get_ase().get_chemical_symbols()
    for key, value in spin.iteritems():
        multiplicity += all_atoms.count(key) * value * 2.0
    multiplicity = int(round(multiplicity))
    multiplicity_dict = {'FORCE_EVAL': {'DFT': {'MULTIPLICITY' :multiplicity}}}
    if multiplicity != 1:
        multiplicity_dict['FORCE_EVAL']['DFT']['UKS'] = True
    return ParameterData(dict =multiplicity_dict)

@wf
def from_cif_to_structuredata(cif, threshold):
    """Helper function that converts CifData object into StructureData"""
    repeat = multiply_unit_cell(cif, threshold=threshold.value)
    structure = cif.get_ase().repeat(repeat)
    return StructureData(ase=structure).store()


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

        # cp2k
        spec.input('cp2k_code', valid_type=Code)
        spec.input("_cp2k_options", valid_type=dict, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)

        # ddec
        spec.input('ddec_code', valid_type=Code)
        spec.input("_ddec_options", valid_type=dict, default=None, required=False)

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
        spec.input('_guess_multiplicity', valid_type=bool, default=False)

        # workflow
        spec.outline(
            cls.init,                    #read pressures, switch on cif charges if _usecharges=True
            cls.run_geo_opt,             #robust 4 stesps cellopt, first expands the uc according to min_cell_size 
            cls.parse_geo_opt,           
            if_(cls.should_use_charges)( #if charges are needed, wfn > (E_DENSITY & core_el) > DDEC charges
                cls.run_point_charges,   
                cls.parse_point_charges,
            ),
            cls.run_geom_zeopp,          #computes sa, vol, povol, res, e chan. Block pore?
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

    def run_geo_opt(self):
        """Optimize the geometry using the robust 4 steps process"""

        # Expand the unit cell so that: min(perpendicular_width) > threshold
        threshold = self.inputs.min_cell_size
        new_structure = from_cif_to_structuredata(cif=self.ctx.structure, threshold=threshold)

        inputs = {
            'code'      : self.inputs.cp2k_code,
            'structure' : new_structure,
            '_options'  : self.inputs._cp2k_options,
            '_label'    : "Cp2kRobustGeoOptWorkChain",
        }

        # Trying to guess the multiplicity of the system
        if self.inputs._guess_multiplicity:
            self.report("Guessing multiplicity")
            self.ctx.cp2k_parameters = guess_multiplicity(new_structure)
            inputs['parameters'] = self.ctx.cp2k_parameters

        # uncomment for the test runs
        # FROM HERE
        #params_dict = ParameterData(dict={
        #        'MOTION':{
        #            'MD':{
        #                'STEPS': 5,
        #                },
        #            'GEO_OPT': {
        #                'MAX_ITER': 5,
        #            },
        #            'CELL_OPT': {
        #                'MAX_ITER': 5,
        #            },
        #        },
        #        }).store()
        #inputs['parameters'] = merge_ParameterData(self.ctx.cp2k_parameters, params_dict)
        # TILL HERE

        # Create the calculation process and launch it
        running = submit(Cp2kRobustGeoOptWorkChain, **inputs)
        self.report("pk: {} | Running Cp2kRobustGeoOptWorkChain to optimize geometry".format(running.pid))
        return ToContext(geo_opt_calc=Outputs(running))

    def parse_geo_opt(self):
        """Extract optimized structure and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.geo_opt_calc['output_structure']

    def should_use_charges(self):
        """Whether it is needed to employ charges."""
        return self.inputs._usecharges

    def run_point_charges(self):
        """Compute the charge-density of a structure that can be later
        used for extracting ddec point charges."""

        inputs = {
            'structure'          : self.ctx.structure,
            'cp2k_code'          : self.inputs.cp2k_code,
            'cp2k_parameters'    : self.ctx.cp2k_parameters,
            '_cp2k_options'      : self.inputs._cp2k_options,
            'cp2k_parent_folder' : self.ctx.geo_opt_calc['remote_folder'],
            # ddec parameters
            'ddec_code'          : self.inputs.ddec_code,
            '_ddec_options'      : self.inputs._ddec_options,
            '_label'             : "DdecCp2kChargesWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(DdecCp2kChargesWorkChain, **inputs)
        self.report("pk: {} | Running DdecCp2kChargesWorkChain to compute the point charges".format(running.pid))
        return ToContext(point_charges_calc=Outputs(running))

    def parse_point_charges(self):
        """Extract structure with charges and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.point_charges_calc['output_structure']

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
