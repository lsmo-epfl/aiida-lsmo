from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.orm.data.base import Float
from aiida.work import workfunction as wf
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, if_, Outputs
from copy import deepcopy

# subworkflows
from aiida_raspa.workflows import RaspaConvergeWorkChain

# calculation objects
ZeoppCalculation = CalculationFactory('zeopp.network')

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')

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

class VolpoKh(WorkChain):
    """Workchain that computes volpo and blocking spheres: if accessible volpo>0
    it also runs a raspa widom calculation for the Henry coefficient
    """
    @classmethod
    def define(cls, spec):
        super(VolpoKh, cls).define(spec)

        # structure
        spec.input('structure', valid_type=CifData)

        # zeopp
        spec.input('zeopp_code', valid_type=Code)
        spec.input("_zeopp_options", valid_type=dict, default=None, required=False)
        spec.input("zeopp_probe_radius", valid_type=Float)
        spec.input("zeopp_atomic_radii", valid_type=SinglefileData, default=None, required=False)

        # raspa
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData)
        spec.input("_raspa_options", valid_type=dict, default=None, required=False)
        spec.input("_raspa_usecharges", valid_type=bool, default=False, required=False)

        # workflow
        spec.outline(
            cls.run_zeopp,             #computes volpo and blocks
            if_(cls.should_run_widom)( #run Widom only if porous
                cls.init_raspa_widom,  #initialize raspa calculation
                cls.run_raspa_widom,   #run raspa calculation
            ),
            cls.return_results,
        )

        spec.dynamic_output()

    def run_zeopp(self):
        """Main function that performs zeo++ VOLPO and block calculations."""
        params = {
                'ha': True,
                #100 samples / Ang^3: accurate for all the structures
                'block': [self.inputs.zeopp_probe_radius.value, 100],
                #100k samples, may need more for structures bigger than 30x30x30
                'volpo': [self.inputs.zeopp_probe_radius.value,
                          self.inputs.zeopp_probe_radius.value, 100000]
        }

        inputs = {
            'code'         : self.inputs.zeopp_code,
            'structure'    : self.inputs.structure,
            'parameters'   : NetworkParameters(dict=params).store(),
            '_options'     : self.inputs._zeopp_options,
            '_label'       : "ZeoppVolpoBlock",
        }

        # Use default zeopp atomic radii only if a .rad file is not specified
        try:
            inputs['atomic_radii'] = self.inputs.zeopp_atomic_radii
            self.report("Zeopp will use atomic radii from the .rad file")
        except:
            self.report("Zeopp will use default atomic radii")

        # Create the calculation process and launch it
        running = submit(ZeoppCalculation.process(), **inputs)
        self.report("pk: {} | Running zeo++ volpo and block calculations".format(running.pid))
        return ToContext(zeopp=Outputs(running))

    def should_run_widom(self):
        """Submit widom calculation only if there is some accessible volume."""
        if self.ctx.zeopp['output_parameters'].dict.POAV_Volume_fraction > 1e-5:
            self.report("Found accessible pore volume: continue")
            return True
        else:
            self.report("NOT Found any accessible pore volume: stop")
            return False

    def init_raspa_widom(self):
        """ Write the ParameterData for a Raspa and use the blocking spheres.
        REMEMBER that we want only one component for widom!
        """
        # Build the default parameters dict
        raspa_widom_default_param = {
                "GeneralSettings":
                {
                "SimulationType"                   : "MonteCarlo",
                'NumberOfInitializationCycles'     : 0,
                "NumberOfCycles"                   : 10000,
                "PrintPropertiesEvery"             : 1000,  # info on henry coeff
                "PrintEvery"                       : 0,     # info on loading
                "Forcefield"                       : "LSMO_UFF-TraPPE",
                "CutOff"                           : 13.0,
                "ExternalTemperature"              : 298.0,
                },
                "Component":
                [{ #TODO: make it as a dictionary
                "MoleculeName"                     : "XXX",
                "MoleculeDefinition"               : "XXX",
                "WidomProbability"                 : 1.0,
                }],
        }
        # Turn on charges if requested
        if self.inputs._raspa_usecharges:
            raspa_widom_default_param['ChargeMethod'] = "Ewald"
            raspa_widom_default_param['EwaldPrecision'] = 1e-6
            raspa_widom_default_param['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
        else:
            raspa_widom_default_param['ChargeMethod'] = "None"

        # Merge the user with default parameters, giving priority to the first
        self.ctx.raspa_parameters = deepcopy(raspa_widom_default_param)
        dict_merge(self.ctx.raspa_parameters,
                   self.inputs.raspa_parameters.get_dict())

        # Create and store the ParameterData
        self.ctx.raspa_parameterdata = ParameterData(dict=self.ctx.raspa_parameters).store()

    def run_raspa_widom(self):
        """Run a Widom calculation in Raspa"""

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.inputs.structure,
            'parameters' : self.ctx.raspa_parameterdata,
            '_options'   : self.inputs._raspa_options,
            '_label'     : "RaspaWidom",
        }

        # Check if there are blocking spheres (reading the header of the file) and use them for Raspa
        bs = self.ctx.zeopp['block']
        with open(bs.get_abs_path() + '/path/' + bs.get_folder_list()[0]) as f:
            self.ctx.number_blocking_spheres = int(f.readline().strip())
        if self.ctx.number_blocking_spheres > 0:
            inputs['block_component_0'] = bs
            self.report("Blocking spheres ({}) are present and used for Raspa".format(self.ctx.number_blocking_spheres))
        else:
            self.report("No blocking spheres found")
        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the Henry coefficients".format(running.pid))

        return ToContext(raspa_widom=Outputs(running))


    def return_results(self):
        """Attach the results and the initial structure to the outputs."""
        result_dict = {}
        result_dict['POAV_Volume_fraction'] = self.ctx.zeopp['output_parameters'].dict.POAV_Volume_fraction
        try:
            result_dict['number_blocking_spheres'] = self.ctx.number_blocking_spheres
        except AttributeError:
            pass

        try:
            result_dict['henry_coefficient_average'] = self.ctx.raspa_widom["component_0"].dict.henry_coefficient_average
            result_dict['henry_coefficient_dev'] = self.ctx.raspa_widom["component_0"].dict.henry_coefficient_dev
            result_dict['henry_coefficient_units'] = self.ctx.raspa_widom["component_0"].dict.henry_coefficient_units
            result_dict['adsorption_energy_average'] = self.ctx.raspa_widom["component_0"].dict.adsorption_energy_widom_average
            result_dict['adsorption_energy_dev'] = self.ctx.raspa_widom["component_0"].dict.adsorption_energy_widom_dev
            result_dict['adsorption_energy_units'] = self.ctx.raspa_widom["component_0"].dict.adsorption_energy_widom_units
        except AttributeError:
            pass

        self.out("results", ParameterData(dict=result_dict).store())
        self.out('blocking_spheres', self.ctx.zeopp['block'])
        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return

# EOF
