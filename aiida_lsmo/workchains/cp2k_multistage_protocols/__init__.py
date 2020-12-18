# -*- coding: utf-8 -*-
"""Protocols for multi-stage CP2K calculations"""
from pathlib import Path
from voluptuous import Schema, Optional, Any, MultipleInvalid, Match
import ruamel.yaml as yaml  # does not convert OFF to False

__all__ = ('PROTOCOL_DIR', 'ISOTHERM_PROTOCOL_SCHEMA')

PROTOCOL_DIR = Path(__file__).resolve().parent

NUMBER = Any(int, float)
ELEMENT = Any(
    'H',
    'He',
    'Li',
    'Be',
    'B',
    'C',
    'N',
    'O',
    'F',
    'Ne',
    'Na',
    'Mg',
    'Al',
    'Si',
    'P',
    'S',
    'Cl',
    'Ar',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'Ge',
    'As',
    'Se',
    'Br',
    'Kr',
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
    'Cs',
    'Ba',
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
    'Fr',
    'Ra',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
    'Rf',
    'Db',
    'Sg',
    'Bh',
    'Hs',
    'Mt',
    'Ds',
    'Rg',
    'Cn',
    'Nh',
    'Fl',
    'Mc',
    'Lv',
    'Ts',
    'Og',
)
SETTINGS_SCHEMA = Schema({
    str: dict,
})
ISOTHERM_PROTOCOL_SCHEMA = Schema(
    {
        'protocol_description':
            str,
        Optional('initial_magnetization', default='oxidation_state'):
            Any('element', 'oxidation_state', 'zero', {ELEMENT: int}, {ELEMENT: dict}),
        'basis_set': {
            ELEMENT: str,
        },
        'pseudopotential': {
            ELEMENT: str,
        },
        'bandgap_thr_ev':
            NUMBER,
        Match(r'settings_\d+'):
            SETTINGS_SCHEMA,
        Match(r'stage_\d+'):
            SETTINGS_SCHEMA,
    },
    required=True)


def load_isotherm_protocol(*, singlefiledata=None, tag=None):
    """Load isotherm protocol from yaml file (with validation)."""

    if singlefiledata is not None:
        with singlefiledata.open() as stream:
            protocol_dict = yaml.safe_load(stream)
    elif tag is not None:
        yaml_file = PROTOCOL_DIR / (tag + '.yaml')
        with open(yaml_file, 'r') as stream:
            protocol_dict = yaml.safe_load(stream)
    else:
        raise ValueError('Provide either path or tag.')

    ISOTHERM_PROTOCOL_SCHEMA(protocol_dict)

    return protocol_dict


### Magnetization / oxidation state handling
INITIAL_MAGNETIZATION_SCHEMA = Schema(
    {
        ELEMENT: {
            'atomic_number': int,
            'magnetization': {
                int: NUMBER
            },
            'default_oxidation': Any(int, None),
            'is_metal': bool
        }
    },
    required=True)

with open(PROTOCOL_DIR / 'magnetization_data' / 'initial_magnetization.yaml') as handle:
    INITIAL_MAGNETIZATION = yaml.safe_load(handle)
    INITIAL_MAGNETIZATION_SCHEMA(INITIAL_MAGNETIZATION)


def set_initial_conditions(atoms, initial_magnetization, oxidation_states=None):  # pylint: disable=too-many-branches
    """Set initial conditions (magnetizations and oxidation states) on ASE atoms instance.

    Note: The oxidation_states argument could be eliminated by passing the oxidation states directly as charges
    on the atoms. This is currently not done because StructureData does not support charges.
    In the future, we may want to support setting initial charges via a modified StructureData.

    :param atoms: ASE atoms instance
    :param initial_magnetization: can be a string, a dictionary of element=>magnetization or a dictionary
        of element=>dictionary (see initial_magnetizations.yaml)
    :returns: atoms instance (modified in-place)
    """
    mode = initial_magnetization

    if mode in ['element', 'oxidation_state']:
        # First use magnetization table
        for atom in atoms:
            oxidation_state = INITIAL_MAGNETIZATION[atom.symbol]['default_oxidation']
            if oxidation_state is None:
                atom.charge = 0
                atom.magmom = 0
            else:
                atom.charge = -oxidation_state
                atom.magmom = INITIAL_MAGNETIZATION[atom.symbol]['magnetization'][oxidation_state]

        # If oxidation state provided, update for specified elements
        if mode == 'oxidation_state':
            for index, element, oxidation_state in zip(oxidation_states['metal_indices'],
                                                       oxidation_states['metal_symbols'],
                                                       oxidation_states['prediction']):
                atoms[index].charge = -oxidation_state
                atoms[index].magmom = INITIAL_MAGNETIZATION[element]['magnetization'][oxidation_state]

    elif mode == 'zero':
        # all zero
        for atom in atoms:
            atom.charge = 0
            atom.magmom = 0

    elif is_valid(mode, Schema({ELEMENT: NUMBER})):
        # simple format, e.g.
        # Fe: 4
        # use only those (don't try to merge with defaults)
        for atom in atoms:
            if atom.symbol not in mode:
                continue
            atom.charge = 0
            atom.magmom = mode[atom.symbol]

    elif is_valid(mode, Schema({ELEMENT: dict})):
        # complex format, e.g.
        # Fe: {'default_oxidation': 2, 'magnetization': { 1: 3, 2: 4, ...} }
        # use only those (don't try to merge with defaults)
        for atom in atoms:
            if atom.symbol not in mode:
                continue
            oxidation_state = mode[atom.symbol]['default_oxidation']
            if oxidation_state is None:
                atom.charge = 0
                atom.magmom = 0
            else:
                atom.charge = -oxidation_state
                atom.magmom = mode[atom.symbol]['magnetization'][oxidation_state]
    else:
        raise ValueError(f"Invalid 'initial_magnetization' field {mode}")

    return tag_kinds(atoms)


def is_valid(data, schema):
    """Return True, if data is valid according to schema"""
    try:
        schema(data)
        return True
    except MultipleInvalid:
        return False


def tag_kinds(atoms):
    """Tag different atom kinds, depending on oxidation state, magnetization, etc.

    E.g. if there are 4 different types of 'Fe', tag them with 1,2,3 and 4.
    """
    # Prepare dictionary [symbol][magnetization] => tag  (tag is 0, if only single oxidation state)
    symbols = set(atoms.get_chemical_symbols())

    def get_kind(atom):
        return f"{{'symbol': {atom.symbol}, 'charge': {atom.charge}, 'magmom': {atom.magmom}}}"

    element_groups = [[atom for atom in atoms if atom.symbol == symbol] for symbol in symbols]
    for group in element_groups:
        kinds = sorted({get_kind(atom) for atom in group})

        if len(kinds) == 1:
            # we don't tag if all atoms of the element have the same properties
            continue

        for atom in group:
            atom.tag = kinds.index(get_kind(atom)) + 1

    return atoms
