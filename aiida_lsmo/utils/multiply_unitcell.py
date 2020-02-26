"""Utilities for unit cell multiplication, typically for cut-off issues."""

from aiida.engine import calcfunction
from aiida.orm import Dict, StructureData


@calcfunction
def check_resize_unit_cell_legacy(struct, threshold):  #pylint: disable=too-many-locals
    """Returns the multiplication factors for the cell vectors to respect, in every direction:
    min(perpendicular_width) > threshold.
    TODO: this has been used for CP2K, make it uniform to the other one used for Raspa (from CifFile).
    """
    from math import cos, sin, sqrt, fabs, ceil
    import numpy as np

    # Parsing structure's cell
    def angle(vect1, vect2):
        return np.arccos(np.dot(vect1, vect2) / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))

    a_len = np.linalg.norm(struct.cell[0])
    b_len = np.linalg.norm(struct.cell[1])
    c_len = np.linalg.norm(struct.cell[2])

    alpha = angle(struct.cell[1], struct.cell[2])
    beta = angle(struct.cell[0], struct.cell[2])
    gamma = angle(struct.cell[0], struct.cell[1])

    # Computing triangular cell matrix
    vol = np.sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2 * cos(alpha) * cos(beta) * cos(gamma))
    cell = np.zeros((3, 3))
    cell[0, :] = [a_len, 0, 0]
    cell[1, :] = [b_len * cos(gamma), b_len * sin(gamma), 0]
    cell[2, :] = [
        c_len * cos(beta), c_len * (cos(alpha) - cos(beta) * cos(gamma)) / (sin(gamma)), c_len * vol / sin(gamma)
    ]
    cell = np.array(cell)

    # Computing perpendicular widths, as implemented in Raspa
    # for the check (simplified for triangular cell matrix)
    axc1 = cell[0, 0] * cell[2, 2]
    axc2 = -cell[0, 0] * cell[2, 1]
    bxc1 = cell[1, 1] * cell[2, 2]
    bxc2 = -cell[1, 0] * cell[2, 2]
    bxc3 = cell[1, 0] * cell[2, 1] - cell[1, 1] * cell[2, 0]
    det = fabs(cell[0, 0] * cell[1, 1] * cell[2, 2])
    perpwidth = np.zeros(3)
    perpwidth[0] = det / sqrt(bxc1**2 + bxc2**2 + bxc3**2)
    perpwidth[1] = det / sqrt(axc1**2 + axc2**2)
    perpwidth[2] = cell[2, 2]

    #prevent from crashing if threshold.value is zero
    if threshold.value == 0:
        thr = 0.1
    else:
        thr = threshold.value

    resize = {
        'nx': int(ceil(thr / perpwidth[0])),
        'ny': int(ceil(thr / perpwidth[1])),
        'nz': int(ceil(thr / perpwidth[2]))
    }
    return Dict(dict=resize)


@calcfunction
def resize_unit_cell(struct, resize):
    """Resize the StructureData according to the resize Dict"""
    resize_tuple = tuple([resize[x] for x in ['nx', 'ny', 'nz']])
    return StructureData(ase=struct.get_ase().repeat(resize_tuple))


def check_resize_unit_cell(cif, threshold):  #pylint: disable=too-many-locals
    """Returns the multiplication factors for the cell vectors to respect, in every direction:
    min(perpendicular_width) > threshold."""
    from math import cos, sin, sqrt, fabs, ceil, pi
    import numpy as np

    deg2rad = pi / 180.

    # Parsing cif
    struct = list(cif.values.dictionary.values())[0]

    a_len = float(struct['_cell_length_a'])
    b_len = float(struct['_cell_length_b'])
    c_len = float(struct['_cell_length_c'])

    alpha = float(struct['_cell_angle_alpha']) * deg2rad
    beta = float(struct['_cell_angle_beta']) * deg2rad
    gamma = float(struct['_cell_angle_gamma']) * deg2rad

    # Computing triangular cell matrix
    vol = np.sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2 * cos(alpha) * cos(beta) * cos(gamma))
    cell = np.zeros((3, 3))
    cell[0, :] = [a_len, 0, 0]
    cell[1, :] = [b_len * cos(gamma), b_len * sin(gamma), 0]
    cell[2, :] = [
        c_len * cos(beta), c_len * (cos(alpha) - cos(beta) * cos(gamma)) / (sin(gamma)), c_len * vol / sin(gamma)
    ]
    cell = np.array(cell)

    # Computing perpendicular widths, as implemented in Raspa
    # for the check (simplified for triangular cell matrix)
    axc1 = cell[0, 0] * cell[2, 2]
    axc2 = -cell[0, 0] * cell[2, 1]
    bxc1 = cell[1, 1] * cell[2, 2]
    bxc2 = -cell[1, 0] * cell[2, 2]
    bxc3 = cell[1, 0] * cell[2, 1] - cell[1, 1] * cell[2, 0]
    det = fabs(cell[0, 0] * cell[1, 1] * cell[2, 2])
    perpwidth = np.zeros(3)
    perpwidth[0] = det / sqrt(bxc1**2 + bxc2**2 + bxc3**2)
    perpwidth[1] = det / sqrt(axc1**2 + axc2**2)
    perpwidth[2] = cell[2, 2]

    #prevent from crashing if threshold value is zero
    thr = max(0.001, threshold)

    return int(ceil(thr / perpwidth[0])), int(ceil(thr / perpwidth[1])), int(ceil(thr / perpwidth[2]))
