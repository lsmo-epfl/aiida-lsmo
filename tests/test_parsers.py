"""Test Cp2k output parsers."""
import os
from aiida_lsmo.parsers.parser_functions import parse_cp2k_output_bsse
from . import DATA_DIR


def test_bsse_parser():
    """Testing BSSE parser."""
    with open(os.path.join(DATA_DIR, 'BSSE_output_v5.1_.out')) as fobj:
        res = parse_cp2k_output_bsse(fobj)
        assert res["exceeded_walltime"] is False
        assert res["energy_description_list"] == [
            "Energy of A with basis set A", "Energy of B with basis set B", "Energy of A with basis set of A+B",
            "Energy of B with basis set of A+B", "Energy of A+B with basis set of A+B"
        ]
        assert res["energy_list"] == [
            -792.146217025347, -37.76185844385889, -792.1474366957273, -37.76259020339316, -829.920698393915
        ]
        assert res["energy_dispersion_list"] == [
            -0.15310795077994, -0.0009680299214, -0.15310795077994, -0.0009680299214, -0.16221221242898
        ]
        assert res["energy"] == -829.920698393915
        assert res["energy_units"] == "a.u."
        assert res["binding_energy_raw"] == -33.14148882373938
        assert res["binding_energy_corr"] == -28.018009583204375
        assert res["binding_energy_bsse"] == -5.123479240535005
        assert res["binding_energy_unit"] == "kJ/mol"
        assert res["binding_energy_dispersion"] == -21.361676400918867
