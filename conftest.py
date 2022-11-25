# -*- coding: utf-8 -*-
"""
For pytest
initialise a test database and profile
"""
import os
import pytest
from aiida_ddec.calculations import DENSITY_DIR_EXTRA, DENSITY_DIR_SYMLINK

from tests import DATA_DIR
from tests.cache_hashing import CustomInputHasher
from examples import DATA_DIR as EXAMPLES_DATA_DIR


@pytest.fixture(scope='function', autouse=True)
def clear_database_auto(clear_database):  # pylint: disable=unused-argument
    """Automatically clear database in between tests."""


@pytest.fixture(scope='function')
def cp2k_code(mock_code_factory):
    """Create mocked "cp2k" code."""
    return mock_code_factory(
        label='cp2k-7.1',
        data_dir_abspath=DATA_DIR,
        entry_point='cp2k',
        hasher=CustomInputHasher,
        # files *not* to copy into the data directory
        ignore_paths=('_aiidasubmit.sh', 'BASIS_MOLOPT', 'GTH_POTENTIALS', 'dftd3.dat', '*.bak*'))


@pytest.fixture(scope='function')
def raspa_code(mock_code_factory):
    """Create mocked "raspa" code."""
    return mock_code_factory(
        label='raspa-e968334',
        data_dir_abspath=DATA_DIR,
        entry_point='raspa',
        # paths *not* to copy into the data directory
        ignore_paths=('_aiidasubmit.sh', 'CrashRestart/*', 'Movies/*', 'VTK/*', 'RestartInitial/*'))


@pytest.fixture(scope='function')
def zeopp_code(mock_code_factory):
    """Create mocked "zeo++" code."""
    return mock_code_factory(
        label='zeopp-0.3',
        data_dir_abspath=DATA_DIR,
        entry_point='zeopp.network',
        # files *not* to copy into the data directory
        ignore_paths=('_aiidasubmit.sh', 'UFF.rad'))


@pytest.fixture(scope='function')
def ddec_code(mock_code_factory):
    """Create mocked "ddec" code."""
    code = mock_code_factory(
        label='chargemol-09_26_2017',
        data_dir_abspath=DATA_DIR,
        entry_point='ddec',
        # files *not* to copy into the data directory
        ignore_paths=('_aiidasubmit.sh', '*.cube', DENSITY_DIR_SYMLINK))

    # Set atomic density directory extra on code
    density_dir = os.environ.get(DENSITY_DIR_EXTRA)
    if not density_dir:
        density_dir = EXAMPLES_DATA_DIR / 'ddec' / 'atomic_densities'
    code.set_extra(DENSITY_DIR_EXTRA, str(density_dir))

    return code
