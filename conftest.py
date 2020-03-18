"""
For pytest
initialise a test database and profile
"""
import os
import pytest
from tests import DATA_DIR

pytest_plugins = ['aiida.manage.tests.pytest_fixtures', 'aiida_testing.mock_code']  # pylint: disable=invalid-name


@pytest.fixture(scope='function', autouse=True)
def clear_database_auto(clear_database):  # pylint: disable=unused-argument
    """Automatically clear database in between tests."""


@pytest.fixture(scope='function')
def cp2k_code(mock_code_factory):
    """
    Create mocked "cp2k" code
    """
    return mock_code_factory(
        label='cp2k-7.1',
        data_dir_abspath=DATA_DIR,
        entry_point='cp2k',
        # files *not* to copy into the data directory
        ignore_files=('_aiidasubmit.sh'))


@pytest.fixture(scope='function')
def raspa_code(mock_code_factory):
    """
    Create mocked "raspa" code
    """
    return mock_code_factory(
        label='raspa-e968334',
        data_dir_abspath=DATA_DIR,
        entry_point='raspa',
        # files *not* to copy into the data directory
        ignore_files=('_aiidasubmit.sh',))


@pytest.fixture(scope='function')
def zeopp_code(mock_code_factory):
    """
    Create mocked "zeo++" code
    """
    return mock_code_factory(
        label='zeopp-0.3',
        data_dir_abspath=DATA_DIR,
        entry_point='zeopp.network',
        # files *not* to copy into the data directory
        ignore_files=('_aiidasubmit.sh', 'UFF.rad'))


@pytest.fixture(scope='function')
def mg_mof74_cifdata():
    """CifData for Mg MOF74 CIF.
    """
    from aiida.orm import CifData  # pylint: disable=import-outside-toplevel
    with open(os.path.join(DATA_DIR, 'Mg_MOF_74.cif'), 'rb') as handle:
        cif = CifData(file=handle)

    return cif
