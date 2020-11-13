===============
Development
===============


Running the tests
++++++++++++++++++

``aiida-lsmo`` uses the `aiida-testing <https://github.com/aiidateam/aiida-testing>`_ package in order to be able to run full integration tests of work chains without the need to have all the computational software (zeo++, raspa, cp2k, ...) installed.

As long as you are not changing the inputs for the simulation codes (as produced by AiiDA), you can run the entire test suite as follows::

    pip install -e .[testing]
    pytest

Updating the test data
+++++++++++++++++++++++

If you are changing the inputs for one or more of the simulation codes, you will need to

 1. Install the corresponding code on your machine
 2. Add an ``.aiida-testing-config.yml`` file to the top-level directory of the repository with a content like

    .. code-block:: yaml

        mock_code:
         # code-label: absolute path
         cp2k-7.1: /path/to/cp2k-7.1-Linux-x86_64.ssmp
         zeopp-0.3: /path/to/zeoplusplus/network
         raspa-e968334: /path/to/RASPA2/src/.libs/simulate

While running the tests, ``aiida-testing`` will then automatically run the simulation code for new inputs as needed and store its outputs in ``tests/data``

Please remember to:

 - Commit the new test data (as we do not install the simulation codes on CI)
 - Do **not** commit your ``.aiida-testing-config.yaml``, since the paths to the simulation codes is only valid on your computer.
