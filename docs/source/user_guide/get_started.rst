===============
Getting started
===============

This plugin is a collection of workchains and calcfunctions that combine the use of multiple codes
(e.g., CP2K, DDEC, Raspa, Zeo++, ...) to achieve advanced automated tasks.

Installation
++++++++++++

Use the following commands to install the plugin::

    git clone https://github.com/yakutovicha/aiida-lsmo .
    cd aiida-lsmo
    pip install -e .

Note: this will install also the related plugins (e.g., `aiida-cp2k`, `aiida-raspa`, ...) if not present already,
but they need to be configured before being used in a workchain.

Usage
+++++

A quick demo of how to submit a workchain::

    verdi daemon start         # make sure the daemon is running
    cd examples
    verdi run run_IsothermWorkChain_HKUST-1.py raspa@localhost zeopp@localhost

Note that the workchain is called as::

    from aiida.plugins import WorkflowFactory
    IsothermWorkChain = WorkflowFactory('lsmo.isotherm')
