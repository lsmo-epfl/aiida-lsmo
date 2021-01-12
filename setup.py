#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Setting up aiida-lsmo for AiiDA"""

import json
from setuptools import setup


def run_setup():
    """Provide static information in setup.json such that it can be discovered automatically."""
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(packages=['aiida_lsmo'],
          long_description=open('README.md').read(),
          long_description_content_type='text/markdown',
          **kwargs)


if __name__ == '__main__':
    run_setup()
