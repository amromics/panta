#!/usr/bin/env python
from setuptools import setup
import pan_genome

setup(
    name='panta',
    version=pan_genome.__version__,
    description='PanTA',
    packages=['pan_genome'],
    install_requires=[],
    #scripts=['pan-genome.py'],
    entry_points={'console_scripts': ['panta=pan_genome.panta:main']},
    zip_safe=True,
    include_package_data=True
)
