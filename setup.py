#!/usr/bin/env python
from setuptools import setup
import panta

setup(
    name='panta',
    version=panta.__version__,
    description='PanTA',
    packages=['panta'],
    install_requires=[],
    #scripts=['pan-genome.py'],
    entry_points={'console_scripts': ['panta=panta.cli:main']},
    zip_safe=True,
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
