#!/usr/bin/env python
from setuptools import setup
import amromics

setup(
    name='amromics',
    version=amromics.__version__,
    description='AMR analysis pipelines',
    packages=['amromics', 'amromics.pipeline', 'amromics.libs','amromics.utils','amromics.db'],
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'pysam'
    ],
    scripts=['amr-analysis.py'],
    #entry_points={'console_scripts': ['amromics=scripts/amr-analysis.py']},
    zip_safe=True,
    include_package_data=True
)
