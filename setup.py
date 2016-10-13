#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Setup script for diffusiond

To install, run:

python setup.py install

"""

from setuptools import setup, find_packages

setup(
    name='diffusiond',
    version='0.1.0',
    description='Heat diffusion daemon',
    long_description='Sub-network finder using heat diffusion simulation.',
    author='Daniel Carlin, Keiichiro Ono, Eric Sage',
    author_email='eric.david.sage@gmail.com',
    url='https://github.com/ericsage/diffusiond',
    license='MIT License',
    scripts=['diffuse.py'],
    packages=['diffusiond'],
    install_requires=[
        'bottle',
        'ndex',
        'numpy',
        'scipy'
    ],
    keywords=['bioinformatics', 'graph', 'network', 'cytoscape'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    include_package_data=True,
)
