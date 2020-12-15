#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

short_descr = "ASTEC package"
readme = open('README').read()

# find packages
pkgs = find_packages('src')

setup_kwds = dict(
    name='ASTEC',
    version="0.0.1",
    description=short_descr,
    long_description=readme,
    author="Gregoire Malandain",
    author_email="gregoire.malandain@inria.fr",
    url='https://github.com/astec-segmentation/astec',
    license='GPL',
    zip_safe=False,

    packages=pkgs,
    
    package_dir={'': 'src'},
    setup_requires=[],
    install_requires=[],
    tests_require=[],
    entry_points={},
    keywords='',
    
    test_suite='pytest',
    )

setup_kwds['entry_points']['console_scripts'] = []

setup(**setup_kwds)
