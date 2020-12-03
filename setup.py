# -*- coding: utf-8 -*-
from __future__ import absolute_import

from setuptools import setup, find_packages

import versioneer

setup(
    name="mof_vac_level",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    url="",
    license="MIT",
    install_requires=["numpy", "scipy", "matplotlib", "tqdm", "multiprocess", "numba"],
    extras_require={
        "testing": ["pytest", "pytest-cov<2.11"],
        "docs": ["sphinx-rtd-theme", "sphinxcontrib-bibtex"],
        "pre-commit": ["pre-commit", "yapf", "prospector", "pylint", "versioneer"],
    },
    author="Kevin M. Jablonka, Maria Fumanal, Berend Smit",
    author_email="kevin.jablonka@epfl.ch",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 1 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
