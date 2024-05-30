#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Project setup with setuptools
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from pathlib import Path

cwd = Path(__file__).parent
long_description = (cwd / "README.md").read_text()

SPECIAL_DEPENDENCIES_FOR_RELEASE = [
    "idaes-pse>=2.4.0,<2.5.0rc0",  # from PyPI
]

SPECIAL_DEPENDENCIES_FOR_PRERELEASE = [
    # update with a tag from the nawi-hub/idaes-pse
    # when a version of IDAES newer than the latest stable release from PyPI
    # will become needed for the watertap development
    "idaes-pse @ git+https://github.com/watertap-org/idaes-pse@2.5.0.dev0.watertap.24.05.09",
]

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name="watertap",
    url="https://github.com/watertap-org/watertap",
    version="1.0.dev0",
    description="WaterTAP modeling library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="NAWI team",
    license="BSD",
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="water systems, chemical engineering, process modeling, filtration, desalination, nawi",
    # just include watertap and everything under it
    packages=find_packages(
        include=("watertap*",),
    ),
    python_requires=">=3.8",
    install_requires=[
        # primary requirements for unit and property models
        # maintainers: switch to SPECIAL_DEPENDENCIES_FOR_RELEASE when cutting a release of watertap
        *SPECIAL_DEPENDENCIES_FOR_PRERELEASE,
        "pyomo>=6.6.1",  # (also needed for units in electrolyte database (edb))
        "pyyaml",  # watertap.core.wt_database
        # the following requirements are for the electrolyte database (edb)
        "pymongo>3",  # database interface
        "fastjsonschema",  # schema validation
        "click",  # command-line tools with Click
        # for parameter_sweep
        "h5py",
        "requests",
        "scipy",
        # for watertap.ui.api_model (though may be generally useful)
        "pydantic >= 2, <3",
        "numpy",
        "importlib-resources",
    ],
    extras_require={
        "testing": [
            "pytest",
            "json-schema-for-humans",
            "mongomock",
            "pandas",
            # treebeardtech/nbmake#121
            "nbmake != 1.5.1",
            "nbconvert",
        ],
        "notebooks": [
            "jupyter",
            "ipykernel",
        ],
        "oli_api": [
            "requests",
            "cryptography",  # for encrypting OLI credentials
        ],
        "dev": [
            "nbsphinx",  # jupyter notebook support for sphinx
            "jinja2<3.1.0",  # see watertap-org/watertap#449
            "Sphinx==7.1.*",  # docs
            "sphinx_rtd_theme",  # docs
            "urllib3 < 2",  # see watertap-org/watertap#1021,
            # other requirements
            "linkify-it-py",
            "json-schema-for-humans",  # pretty JSON schema in HTML
            "black",  # code formatting
            # other requirements
            "pytest",  # test framework
            "pytest-cov",  # code coverage
            "mongomock",  # mongodb mocking for testing
            # treebeardtech/nbmake#121
            "nbmake != 1.5.1",
        ],
    },
    package_data={  # Optional
        "": [
            "*.json",
            "*.yaml",
            "*.yml",
            "*.csv",
            "*.png",
        ],
        "watertap.tools.oli_api.tests": [
            "test.dbs",
        ],
    },
    entry_points={
        # add edb CLI commands
        "console_scripts": [
            "edb = watertap.edb.commands:command_base",
        ],
        "watertap.flowsheets": [
            "nf = watertap.examples.flowsheets.nf_dspmde.nf_ui",
            "nf_with_bypass = watertap.examples.flowsheets.nf_dspmde.nf_with_bypass_ui",
            "bsm2 = watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2_ui",
            "bsm2_P_extension = watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2_P_extension_ui",
            "dye_desalination = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination_ui",
            "mvc = watertap.examples.flowsheets.mvc.mvc_single_stage_ui",
            "RO = watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery_ui",
            "OARO = watertap.examples.flowsheets.oaro.oaro_multi_ui",
            "GAC = watertap.examples.flowsheets.gac.gac_ui",
        ],
    },
)
