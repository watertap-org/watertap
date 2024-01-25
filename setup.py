#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
    "idaes-pse>=2.3.0,<2.4.0rc0",  # from PyPI
]

SPECIAL_DEPENDENCIES_FOR_PRERELEASE = [
    # update with a tag from the nawi-hub/idaes-pse
    # when a version of IDAES newer than the latest stable release from PyPI
    # will become needed for the watertap development
    "idaes-pse>=2.3.0,<2.4.0rc0",
]

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name="watertap",
    url="https://github.com/watertap-org/watertap",
    version="0.12.dev0",
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
        "pydantic < 2",
        "numpy",
        "importlib-resources",
    ],
    extras_require={
        "testing": [
            "pytest",
            "json-schema-for-humans",
            "mongomock",
            "pandas",
            "nbmake",
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
            "nbmake",
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
            "metab = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab_ui",
            "suboxic_ASM = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.suboxic_activated_sludge_process.suboxic_ASM_ui",
            "Magprex = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1575_magprex.magprex_ui",
            "biomembrane_filtration = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.biomembrane_filtration.biomembrane_filtration_ui",
            "ENR = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.electrochemical_nutrient_removal.electrochemical_nutrient_removal_ui",
            "CANDO_P = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1595_photothermal_membrane_candoP.amo_1595_ui",
            "supercritical_sludge_to_gas = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.supercritical_sludge_to_gas.supercritical_sludge_to_gas_ui",
            "PAA = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.peracetic_acid_disinfection.peracetic_acid_disinfection_ui",
            "AMO 1690 = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1690.amo_1690_ui",
            "HRCS = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1575_hrcs.hrcs_ui",
            "groundwater_treatment = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.groundwater_treatment.groundwater_treatment_ui",
            "dye_desalination = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination_ui",
            "swine_wwt = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.swine_wwt.swine_wwt_ui",
            "GLSD anaerobic digestion = watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.GLSD_anaerobic_digester.GLSD_anaerobic_digestion_ui",
            "mvc = watertap.examples.flowsheets.mvc.mvc_single_stage_ui",
            "RO = watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery_ui",
            "GAC = watertap.examples.flowsheets.gac.gac_ui",
        ],
    },
)
