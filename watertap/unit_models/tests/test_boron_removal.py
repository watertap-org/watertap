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
import pytest
import idaes.core.util.scaling as iscale
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.boron_removal import BoronRemoval
from pyomo.environ import (
    ConcreteModel,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.exceptions import ConfigurationError
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
from watertap.core.solvers import get_solver
import re

__author__ = "Austin Ladshaw"

solver = get_solver()


# Helper function for multiple test setup
def model_setup(
    m,
    chem_add=1.8e-5,
    state={
        "H2O": 100,
        "H_+": 1e-7,
        "OH_-": 1e-7,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-4,
        "HCO3_-": 1e-4,
    },
):

    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.unit.inlet.flow_mol_phase_comp[idx].fix(state[j])
    m.fs.unit.caustic_dose_rate.fix(chem_add)
    m.fs.unit.reactor_volume.fix(1)


# Helper function to automate scaling
def scaling_setup(
    m,
    state={
        "H2O": 100,
        "H_+": 5e-5,
        "OH_-": 5e-5,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-4,
        "HCO3_-": 1e-4,
    },
):
    # Set some scaling factors and look for 'bad' scaling
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / state[j], index=("Liq", j)
            )

    iscale.calculate_scaling_factors(m.fs)


# -----------------------------------------------------------------------------
# Start test class
def min_boron_removal_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # create dict to define ions (the prop pack requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
        "mw_data": {
            "H2O": 18e-3,
            "B[OH]3": 61.83e-3,
            "B[OH]4_-": 78.83e-3,
            "Na_+": 23e-3,
        },
        "charge": {
            "B[OH]4_-": -1,
            "Na_+": 1,
            "B[OH]3": 0,
        },
    }

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)

    map = {
        "boron_name": "B[OH]3",  # [is required]
        "borate_name": "B[OH]4_-",  # [is required]
        "caustic_additive": {
            "cation_name": "Na_+",  # [is optional]
            "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
            "moles_cation_per_additive": 1,  # [is required]
        },
    }
    m.fs.unit = BoronRemoval(
        property_package=m.fs.properties, chemical_mapping_data=map
    )
    state = {
        "H2O": 100,
        "H_+": 1e-7,
        "OH_-": 1e-7,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-15,
        "HCO3_-": 1e-4,
    }

    # For this test, we don't actually know how much
    #   of the base to add. Instead, we have a target
    #   exit flow of boron (see below). We will initially
    #   guess that we want 5 mg/L of additive, but the
    #   actual solution is 10 mg/L.
    model_setup(m, chem_add=0.9e-5, state=state)

    # Modified this test to fix the desired outlet flow
    #   of Boron and unfix the dosage needed to get that outlet
    m.fs.unit.outlet.flow_mol_phase_comp[(0, "Liq", "B[OH]3")].fix(1.98677e-5)
    m.fs.unit.caustic_dose_rate.unfix()

    scaling_setup(m)

    return m


class TestBoronRemoval_IonPropPack_Min(UnitTestHarness):
    def configure(self):
        m = min_boron_removal_model()

        self.default_small = 1e-8

        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]3"]
        ] = 1.98677e-5
        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]4_-"]
        ] = 1.81133e-4
        self.unit_solutions[m.fs.unit.outlet_pH()] = 10.171
        self.unit_solutions[m.fs.unit.outlet_pOH()] = 3.8257
        self.unit_solutions[m.fs.unit.caustic_dose_rate[0]] = 1.8e-5
        return m


def alk_boron_removal_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # create dict to define ions (the prop pack requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+", "HCO3_-"],
        "mw_data": {
            "H2O": 18e-3,
            "B[OH]3": 61.83e-3,
            "B[OH]4_-": 78.83e-3,
            "Na_+": 23e-3,
            "HCO3_-": 61e-3,
        },
        "charge": {"B[OH]3": 0, "B[OH]4_-": -1, "Na_+": 1, "HCO3_-": -1},
    }

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)

    map = {
        "boron_name": "B[OH]3",  # [is required]
        "borate_name": "B[OH]4_-",  # [is required]
        "caustic_additive": {
            "cation_name": "Na_+",  # [is optional]
            "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
            "moles_cation_per_additive": 1,  # [is required]
        },
    }
    m.fs.unit = BoronRemoval(
        property_package=m.fs.properties, chemical_mapping_data=map
    )

    model_setup(m, chem_add=1.8e-4)
    scaling_setup(m)

    return m


class TestBoronRemoval_IonPropPack_with_ResAlk(UnitTestHarness):
    def configure(self):
        m = alk_boron_removal_model()

        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]3"]
        ] = 1.36911e-6
        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]4_-"]
        ] = 1.99642e-4
        self.unit_solutions[m.fs.unit.outlet_pH()] = 11.375
        self.unit_solutions[m.fs.unit.outlet_pOH()] = 2.6218
        return m


def base_boron_removal_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # create dict to define ions (the prop pack requires this)
    ion_dict = {
        "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
        "mw_data": {
            "H2O": 18e-3,
            "B[OH]3": 61.83e-3,
            "B[OH]4_-": 78.83e-3,
            "Na_+": 23e-3,
        },
        "charge": {
            "B[OH]3": 0,
            "B[OH]4_-": -1,
            "Na_+": 1,
        },
    }

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)

    map = {
        "boron_name": "B[OH]3",  # [is required]
        "borate_name": "B[OH]4_-",  # [is required]
        "caustic_additive": {
            "cation_name": "Na_+",  # [is optional]
            "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
            "moles_cation_per_additive": 1,  # [is required]
        },
    }
    m.fs.unit = BoronRemoval(
        property_package=m.fs.properties, chemical_mapping_data=map
    )

    model_setup(m, chem_add=0)
    scaling_setup(m)

    return m


class TestBoronRemoval_IonPropPack_with_ResBase(UnitTestHarness):
    def configure(self):
        m = base_boron_removal_model()

        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]3"]
        ] = 1.20642e-4
        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", "B[OH]4_-"]
        ] = 8.03579e-5
        self.unit_solutions[m.fs.unit.outlet_pH()] = 9.0343
        self.unit_solutions[m.fs.unit.outlet_pOH()] = 4.9621
        return m


# -----------------------------------------------------------------------------
# Start test class with bad config
class TestBoronRemoval_BadConfigs:
    @pytest.fixture(scope="class")
    def boron_removal_bad_configs(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "H_+", "OH_-", "Na_+"],
            "mw_data": {
                "H2O": 18e-3,
                "B[OH]3": 61.83e-3,
                "B[OH]4_-": 78.83e-3,
                "H_+": 1e-3,
                "OH_-": 17e-3,
                "Na_+": 23e-3,
            },
            "charge": {
                "B[OH]4_-": -1,
            },
        }

        # attach prop pack to flowsheet
        m.fs.properties = MCASParameterBlock(**ion_dict, ignore_neutral_charge=True)

        return m

    @pytest.mark.unit
    def test_build_failures(self, boron_removal_bad_configs):
        m = boron_removal_bad_configs

        # Invalid name of boron
        map = {
            "boron_name": "Boron",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'boron_name' {Boron} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of borate
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "Borate",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'borate_name' {Borate} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Empty dict
        map = {}
        error_msg = "Did not provide a 'dict' for 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of hydroxide
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "hydroxide_name": "Raccoon_City",
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'hydroxide_name' {Raccoon_City} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of protons
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "proton_name": "Pulled_Pork",
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'proton_name' {Pulled_Pork} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of cation
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Im_Batman",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'cation_name' {Im_Batman} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Missing information (borate_name)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "caustic_additive": {
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Missing some required information in 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Missing information (mw_additive)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Missing some required information in 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Improper data type (mw_additive must be a tuple with value and units)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is required]
                "mw_additive": 40,  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Did not provide a tuple for 'mw_additive' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )
