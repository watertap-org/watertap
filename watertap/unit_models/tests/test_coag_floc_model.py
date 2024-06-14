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
from watertap.property_models.unit_specific.coagulation_prop_pack import (
    CoagulationParameterBlock,
)
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.coag_floc_model import CoagulationFlocculation
from pyomo.environ import (
    ConcreteModel,
    value,
    units as pyunits,
)
from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
from watertap.core.solvers import get_solver
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import re

__author__ = "Austin Ladshaw"

solver = get_solver()


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = CoagulationParameterBlock()
    ## NOTE: These values provided are just DUMMY values for the purposes
    #        of testing. They are not meant to be representative of any
    #        particular chemicals or real-world additives.
    chem_dict = {
        "Alum": {
            "parameter_data": {
                "mw_additive": (200, pyunits.g / pyunits.mol),
                "moles_salt_per_mole_additive": 3,
                "mw_salt": (100, pyunits.g / pyunits.mol),
            }
        },
        "Poly": {
            "parameter_data": {
                "mw_additive": (25, pyunits.g / pyunits.mol),
                "moles_salt_per_mole_additive": 0,
                "mw_salt": (23, pyunits.g / pyunits.mol),
            }
        },
    }
    m.fs.unit = CoagulationFlocculation(
        property_package=m.fs.properties, chemical_additives=chem_dict
    )

    m.fs.unit.fix_tss_turbidity_relation_defaults()
    m.fs.unit.initial_turbidity_ntu.fix()
    m.fs.unit.final_turbidity_ntu.fix(5)
    m.fs.unit.chemical_doses[0, "Alum"].fix(10)
    m.fs.unit.chemical_doses[0, "Poly"].fix(5)

    # set the inlet streams
    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TSS"].fix(0.01)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Sludge"].fix(0.0)

    # set performance vars
    m.fs.unit.rapid_mixing_retention_time[0].fix(60)
    m.fs.unit.num_rapid_mixing_basins.fix(4)
    m.fs.unit.rapid_mixing_vel_grad[0].fix(750)

    m.fs.unit.floc_retention_time[0].fix(1800)
    m.fs.unit.single_paddle_length.fix(4)
    m.fs.unit.single_paddle_width.fix(0.5)
    m.fs.unit.paddle_rotational_speed[0].fix(0.03)

    m.fs.unit.paddle_drag_coef[0].fix(1.5)
    m.fs.unit.vel_fraction.fix(0.7)
    m.fs.unit.num_paddle_wheels.fix(4)
    m.fs.unit.num_paddles_per_wheel.fix(4)

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_in[0.0].mass_frac_phase_comp[
            "Liq", "Sludge"
        ],
        1e12,
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_no_chemicals():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = CoagulationParameterBlock()
    m.fs.unit = CoagulationFlocculation(property_package=m.fs.properties)

    m.fs.unit.fix_tss_turbidity_relation_defaults()
    m.fs.unit.initial_turbidity_ntu.fix(5000)
    m.fs.unit.final_turbidity_ntu.fix(100)

    tss_in = value(m.fs.unit.compute_inlet_tss_mass_flow(0))

    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TSS"].fix(tss_in)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Sludge"].fix(0.0)

    # set performance vars
    m.fs.unit.rapid_mixing_retention_time[0].fix(60)
    m.fs.unit.num_rapid_mixing_basins.fix(4)
    m.fs.unit.rapid_mixing_vel_grad[0].fix(750)

    m.fs.unit.floc_retention_time[0].fix(1800)
    m.fs.unit.single_paddle_length.fix(4)
    m.fs.unit.single_paddle_width.fix(0.5)
    m.fs.unit.paddle_rotational_speed[0].fix(0.03)

    m.fs.unit.paddle_drag_coef[0].fix(1.5)
    m.fs.unit.vel_fraction.fix(0.7)
    m.fs.unit.num_paddle_wheels.fix(4)
    m.fs.unit.num_paddles_per_wheel.fix(4)

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_in[0.0].mass_frac_phase_comp[
            "Liq", "Sludge"
        ],
        1e12,
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestCoagFloc(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]] = 1
        self.unit_solutions[
            m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "Sludge"]
        ] = 0.009990636
        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TDS"]] = (
            0.01001510
        )
        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TSS"]] = (
            9.36352627e-6
        )

        return m


class TestCoagFlocNoChemicals(UnitTestHarness):
    def configure(self):
        m = build_no_chemicals()

        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]] = 1
        self.unit_solutions[
            m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "Sludge"]
        ] = 0.009112749
        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TDS"]] = (
            0.01
        )
        self.unit_solutions[m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TSS"]] = (
            0.0001872506
        )

        return m


class TestCoagFlocErrorLog:
    @pytest.mark.unit
    def test_dictionary_error(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = CoagulationParameterBlock()

        bad_dict1 = {
            "Alum": {
                "foo_bar": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }

        with pytest.raises(
            ConfigurationError,
            match="Did not provide a 'parameter_data' for chemical",
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict1
            )

        bad_dict2 = {
            "Alum": {
                "parameter_data": {
                    "foo_bar": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a 'mw_additive' for chemical"
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict2
            )

        bad_dict3 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": "foo-bar",
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError,
            match="Did not provide a number for 'moles_salt_per_mole_additive'",
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict3
            )

        bad_dict4 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "foo-bar": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError,
            match="Did not provide a 'moles_salt_per_mole_additive' for chemical",
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict4
            )

        bad_dict5 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "foo-bar": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a 'mw_salt' for chemical"
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict5
            )

        bad_dict6 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": "not a tuple",
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a tuple for 'mw_additive'"
        ):
            m.fs.unit = CoagulationFlocculation(
                property_package=m.fs.properties, chemical_additives=bad_dict6
            )

    @pytest.mark.unit
    def test_property_error(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        error_msg = (
            "Coagulation-Flocculation model MUST contain ('Liq','TDS') "
            "as a component, but the property package has only specified "
            "the following components [('Liq', 'H2O'), ('Liq', 'NaCl')]"
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.properties = NaClParameterBlock()
            m.fs.unit = CoagulationFlocculation(property_package=m.fs.properties)

        error_msg = (
            "Coagulation-Flocculation model MUST contain ('Liq','Sludge') "
            "as a component, but the property package has only specified "
            "the following components [('Liq', 'H2O'), ('Liq', 'TDS')]"
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.properties = SeawaterParameterBlock()
            m.fs.unit = CoagulationFlocculation(property_package=m.fs.properties)
