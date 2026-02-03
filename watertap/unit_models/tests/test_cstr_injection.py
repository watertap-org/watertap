#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for CSTR unit model with injection.
Authors: Andrew Lee, Adam Atia, Vibhav Dabadghao
"""

from io import StringIO
import platform
import pytest
from pyomo.environ import (
    ConcreteModel,
    units,
    value,
    Objective,
    Suffix,
    TransformationFactory,
)
from idaes.core import (
    FlowsheetBlock,
)
import re
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)
from idaes.core.scaling.scaling_base import ScalerBase
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from idaes.core.util.exceptions import ConfigurationError
from watertap.core.solvers import get_solver

from watertap.unit_models.cstr_injection import (
    CSTR_Injection,
    ElectricityConsumption,
    CSTR_InjectionScaler,
)
from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting
from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
    ASM1PropertiesScaler,
)
from watertap.property_models.unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
    ASM1ReactionScaler,
)
from watertap.property_models.unit_specific.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)

from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

is_reference_platform = (
    platform.system() == "Windows"
    and platform.python_version_tuple()[0] == "3"
    and platform.python_version_tuple()[1] == "11"
)

reference_platform_only = pytest.mark.xfail(
    condition=(not is_reference_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the reference platform (Python 3.11 on Windows)",
)


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_equilibrium_reactions=False,
        has_heat_transfer=True,
        has_heat_of_reaction=True,
        has_pressure_change=True,
    )

    m.fs.unit.inlet.flow_vol.fix(1.0e-03)
    m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.unit.inlet.temperature.fix(303.15)
    m.fs.unit.inlet.pressure.fix(101325.0)

    m.fs.unit.volume.fix(1.5e-03)
    m.fs.unit.heat_duty.fix(0)
    m.fs.unit.deltaP.fix(0)
    m.fs.unit.injection.fix(0)

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].conc_mol_comp["H2O"], 1e-4
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].flow_vol, 1e3
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


def build_ASM1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM1ParameterBlock()
    m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )

    m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
    m.fs.unit.inlet.temperature.fix(308.15 * units.K)
    m.fs.unit.inlet.pressure.fix(1 * units.atm)
    m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
    m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

    m.fs.unit.volume.fix(500)
    m.fs.unit.injection.fix(0)
    m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].conc_mass_comp["X_P"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_S"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_S"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_BH"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_P"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_O"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_NH"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ND"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_ND"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ALK"], 1e3
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R4"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R6"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R7"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R8"], 1e5
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R1"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R4"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R6"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R7"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R8"], 1e5
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestCSTR_injection(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[m.fs.unit.outlet.pressure[0]] = 101325.0
        self.unit_solutions[m.fs.unit.outlet.temperature[0]] = 304.09
        self.unit_solutions[m.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]] = 20.32
        self.unit_solutions[m.fs.unit.electricity_consumption[0]] = 0.03888
        self.unit_solutions[m.fs.unit.hydraulic_retention_time[0]] = 1.5
        self.unit_solutions[m.fs.unit.control_volume.heat_of_reaction[0]] = 3904.51

        # Conservation checks
        self.unit_solutions[m.fs.unit.outlet.flow_vol[0]] = value(
            m.fs.unit.inlet.flow_vol[0]
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0]
                * sum(
                    m.fs.unit.inlet.conc_mol_comp[0, j]
                    for j in m.fs.properties.component_list
                ),
                "out": m.fs.unit.outlet.flow_vol[0]
                * sum(
                    m.fs.unit.inlet.conc_mol_comp[0, j]
                    for j in m.fs.properties.component_list
                ),
            },
            "Check 2": {
                "in": (
                    m.fs.unit.inlet.flow_vol[0]
                    * m.fs.properties.dens_mol
                    * m.fs.properties.cp_mol
                    * (m.fs.unit.inlet.temperature[0] - m.fs.properties.temperature_ref)
                ),
                "out": m.fs.unit.outlet.flow_vol[0]
                * m.fs.properties.dens_mol
                * m.fs.properties.cp_mol
                * (m.fs.unit.outlet.temperature[0] - m.fs.properties.temperature_ref)
                - m.fs.unit.control_volume.heat_of_reaction[0],
            },
        }

        return m


class TestCosting_Saponification(UnitTestHarness):
    def configure(self):
        m = build()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        m.objective = Objective(expr=m.fs.costing.LCOW)

        iscale.calculate_scaling_factors(m.fs.unit)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-2)

        m.fs.unit.initialize()

        solver.solve(m)

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 7.75429 * 2
        self.unit_solutions[m.fs.costing.LCOW] = 0.00082698

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0],
                "out": m.fs.unit.outlet.flow_vol[0],
            },
        }

        return m


class TestSaponification(object):
    @pytest.mark.unit
    def test_get_performance_contents(self):

        sapon = build()
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Volume": sapon.fs.unit.volume[0],
                "Injection [('Liq', 'H2O')]": sapon.fs.unit.injection[0, "Liq", "H2O"],
                "Injection [('Liq', 'NaOH')]": sapon.fs.unit.injection[
                    0, "Liq", "NaOH"
                ],
                "Injection [('Liq', 'EthylAcetate')]": sapon.fs.unit.injection[
                    0, "Liq", "EthylAcetate"
                ],
                "Injection [('Liq', 'SodiumAcetate')]": sapon.fs.unit.injection[
                    0, "Liq", "SodiumAcetate"
                ],
                "Injection [('Liq', 'Ethanol')]": sapon.fs.unit.injection[
                    0, "Liq", "Ethanol"
                ],
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
                "Electricity Consumption": sapon.fs.unit.electricity_consumption[0],
            }
        }


class TestCSTR_injection_ASM1(UnitTestHarness):
    def configure(self):
        m = build_ASM1()

        iscale.set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0].conc_mass_comp["X_P"], 1e5
        )

        self.unit_solutions[m.fs.unit.outlet.pressure[0]] = 101325.0
        self.unit_solutions[m.fs.unit.outlet.temperature[0]] = 308.15
        self.unit_solutions[m.fs.unit.outlet.conc_mass_comp[0, "S_O"]] = 6.258e-3
        self.unit_solutions[m.fs.unit.electricity_consumption[0]] = 18.3765
        self.unit_solutions[m.fs.unit.hydraulic_retention_time[0]] = 2092.2123
        self.unit_solutions[m.fs.unit.KLa] = 8.2694

        # Conservation checks
        self.unit_solutions[m.fs.unit.outlet.flow_vol[0]] = value(
            m.fs.unit.inlet.flow_vol[0]
        )

        self.default_relative_tolerance = 1e-2

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0]
                * sum(
                    m.fs.unit.inlet.conc_mass_comp[0, j]
                    for j in m.fs.properties.solute_set
                ),
                "out": m.fs.unit.outlet.flow_vol[0]
                * sum(
                    m.fs.unit.outlet.conc_mass_comp[0, j]
                    for j in m.fs.properties.solute_set
                )
                + sum(
                    m.fs.unit.control_volume.rate_reaction_generation[0, "Liq", j]
                    for j in m.fs.properties.solute_set
                ),
            },
        }

        return m


class TestCosting(UnitTestHarness):
    def configure(self):
        m = build_ASM1()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        m.objective = Objective(expr=m.fs.costing.LCOW)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-7)
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0].conc_mass_comp["X_P"], 1e5
        )

        iscale.calculate_scaling_factors(m.fs.unit)

        m.fs.unit.initialize()

        solver.solve(m)

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 613502.91662
        self.unit_solutions[m.fs.costing.LCOW] = 0.0132455

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0],
                "out": m.fs.unit.outlet.flow_vol[0],
            },
        }

        return m


@pytest.mark.build
@pytest.mark.unit
def test_with_asm2d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM2dParameterBlock()
    m.fs.reactions = ASM2dReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )


@pytest.mark.build
@pytest.mark.unit
def test_with_mod_asm2d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ModifiedASM2dParameterBlock()
    m.fs.reactions = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )


@pytest.mark.unit
def test_error_without_oxygen():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ADM1ParameterBlock()
    m.fs.reactions = ADM1ReactionParameterBlock(property_package=m.fs.properties)

    # Expect exception if has_aeration=True but S_O or S_O2 not listed in component_list of prop package.
    with pytest.raises(
        ConfigurationError,
        match="has_aeration was set to True, but the property package has neither 'S_O' nor 'S_O2' in its list of components.",
    ):
        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
        )


class TestCSTR_InjectionScaler:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()
        m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
            electricity_consumption=ElectricityConsumption.calculated,
        )

        m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
        m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

        m.fs.unit.volume.fix(500)
        m.fs.unit.injection.fix(0)
        m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

        return m

    @pytest.mark.component
    def test_variable_scaling_routine(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, CSTR_InjectionScaler)

        scaler.variable_scaling_routine(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.control_volume.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        # Scaling factors for FTP
        assert len(sfx_in) == 3

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.control_volume.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        # Scaling factors for FTP
        assert len(sfx_out) == 3

        # Reaction block
        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for rxn rate
        assert len(sfx_rxn) == 8

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.control_volume.scaling_factor
        assert isinstance(sfx_cv, Suffix)
        # Scaling factors for volume and oxygen mass transfer
        assert len(sfx_cv) == 2

    @pytest.mark.component
    def test_constraint_scaling_routine(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, CSTR_InjectionScaler)

        scaler.constraint_scaling_routine(model.fs.unit)

        assert not hasattr(
            model.fs.unit.control_volume.properties_out[0], "scaling_factor"
        )

        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for rate expressions
        assert len(sfx_rxn) == 8

        sfx_unit = model.fs.unit.scaling_factor
        assert isinstance(sfx_unit, Suffix)
        # Scaling factors for retention time, mass transfer, and other unit model constraints
        assert len(sfx_unit) == 11

    @pytest.mark.component
    def test_scale_model(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, CSTR_InjectionScaler)

        scaler.scale_model(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.control_volume.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        # Scaling factors for FTP
        assert len(sfx_in) == 3

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.control_volume.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        # Scaling factors for FTP
        assert len(sfx_out) == 3

        # Reaction block
        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for rxn rates and rate expressions
        assert len(sfx_rxn) == 16

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.control_volume.scaling_factor
        assert isinstance(sfx_cv, Suffix)
        # Scaling factors for volume, mass transfer, and other control volume variables/constraints
        assert len(sfx_cv) == 32

        sfx_unit = model.fs.unit.scaling_factor
        assert isinstance(sfx_unit, Suffix)
        # Scaling factors for HRT, mass transfer, electricity consumption, and other variables/constraints
        assert len(sfx_unit) == 13

    @pytest.mark.integration
    def test_example_case_iscale(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()
        m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
            electricity_consumption=ElectricityConsumption.calculated,
        )

        m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
        m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

        m.fs.unit.volume.fix(500)
        m.fs.unit.injection.fix(0)
        m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

        # Set scaling factors for badly scaled variables
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].conc_mass_comp["X_P"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_S"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_S"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_BH"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_P"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_O"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_NH"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ND"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_ND"], 1e3
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ALK"], 1e3
        )

        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_extent[0.0, "R4"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_extent[0.0, "R6"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_extent[0.0, "R7"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.rate_reaction_extent[0.0, "R8"], 1e5
        )

        iscale.set_scaling_factor(
            m.fs.unit.control_volume.reactions[0.0].reaction_rate["R1"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.reactions[0.0].reaction_rate["R4"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.reactions[0.0].reaction_rate["R6"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.reactions[0.0].reaction_rate["R7"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.control_volume.reactions[0.0].reaction_rate["R8"], 1e5
        )

        iscale.calculate_scaling_factors(m.fs.unit)

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.0843930927394e13, rel=1e-3
        )

    @pytest.mark.integration
    def test_example_case_scaler(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()
        m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

        m.fs.unit = CSTR_Injection(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_aeration=True,
            electricity_consumption=ElectricityConsumption.calculated,
        )

        m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
        m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

        m.fs.unit.volume.fix(500)
        m.fs.unit.injection.fix(0)
        m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

        # Set scaling factors for badly scaled variables
        sb = ScalerBase()
        sb.set_variable_scaling_factor(m.fs.unit.hydraulic_retention_time[0], 1e-3)

        scaler = CSTR_InjectionScaler()
        scaler.scale_model(
            m.fs.unit,
            submodel_scalers={
                m.fs.unit.control_volume.properties_in: ASM1PropertiesScaler,
                m.fs.unit.control_volume.properties_out: ASM1PropertiesScaler,
                m.fs.unit.control_volume.reactions: ASM1ReactionScaler,
            },
        )

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.1526931e7, rel=1e-3
        )


def build_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM1ParameterBlock()
    m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = CSTR_Injection(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_aeration=True,
        electricity_consumption=ElectricityConsumption.calculated,
    )

    m.fs.unit.inlet.flow_vol.fix(20648 * units.m**3 / units.day)
    m.fs.unit.inlet.temperature.fix(308.15 * units.K)
    m.fs.unit.inlet.pressure.fix(1 * units.atm)
    m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)
    m.fs.unit.inlet.alkalinity.fix(7 * units.mol / units.m**3)

    m.fs.unit.volume.fix(500)
    m.fs.unit.injection.fix(0)
    m.fs.unit.injection[0, "Liq", "S_O"].fix(2e-3)

    solver = get_solver()
    solver.solve(m)

    return m


def scale_vars_with_scalers(m):
    scaler = CSTR_InjectionScaler()
    scaler.scale_model(
        m.fs.unit,
        submodel_scalers={
            m.fs.unit.control_volume.properties_in: ASM1PropertiesScaler,
            m.fs.unit.control_volume.properties_out: ASM1PropertiesScaler,
            m.fs.unit.control_volume.reactions: ASM1ReactionScaler,
        },
    )


def scale_vars_with_iscale(m):
    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.properties_out[0.0].conc_mass_comp["X_P"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_S"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_S"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_BH"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_P"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_O"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_NH"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ND"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "X_ND"], 1e3
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_generation[0.0, "Liq", "S_ALK"], 1e3
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R4"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R6"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R7"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R8"], 1e5
    )

    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R1"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R4"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R6"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R7"], 1e5
    )
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.reactions[0.0].reaction_rate["R8"], 1e5
    )

    iscale.calculate_scaling_factors(m.fs.unit)


def perturb_solution(m):
    m.fs.unit.inlet.flow_vol.fix(20648 * 0.9 * units.m**3 / units.day)
    m.fs.unit.volume.fix(500 * 0.85)


Scaling Profile Report
----------------------------------------------------------------------------
Scaling Method           || User Scaling           || Perfect Scaling
Unscaled                 || 1.826E+16 | Solved 4   ||
Vars Only                || 2.740E+17 | Solved 4   || 2.014E+21 | Solved 4  
Harmonic                 || 2.740E+17 | Solved 4   || 4.443E+22 | Solved 18 
Inverse Sum              || 2.740E+17 | Solved 4   || 2.399E+14 | Solved 4  
Inverse Root Sum Squares || 2.740E+17 | Solved 4   || 3.412E+14 | Solved 4  
Inverse Maximum          || 2.740E+17 | Solved 4   || 4.809E+14 | Solved 4  
Inverse Minimum          || 2.740E+17 | Solved 4   || 4.455E+22 | Solved 18 
Nominal L1 Norm          || 2.740E+17 | Solved 4   || 2.841E+14 | Solved 4  
Nominal L2 Norm          || 2.740E+17 | Solved 4   || 3.755E+14 | Solved 4  
Actual L1 Norm           || 2.740E+17 | Solved 4   || 5.461E+13 | Solved 4  
Actual L2 Norm           || 2.740E+17 | Solved 4   || 6.491E+13 | Solved 4  
"""

    assert assert_solve_tolerance(stream.getvalue(), expected, tolerance=3)


@pytest.mark.requires_idaes_solver
@pytest.mark.unit
# @reference_platform_only
def test_scaling_profiler_with_iscale():
    sp = ScalingProfiler(
        build_model=build_model,
        user_scaling=scale_vars_with_iscale,
        perturb_state=perturb_solution,
    )

    stream = StringIO()

    sp.report_scaling_profiles(stream=stream)

    expected = """
Scaling Profile Report
----------------------------------------------------------------------------
Scaling Method           || User Scaling           || Perfect Scaling
Unscaled                 || 1.826E+16 | Solved 4   ||
Vars Only                || 8.948E+12 | Solved 4   || 2.014E+21 | Solved 4  
Harmonic                 || 1.044E+17 | Solved 59  || 4.443E+22 | Solved 18 
Inverse Sum              || 5.247E+17 | Failed 50  || 2.399E+14 | Solved 4  
Inverse Root Sum Squares || 5.220E+17 | Failed 55  || 3.412E+14 | Solved 4  
Inverse Maximum          || 5.208E+17 | Failed 52  || 4.809E+14 | Solved 4  
Inverse Minimum          || 2.103E+17 | Solved 65  || 4.455E+22 | Solved 18 
Nominal L1 Norm          || 7.817E+09 | Solved 4   || 2.841E+14 | Solved 4  
Nominal L2 Norm          || 1.278E+10 | Solved 4   || 3.755E+14 | Solved 4  
Actual L1 Norm           || 3.950E+09 | Solved 3   || 5.461E+13 | Solved 4  
Actual L2 Norm           || 4.339E+09 | Solved 3   || 6.491E+13 | Solved 4  
"""

    assert assert_solve_tolerance(stream.getvalue(), expected, tolerance=3)


def assert_solve_tolerance(output, expected, tolerance=0):
    output_lines = output.strip().splitlines()
    expected_lines = expected.strip().splitlines()

    assert len(output_lines) == len(expected_lines)

    for output_line, expected_line in zip(output_lines, expected_lines):
        # find the groups that match 'Solved <number>'
        output_solved = re.search(r"Solved (\d+)", output_line)
        expected_solved = re.search(r"Solved (\d+)", expected_line)

        if output_solved and expected_solved:
            # get the number after 'Solved'
            output_num = int(output_solved.group(1))
            expected_num = int(expected_solved.group(1))

            # assert that the numbers after 'Solved' are within the tolerance
            assert abs(output_num - expected_num) <= tolerance, (
                f"Solved number mismatch: output {output_num} "
                f"vs expected {expected_num} within tolerance {tolerance}"
            )

        # checks if the non-numeric parts of the statemnt and the actual values are equal
        assert re.sub(r"Solved \d+", "Solved X", output_line) == re.sub(
            r"Solved \d+", "Solved X", expected_line
        ), f"Line mismatch: '{output_line}' vs '{expected_line}'"
# TODO Replace these scaling profiler tests with more detailed convergence analysis
# @pytest.mark.requires_idaes_solver
# @pytest.mark.unit
# def test_scaling_profiler_with_scalers():
#     sp = ScalingProfiler(
#         build_model=build_model,
#         user_scaling=scale_vars_with_scalers,
#         perturb_state=perturb_solution,
#     )

#     stream = StringIO()

#     sp.report_scaling_profiles(stream=stream)

#     expected = """
# ============================================================================
# Scaling Profile Report
# ----------------------------------------------------------------------------
# Scaling Method           || User Scaling           || Perfect Scaling
# Unscaled                 || 1.826E+16 | Solved 4   ||
# Vars Only                || 2.740E+17 | Solved 4   || 2.014E+21 | Solved 4
# Harmonic                 || 2.740E+17 | Solved 4   || 4.443E+22 | Solved 18
# Inverse Sum              || 2.740E+17 | Solved 4   || 2.399E+14 | Solved 4
# Inverse Root Sum Squares || 2.740E+17 | Solved 4   || 3.412E+14 | Solved 4
# Inverse Maximum          || 2.740E+17 | Solved 4   || 4.809E+14 | Solved 4
# Inverse Minimum          || 2.740E+17 | Solved 4   || 4.455E+22 | Solved 18
# Nominal L1 Norm          || 2.740E+17 | Solved 4   || 2.841E+14 | Solved 4
# Nominal L2 Norm          || 2.740E+17 | Solved 4   || 3.755E+14 | Solved 4
# Actual L1 Norm           || 2.740E+17 | Solved 4   || 5.461E+13 | Solved 4
# Actual L2 Norm           || 2.740E+17 | Solved 4   || 6.491E+13 | Solved 4
# ============================================================================
# """

#     assert stream.getvalue() == expected


# @pytest.mark.requires_idaes_solver
# @pytest.mark.unit
# def test_scaling_profiler_with_iscale():
#     sp = ScalingProfiler(
#         build_model=build_model,
#         user_scaling=scale_vars_with_iscale,
#         perturb_state=perturb_solution,
#     )

#     stream = StringIO()

#     sp.report_scaling_profiles(stream=stream)

#     expected = """
# ============================================================================
# Scaling Profile Report
# ----------------------------------------------------------------------------
# Scaling Method           || User Scaling           || Perfect Scaling
# Unscaled                 || 1.826E+16 | Solved 4   ||
# Vars Only                || 8.948E+12 | Solved 4   || 2.014E+21 | Solved 4
# Harmonic                 || 1.044E+17 | Solved 57  || 4.443E+22 | Solved 18
# Inverse Sum              || 5.247E+17 | Failed 50  || 2.399E+14 | Solved 4
# Inverse Root Sum Squares || 5.220E+17 | Failed 55  || 3.412E+14 | Solved 4
# Inverse Maximum          || 5.208E+17 | Failed 52  || 4.809E+14 | Solved 4
# Inverse Minimum          || 2.103E+17 | Solved 65  || 4.455E+22 | Solved 18
# Nominal L1 Norm          || 7.817E+09 | Solved 4   || 2.841E+14 | Solved 4
# Nominal L2 Norm          || 1.278E+10 | Solved 4   || 3.755E+14 | Solved 4
# Actual L1 Norm           || 3.950E+09 | Solved 3   || 5.461E+13 | Solved 4
# Actual L2 Norm           || 4.339E+09 | Solved 3   || 6.491E+13 | Solved 4
# ============================================================================
# """

#     assert stream.getvalue() == expected
