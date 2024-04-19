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
Tests for CSTR unit model with injection.
Authors: Andrew Lee, Vibhav Dabadghao
"""

import pytest
from pyomo.environ import (
    ConcreteModel,
    units,
    value,
    assert_optimal_termination,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    initialization_tester,
)
from idaes.core.util.exceptions import ConfigurationError
from watertap.core.solvers import get_solver
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from watertap.unit_models.aeration_tank import AerationTank, ElectricityConsumption

from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting
from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)

from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM1ParameterBlock()
    m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = AerationTank(
        property_package=m.fs.properties, reaction_package=m.fs.reactions
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 16

    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions
    assert m.fs.unit.config.electricity_consumption == ElectricityConsumption.calculated
    assert m.fs.unit.config.has_aeration


class TestAeration_withASM1(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()
        m.fs.reactions = ASM1ReactionParameterBlock(property_package=m.fs.properties)

        m.fs.unit = AerationTank(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
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

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        assert hasattr(model.fs.unit, "inlet")
        assert len(model.fs.unit.inlet.vars) == 5
        assert hasattr(model.fs.unit.inlet, "flow_vol")
        assert hasattr(model.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")
        assert hasattr(model.fs.unit.inlet, "alkalinity")

        assert hasattr(model.fs.unit, "outlet")
        assert len(model.fs.unit.inlet.vars) == 5
        assert hasattr(model.fs.unit.inlet, "flow_vol")
        assert hasattr(model.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")
        assert hasattr(model.fs.unit.inlet, "alkalinity")

        assert hasattr(model.fs.unit, "cstr_performance_eqn")
        assert hasattr(model.fs.unit, "volume")
        assert hasattr(model.fs.unit, "hydraulic_retention_time")
        assert hasattr(model.fs.unit, "KLa")
        assert hasattr(model.fs.unit, "S_O_eq")

        assert hasattr(model.fs.unit, "injection")
        for k in model.fs.unit.injection:
            assert (
                model.fs.unit.injection[k]
                == model.fs.unit.control_volume.mass_transfer_term[k]
            )

        assert number_variables(model) == 102
        assert number_total_constraints(model) == 49
        assert number_unused_variables(model) == 3

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)
        assert_units_equivalent(model.fs.unit.volume[0], units.m**3)
        assert_units_equivalent(model.fs.unit.electricity_consumption[0], units.kW)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(308.15, abs=1e-2) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(2.197e-3, abs=1e-2) == value(
            model.fs.unit.outlet.conc_mass_comp[0, "S_O"]
        )
        assert pytest.approx(18.3765, abs=1e-3) == value(
            model.fs.unit.electricity_consumption[0]
        )
        assert pytest.approx(2092.2123, abs=1e-3) == value(
            model.fs.unit.hydraulic_retention_time[0]
        )
        assert pytest.approx(8.2694, abs=1e-3) == value(model.fs.unit.KLa)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_vol[0] - model.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert abs(
            value(
                model.fs.unit.inlet.flow_vol[0]
                * sum(
                    model.fs.unit.inlet.conc_mass_comp[0, j]
                    for j in model.fs.properties.solute_set
                )
                - model.fs.unit.outlet.flow_vol[0]
                * sum(
                    model.fs.unit.outlet.conc_mass_comp[0, j]
                    for j in model.fs.properties.solute_set
                )
                + sum(
                    model.fs.unit.control_volume.rate_reaction_generation[0, "Liq", j]
                    for j in model.fs.properties.solute_set
                )
                <= 1e-6
            )
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, model):
        m = model

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        solver = get_solver()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(613502.91662, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(0.0132455, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.unit
    def test_report(self, model):
        model.fs.unit.report()


@pytest.mark.build
@pytest.mark.unit
def test_with_asm2d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = ASM2dParameterBlock()
    m.fs.reactions = ASM2dReactionParameterBlock(property_package=m.fs.properties)

    m.fs.unit = AerationTank(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
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

    m.fs.unit = AerationTank(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
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
        m.fs.unit = AerationTank(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
        )
