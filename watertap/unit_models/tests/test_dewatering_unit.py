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
Tests for dewatering unit example.
"""

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
)

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)

from idaes.models.unit_models.separator import SplittingType

from pyomo.environ import (
    units,
)

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import (
    initialization_tester,
)

from idaes.core.util.exceptions import (
    ConfigurationError,
)

from watertap.unit_models.dewatering import DewateringUnit, ActivatedSludgeModelType
from watertap.property_models.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)

from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from pyomo.util.check_units import assert_units_consistent
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.dewatering import (
    cost_dewatering,
    cost_centrifuge,
    cost_filter_belt_press,
    cost_filter_plate_press,
    DewateringType
)

__author__ = "Alejandro Garciadiego, Adam Atia"

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM1ParameterBlock()

    m.fs.unit = DewateringUnit(property_package=m.fs.props)

    assert len(m.fs.unit.config) == 16

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.activated_sludge_model == ActivatedSludgeModelType.ASM1

    assert "underflow" in m.fs.unit.config.outlet_list
    assert "overflow" in m.fs.unit.config.outlet_list
    assert SplittingType.componentFlow is m.fs.unit.config.split_basis


@pytest.mark.unit
def test_list_error():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM1ParameterBlock()

    with pytest.raises(
        ConfigurationError,
        match="fs.unit encountered unrecognised "
        "outlet_list. This should not "
        "occur - please use overflow "
        "and underflow as outlets.",
    ):
        m.fs.unit = DewateringUnit(
            property_package=m.fs.props,
            outlet_list=["outlet1", "outlet2"],
        )


# -----------------------------------------------------------------------------
class TestDu(object):
    @pytest.fixture(scope="class")
    def du(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM1ParameterBlock()

        m.fs.unit = DewateringUnit(property_package=m.fs.props)

        m.fs.unit.inlet.flow_vol.fix(178.4674 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)

        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(130.867 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(258.5789 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(
            17216.2434 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(2611.4843 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(626.0652 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(
            1442.7882 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(0.54323 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(100.8668 * units.mg / units.liter)
        m.fs.unit.inlet.alkalinity.fix(97.8459 * units.mol / units.m**3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, du):

        assert hasattr(du.fs.unit, "inlet")
        assert len(du.fs.unit.inlet.vars) == 5
        assert hasattr(du.fs.unit.inlet, "flow_vol")
        assert hasattr(du.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(du.fs.unit.inlet, "temperature")
        assert hasattr(du.fs.unit.inlet, "pressure")
        assert hasattr(du.fs.unit.inlet, "alkalinity")

        assert hasattr(du.fs.unit, "underflow")
        assert len(du.fs.unit.underflow.vars) == 5
        assert hasattr(du.fs.unit.underflow, "flow_vol")
        assert hasattr(du.fs.unit.underflow, "conc_mass_comp")
        assert hasattr(du.fs.unit.underflow, "temperature")
        assert hasattr(du.fs.unit.underflow, "pressure")
        assert hasattr(du.fs.unit.underflow, "alkalinity")

        assert hasattr(du.fs.unit, "overflow")
        assert len(du.fs.unit.overflow.vars) == 5
        assert hasattr(du.fs.unit.overflow, "flow_vol")
        assert hasattr(du.fs.unit.overflow, "conc_mass_comp")
        assert hasattr(du.fs.unit.overflow, "temperature")
        assert hasattr(du.fs.unit.overflow, "pressure")
        assert hasattr(du.fs.unit.overflow, "alkalinity")

        assert number_variables(du) == 76
        assert number_total_constraints(du) == 60
        assert number_unused_variables(du) == 0

    @pytest.mark.unit
    def test_dof(self, du):
        assert degrees_of_freedom(du) == 0

    @pytest.mark.unit
    def test_units(self, du):
        assert_units_consistent(du)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, du):

        iscale.calculate_scaling_factors(du)
        initialization_tester(du)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, du):
        solver = get_solver()
        results = solver.solve(du)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, du):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            du.fs.unit.overflow.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            du.fs.unit.overflow.temperature[0]
        )
        assert pytest.approx(0.001954, rel=1e-3) == value(
            du.fs.unit.overflow.flow_vol[0]
        )
        assert pytest.approx(0.1308, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.2585, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "S_S"]
        )
        assert pytest.approx(0.3638, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0.0552, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "X_S"]
        )
        assert value(du.fs.unit.overflow.conc_mass_comp[0, "X_BH"]) <= 1e-6
        assert value(du.fs.unit.overflow.conc_mass_comp[0, "X_BA"]) <= 1e-6
        assert pytest.approx(0.01323, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "X_P"]
        )
        assert value(du.fs.unit.overflow.conc_mass_comp[0, "S_O"]) <= 1e-6
        assert value(du.fs.unit.overflow.conc_mass_comp[0, "S_NO"]) <= 1e-6
        assert pytest.approx(1.4427, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "S_NH"]
        )
        assert pytest.approx(0.000543, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "S_ND"]
        )
        assert pytest.approx(0.00213, rel=1e-3) == value(
            du.fs.unit.overflow.conc_mass_comp[0, "X_ND"]
        )
        assert pytest.approx(0.09784, rel=1e-3) == value(
            du.fs.unit.overflow.alkalinity[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, du):
        assert (
            abs(
                value(
                    du.fs.unit.inlet.flow_vol[0] * du.fs.props.dens_mass
                    - du.fs.unit.overflow.flow_vol[0] * du.fs.props.dens_mass
                    - du.fs.unit.underflow.flow_vol[0] * du.fs.props.dens_mass
                )
            )
            <= 1e-6
        )
        for i in du.fs.props.solute_set:
            assert (
                abs(
                    value(
                        du.fs.unit.inlet.flow_vol[0]
                        * du.fs.unit.inlet.conc_mass_comp[0, i]
                        - du.fs.unit.overflow.flow_vol[0]
                        * du.fs.unit.overflow.conc_mass_comp[0, i]
                        - du.fs.unit.underflow.flow_vol[0]
                        * du.fs.unit.underflow.conc_mass_comp[0, i]
                    )
                )
                <= 1e-6
            )

    @pytest.mark.unit
    def test_report(self, du):
        du.fs.unit.report()


class TestDUASM2d(object):
    @pytest.fixture(scope="class")
    def du_asm2d(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()

        m.fs.unit = DewateringUnit(
            property_package=m.fs.props,
            activated_sludge_model=ActivatedSludgeModelType.ASM2D,
        )

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        m.fs.unit.inlet.flow_vol.fix(300 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O2"].fix(7.9707 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_N2"].fix(29.0603 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH4"].fix(8.0209 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO3"].fix(6.6395 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_PO4"].fix(7.8953 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_F"].fix(0.4748 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_A"].fix(0.0336 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(1695.7695 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(68.2975 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_H"].fix(1855.5067 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(
            214.5319 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(63.5316 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(2.7381 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(
            118.3582 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_MeOH"].fix(EPS * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_MeP"].fix(EPS * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_TSS"].fix(
            3525.429 * units.mg / units.liter
        )
        m.fs.unit.inlet.alkalinity[0].fix(4.6663 * units.mmol / units.liter)

        return m

    @pytest.mark.unit
    def test_dof(self, du_asm2d):
        assert degrees_of_freedom(du_asm2d) == 0

    @pytest.mark.unit
    def test_units(self, du_asm2d):
        assert_units_consistent(du_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, du_asm2d):

        iscale.calculate_scaling_factors(du_asm2d)
        initialization_tester(du_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, du_asm2d):
        solver = get_solver()
        results = solver.solve(du_asm2d)
        assert_optimal_termination(results)


class TestDUModifiedASM2d(object):
    @pytest.fixture(scope="class")
    def du_mod_asm2d(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()

        m.fs.unit = DewateringUnit(
            property_package=m.fs.props,
            activated_sludge_model=ActivatedSludgeModelType.modified_ASM2D,
        )

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        m.fs.unit.inlet.flow_vol.fix(300 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O2"].fix(7.9707 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_N2"].fix(29.0603 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH4"].fix(8.0209 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO3"].fix(6.6395 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_PO4"].fix(7.8953 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_F"].fix(0.4748 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_A"].fix(0.0336 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_K"].fix(7 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_Mg"].fix(6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(10 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(1695.7695 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(68.2975 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_H"].fix(1855.5067 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(
            214.5319 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(63.5316 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(2.7381 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(
            118.3582 * units.mg / units.liter
        )

        return m

    @pytest.mark.unit
    def test_dof(self, du_mod_asm2d):
        assert degrees_of_freedom(du_mod_asm2d) == 0

    @pytest.mark.unit
    def test_units(self, du_mod_asm2d):
        assert_units_consistent(du_mod_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, du_mod_asm2d):

        iscale.calculate_scaling_factors(du_mod_asm2d)
        initialization_tester(du_mod_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, du_mod_asm2d):
        solver = get_solver()
        results = solver.solve(du_mod_asm2d)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_du_costing(self, du):
        m = du

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            # costing_method=cost_dewatering,
        )
        m.fs.costing.cost_process()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1390570, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
