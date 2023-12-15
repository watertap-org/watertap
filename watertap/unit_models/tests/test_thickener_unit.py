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
Tests for thickener unit example.
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
from watertap.unit_models.thickener import Thickener, ActivatedSludgeModelType
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
from idaes.core import UnitModelCostingBlock

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

    m.fs.unit = Thickener(property_package=m.fs.props)

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
        m.fs.unit = Thickener(
            property_package=m.fs.props, outlet_list=["outlet1", "outlet2"]
        )


# -----------------------------------------------------------------------------
class TestThickASM1(object):
    @pytest.fixture(scope="class")
    def tu(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM1ParameterBlock()

        m.fs.unit = Thickener(property_package=m.fs.props)

        m.fs.unit.inlet.flow_vol.fix(300 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)

        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(28.0643 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(0.67336 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(3036.2175 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(63.2392 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(
            4442.8377 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(332.5958 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(1922.8108 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(1.3748 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(9.1948 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(0.15845 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(0.55943 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(4.7411 * units.mg / units.liter)
        m.fs.unit.inlet.alkalinity.fix(4.5646 * units.mol / units.m**3)

        m.fs.unit.hydraulic_retention_time.fix()
        m.fs.unit.diameter.fix()

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, tu):

        assert hasattr(tu.fs.unit, "inlet")
        assert len(tu.fs.unit.inlet.vars) == 5
        assert hasattr(tu.fs.unit.inlet, "flow_vol")
        assert hasattr(tu.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(tu.fs.unit.inlet, "temperature")
        assert hasattr(tu.fs.unit.inlet, "pressure")
        assert hasattr(tu.fs.unit.inlet, "alkalinity")

        assert hasattr(tu.fs.unit, "underflow")
        assert len(tu.fs.unit.underflow.vars) == 5
        assert hasattr(tu.fs.unit.underflow, "flow_vol")
        assert hasattr(tu.fs.unit.underflow, "conc_mass_comp")
        assert hasattr(tu.fs.unit.underflow, "temperature")
        assert hasattr(tu.fs.unit.underflow, "pressure")
        assert hasattr(tu.fs.unit.underflow, "alkalinity")

        assert hasattr(tu.fs.unit, "overflow")
        assert len(tu.fs.unit.overflow.vars) == 5
        assert hasattr(tu.fs.unit.overflow, "flow_vol")
        assert hasattr(tu.fs.unit.overflow, "conc_mass_comp")
        assert hasattr(tu.fs.unit.overflow, "temperature")
        assert hasattr(tu.fs.unit.overflow, "pressure")
        assert hasattr(tu.fs.unit.overflow, "alkalinity")

        assert number_variables(tu) == 87
        assert number_total_constraints(tu) == 63
        assert number_unused_variables(tu) == 6

    @pytest.mark.unit
    def test_dof(self, tu):
        assert degrees_of_freedom(tu) == 0

    @pytest.mark.unit
    def test_units(self, tu):
        assert_units_consistent(tu)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, tu):

        iscale.calculate_scaling_factors(tu)
        initialization_tester(tu)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, tu):
        solver = get_solver()
        results = solver.solve(tu)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, tu):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            tu.fs.unit.overflow.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            tu.fs.unit.overflow.temperature[0]
        )
        assert pytest.approx(0.003115, rel=1e-3) == value(
            tu.fs.unit.overflow.flow_vol[0]
        )
        assert pytest.approx(0.02806, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.000673, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_S"]
        )
        assert pytest.approx(0.06768, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0.001409, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_S"]
        )
        assert pytest.approx(0.04286, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_P"]
        )
        assert pytest.approx(0.099046, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_BH"]
        )
        assert pytest.approx(0.007414, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_BA"]
        )
        assert pytest.approx(0.001374, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_O"]
        )
        assert pytest.approx(0.009194, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_NO"]
        )
        assert pytest.approx(0.0001584, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_NH"]
        )
        assert pytest.approx(0.0005594, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "S_ND"]
        )
        assert pytest.approx(0.0001056, rel=1e-3) == value(
            tu.fs.unit.overflow.conc_mass_comp[0, "X_ND"]
        )
        assert pytest.approx(0.004564, rel=1e-3) == value(
            tu.fs.unit.overflow.alkalinity[0]
        )
        assert pytest.approx(3.82, rel=1e-3) == value(tu.fs.unit.height)
        assert pytest.approx(300, rel=1e-3) == value(tu.fs.unit.volume[0])
        assert pytest.approx(78.54, rel=1e-3) == value(tu.fs.unit.surface_area)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, tu):
        assert (
            abs(
                value(
                    tu.fs.unit.inlet.flow_vol[0] * tu.fs.props.dens_mass
                    - tu.fs.unit.overflow.flow_vol[0] * tu.fs.props.dens_mass
                    - tu.fs.unit.underflow.flow_vol[0] * tu.fs.props.dens_mass
                )
            )
            <= 1e-6
        )
        for i in tu.fs.props.solute_set:
            assert (
                abs(
                    value(
                        tu.fs.unit.inlet.flow_vol[0]
                        * tu.fs.unit.inlet.conc_mass_comp[0, i]
                        - tu.fs.unit.overflow.flow_vol[0]
                        * tu.fs.unit.overflow.conc_mass_comp[0, i]
                        - tu.fs.unit.underflow.flow_vol[0]
                        * tu.fs.unit.underflow.conc_mass_comp[0, i]
                    )
                )
                <= 1e-6
            )

    @pytest.mark.unit
    def test_report(self, tu):
        tu.fs.unit.report()


class TestThickASM2d(object):
    @pytest.fixture(scope="class")
    def tu_asm2d(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()

        m.fs.unit = Thickener(
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

        m.fs.unit.hydraulic_retention_time.fix()
        m.fs.unit.diameter.fix()

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, tu_asm2d):

        assert hasattr(tu_asm2d.fs.unit, "inlet")
        assert len(tu_asm2d.fs.unit.inlet.vars) == 5
        assert hasattr(tu_asm2d.fs.unit.inlet, "flow_vol")
        assert hasattr(tu_asm2d.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(tu_asm2d.fs.unit.inlet, "temperature")
        assert hasattr(tu_asm2d.fs.unit.inlet, "pressure")
        assert hasattr(tu_asm2d.fs.unit.inlet, "alkalinity")

        assert hasattr(tu_asm2d.fs.unit, "underflow")
        assert len(tu_asm2d.fs.unit.underflow.vars) == 5
        assert hasattr(tu_asm2d.fs.unit.underflow, "flow_vol")
        assert hasattr(tu_asm2d.fs.unit.underflow, "conc_mass_comp")
        assert hasattr(tu_asm2d.fs.unit.underflow, "temperature")
        assert hasattr(tu_asm2d.fs.unit.underflow, "pressure")
        assert hasattr(tu_asm2d.fs.unit.underflow, "alkalinity")

        assert hasattr(tu_asm2d.fs.unit, "overflow")
        assert len(tu_asm2d.fs.unit.overflow.vars) == 5
        assert hasattr(tu_asm2d.fs.unit.overflow, "flow_vol")
        assert hasattr(tu_asm2d.fs.unit.overflow, "conc_mass_comp")
        assert hasattr(tu_asm2d.fs.unit.overflow, "temperature")
        assert hasattr(tu_asm2d.fs.unit.overflow, "pressure")
        assert hasattr(tu_asm2d.fs.unit.overflow, "alkalinity")

        assert number_variables(tu_asm2d) == 111
        assert number_total_constraints(tu_asm2d) == 87
        assert number_unused_variables(tu_asm2d) == 0

    @pytest.mark.unit
    def test_dof(self, tu_asm2d):
        assert degrees_of_freedom(tu_asm2d) == 0

    @pytest.mark.unit
    def test_units(self, tu_asm2d):
        assert_units_consistent(tu_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, tu_asm2d):

        iscale.calculate_scaling_factors(tu_asm2d)
        initialization_tester(tu_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, tu_asm2d):
        solver = get_solver()
        results = solver.solve(tu_asm2d)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, tu_asm2d):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.temperature[0]
        )
        assert pytest.approx(0.00330, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.flow_vol[0]
        )
        assert pytest.approx(3.36066764339686e-05, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_A"]
        )
        assert pytest.approx(0.0004748063808766291, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_F"]
        )
        assert pytest.approx(0.03, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.029060299999999997, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_N2"]
        )
        assert pytest.approx(0.00802090132578769, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_NH4"]
        )
        tu_asm2d.fs.unit.overflow.conc_mass_comp.pprint()
        assert pytest.approx(0.006639502251179597, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_NO3"]
        )
        assert pytest.approx(0.007970701359416384, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_O2"]
        )
        assert pytest.approx(0.007895301409926407, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_PO4"]
        )
        assert pytest.approx(0.0024900686245326892, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_AUT"]
        )
        assert pytest.approx(0.03903683632802943, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_H"]
        )
        assert pytest.approx(0.03567622602578816, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0, abs=1e-6) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_MeOH"]
        )
        assert pytest.approx(0, abs=1e-6) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_MeP"]
        )
        assert pytest.approx(0.004513405145476244, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PAO"]
        )
        assert pytest.approx(5.761182100796471e-05, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PHA"]
        )
        assert pytest.approx(0.001336607105435848, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PP"]
        )
        assert pytest.approx(0.0014368738054293484, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_S"]
        )
        assert pytest.approx(0.0741692793990388, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_TSS"]
        )
        assert pytest.approx(0.004666303573014916, rel=1e-3) == value(
            tu_asm2d.fs.unit.overflow.alkalinity[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, tu_asm2d):
        assert (
            abs(
                value(
                    tu_asm2d.fs.unit.inlet.flow_vol[0] * tu_asm2d.fs.props.dens_mass
                    - tu_asm2d.fs.unit.overflow.flow_vol[0]
                    * tu_asm2d.fs.props.dens_mass
                    - tu_asm2d.fs.unit.underflow.flow_vol[0]
                    * tu_asm2d.fs.props.dens_mass
                )
            )
            <= 1e-6
        )
        for i in tu_asm2d.fs.props.solute_set:
            assert (
                abs(
                    value(
                        tu_asm2d.fs.unit.inlet.flow_vol[0]
                        * tu_asm2d.fs.unit.inlet.conc_mass_comp[0, i]
                        - tu_asm2d.fs.unit.overflow.flow_vol[0]
                        * tu_asm2d.fs.unit.overflow.conc_mass_comp[0, i]
                        - tu_asm2d.fs.unit.underflow.flow_vol[0]
                        * tu_asm2d.fs.unit.underflow.conc_mass_comp[0, i]
                    )
                )
                <= 1e-6
            )

    @pytest.mark.unit
    def test_report(self, tu_asm2d):
        tu_asm2d.fs.unit.report()


class TestThickModifiedASM2d(object):
    @pytest.fixture(scope="class")
    def tu_mod_asm2d(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()

        m.fs.unit = Thickener(
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

        m.fs.unit.hydraulic_retention_time.fix()
        m.fs.unit.diameter.fix()

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, tu_mod_asm2d):

        assert hasattr(tu_mod_asm2d.fs.unit, "inlet")
        assert len(tu_mod_asm2d.fs.unit.inlet.vars) == 4
        assert hasattr(tu_mod_asm2d.fs.unit.inlet, "flow_vol")
        assert hasattr(tu_mod_asm2d.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(tu_mod_asm2d.fs.unit.inlet, "temperature")
        assert hasattr(tu_mod_asm2d.fs.unit.inlet, "pressure")

        assert hasattr(tu_mod_asm2d.fs.unit, "underflow")
        assert len(tu_mod_asm2d.fs.unit.underflow.vars) == 4
        assert hasattr(tu_mod_asm2d.fs.unit.underflow, "flow_vol")
        assert hasattr(tu_mod_asm2d.fs.unit.underflow, "conc_mass_comp")
        assert hasattr(tu_mod_asm2d.fs.unit.underflow, "temperature")
        assert hasattr(tu_mod_asm2d.fs.unit.underflow, "pressure")

        assert hasattr(tu_mod_asm2d.fs.unit, "overflow")
        assert len(tu_mod_asm2d.fs.unit.overflow.vars) == 4
        assert hasattr(tu_mod_asm2d.fs.unit.overflow, "flow_vol")
        assert hasattr(tu_mod_asm2d.fs.unit.overflow, "conc_mass_comp")
        assert hasattr(tu_mod_asm2d.fs.unit.overflow, "temperature")
        assert hasattr(tu_mod_asm2d.fs.unit.overflow, "pressure")

        assert number_variables(tu_mod_asm2d) == 115
        assert number_total_constraints(tu_mod_asm2d) == 86
        assert number_unused_variables(tu_mod_asm2d) == 0

    @pytest.mark.unit
    def test_dof(self, tu_mod_asm2d):
        assert degrees_of_freedom(tu_mod_asm2d) == 0

    @pytest.mark.unit
    def test_units(self, tu_mod_asm2d):
        assert_units_consistent(tu_mod_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, tu_mod_asm2d):

        iscale.calculate_scaling_factors(tu_mod_asm2d)
        initialization_tester(tu_mod_asm2d)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, tu_mod_asm2d):
        solver = get_solver()
        results = solver.solve(tu_mod_asm2d)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, tu_mod_asm2d):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.temperature[0]
        )
        assert pytest.approx(0.0033139, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.flow_vol[0]
        )
        assert pytest.approx(3.36066764339686e-05, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_A"]
        )
        assert pytest.approx(0.0004748063808766291, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_F"]
        )
        assert pytest.approx(0.03, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.029060299999999997, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_N2"]
        )
        assert pytest.approx(0.00802090132578769, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_NH4"]
        )
        tu_mod_asm2d.fs.unit.overflow.conc_mass_comp.pprint()
        assert pytest.approx(0.006639502251179597, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_NO3"]
        )
        assert pytest.approx(0.007970701359416384, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_O2"]
        )
        assert pytest.approx(0.007895301409926407, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_PO4"]
        )
        assert pytest.approx(0.006999, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_K"]
        )
        assert pytest.approx(0.00599999, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_Mg"]
        )
        assert pytest.approx(0.0099999999, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "S_IC"]
        )
        assert pytest.approx(0.0024802, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_AUT"]
        )
        assert pytest.approx(0.0388828, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_H"]
        )
        assert pytest.approx(0.0355354585, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0.00449559, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PAO"]
        )
        assert pytest.approx(5.737786825381196e-05, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PHA"]
        )
        assert pytest.approx(0.001331327, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_PP"]
        )
        assert pytest.approx(0.0014311986, rel=1e-3) == value(
            tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, "X_S"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, tu_mod_asm2d):
        assert (
            abs(
                value(
                    tu_mod_asm2d.fs.unit.inlet.flow_vol[0]
                    * tu_mod_asm2d.fs.props.dens_mass
                    - tu_mod_asm2d.fs.unit.overflow.flow_vol[0]
                    * tu_mod_asm2d.fs.props.dens_mass
                    - tu_mod_asm2d.fs.unit.underflow.flow_vol[0]
                    * tu_mod_asm2d.fs.props.dens_mass
                )
            )
            <= 1e-6
        )
        for i in tu_mod_asm2d.fs.props.solute_set:
            assert (
                abs(
                    value(
                        tu_mod_asm2d.fs.unit.inlet.flow_vol[0]
                        * tu_mod_asm2d.fs.unit.inlet.conc_mass_comp[0, i]
                        - tu_mod_asm2d.fs.unit.overflow.flow_vol[0]
                        * tu_mod_asm2d.fs.unit.overflow.conc_mass_comp[0, i]
                        - tu_mod_asm2d.fs.unit.underflow.flow_vol[0]
                        * tu_mod_asm2d.fs.unit.underflow.conc_mass_comp[0, i]
                    )
                )
                <= 1e-6
            )

    @pytest.mark.unit
    def test_report(self, tu_mod_asm2d):
        tu_mod_asm2d.fs.unit.report()


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_costing():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ModifiedASM2dParameterBlock()

    m.fs.unit = Thickener(
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
    m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(214.5319 * units.mg / units.liter)
    m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(63.5316 * units.mg / units.liter)
    m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(2.7381 * units.mg / units.liter)
    m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(118.3582 * units.mg / units.liter)

    m.fs.unit.hydraulic_retention_time.fix()
    m.fs.unit.diameter.fix()

    m.fs.costing = WaterTAPCosting()

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)

    assert_optimal_termination(results)
    assert_units_consistent(m)
    assert hasattr(m.fs.costing, "thickener")
    assert value(m.fs.costing.thickener.capital_a_parameter) == 4729.8
    assert value(m.fs.costing.thickener.capital_b_parameter) == 37068

    # Check solutions
    assert pytest.approx(220675.79 * 2, rel=1e-5) == value(
        m.fs.unit.costing.capital_cost
    )
    assert pytest.approx(220675.79, rel=1e-5) == value(
        units.convert(
            (4729.8 * value(units.convert(10 * units.m, to_units=units.feet)) + 37068)
            * units.USD_2007,
            to_units=m.fs.costing.base_currency,
        )
    )
    assert pytest.approx(12.5 * 0.01255, rel=1e-5) == value(
        m.fs.unit.electricity_consumption[0]
    )
