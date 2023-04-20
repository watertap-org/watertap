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
Tests for dewatering unnit example.
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
    EnergyBalanceType,
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
    large_residuals_set
)

from idaes.core.util.scaling import (
    unscaled_variables_generator,
)

from idaes.core.util.testing import initialization_tester

from watertap.unit_models.dewatering import Dewatering_Unit
from watertap.property_models.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)

from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM1ParameterBlock()

    m.fs.unit = Dewatering_Unit(
        property_package=m.fs.props,
        outlet_list=["underflow", "overflow"],
        split_basis=SplittingType.componentFlow,
    )

    m.display()
    assert len(m.fs.unit.config) == 15

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert "underflow" in m.fs.unit.config.outlet_list
    assert "overflow" in m.fs.unit.config.outlet_list
    assert SplittingType.componentFlow is m.fs.unit.config.split_basis


# -----------------------------------------------------------------------------
class TestAdm(object):
    @pytest.fixture(scope="class")
    def adm(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM1ParameterBlock()

        m.fs.unit = Dewatering_Unit(
        property_package=m.fs.props,
        outlet_list=["underflow", "overflow"],
        split_basis=SplittingType.componentFlow,
    )

        m.fs.unit.inlet.flow_vol.fix(178.4674 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)

        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(130.867 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(258.5789 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(
            17216.2434 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(
            2611.4843 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(
            1e-6 * units.mg / units.liter
        )
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(626.0652 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(1e-6 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(1442.7882 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(0.54323 * units.mg / units.liter)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(100.8668 * units.mg / units.liter)
        m.fs.unit.inlet.alkalinity.fix(97.8459 * units.mol / units.m**3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, adm):

        assert hasattr(adm.fs.unit, "inlet")
        assert len(adm.fs.unit.inlet.vars) == 5
        assert hasattr(adm.fs.unit.inlet, "flow_vol")
        assert hasattr(adm.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(adm.fs.unit.inlet, "temperature")
        assert hasattr(adm.fs.unit.inlet, "pressure")
        assert hasattr(adm.fs.unit.inlet, "alkalinity")

        assert hasattr(adm.fs.unit, "underflow")
        assert len(adm.fs.unit.underflow.vars) == 5
        assert hasattr(adm.fs.unit.underflow, "flow_vol")
        assert hasattr(adm.fs.unit.underflow, "conc_mass_comp")
        assert hasattr(adm.fs.unit.underflow, "temperature")
        assert hasattr(adm.fs.unit.underflow, "pressure")
        assert hasattr(adm.fs.unit.underflow, "alkalinity")

        assert hasattr(adm.fs.unit, "overflow")
        assert len(adm.fs.unit.overflow.vars) == 5
        assert hasattr(adm.fs.unit.overflow, "flow_vol")
        assert hasattr(adm.fs.unit.overflow, "conc_mass_comp")
        assert hasattr(adm.fs.unit.overflow, "temperature")
        assert hasattr(adm.fs.unit.overflow, "pressure")
        assert hasattr(adm.fs.unit.overflow, "alkalinity")

        assert number_variables(adm) == 76
        assert number_total_constraints(adm) == 60
        assert number_unused_variables(adm) == 0

    @pytest.mark.unit
    def test_dof(self, adm):
        assert degrees_of_freedom(adm) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, adm):
        initialization_tester(adm)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, adm):
        solver = get_solver()
        results = solver.solve(adm)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, adm):
        assert pytest.approx(101325.0, rel=1e-3) == value(
            adm.fs.unit.overflow.pressure[0]
        )
        assert pytest.approx(308.15, rel=1e-3) == value(
            adm.fs.unit.overflow.temperature[0]
        )
        assert pytest.approx(0.001954, rel=1e-3) == value(
            adm.fs.unit.overflow.flow_vol[0]
        )
        assert pytest.approx(0.1308, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.2585, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_S"]
        )
        assert pytest.approx(0.3638, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(0.0552, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_S"]
        )
        assert value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_BH"]) <= 1e-6
        assert value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_BA"]) <= 1e-6
        
        assert pytest.approx(0.01323, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_P"]
        )
        assert value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_O"]) <= 1e-6
        assert value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_NO"]) <= 1e-6
        
        assert pytest.approx(1.44278, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_NH"]
        )
        assert pytest.approx(0.0005439, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "S_ND"]
        )
        assert pytest.approx(0.002132, rel=1e-3) == value(
            adm.fs.unit.overflow.conc_mass_comp[0, "X_ND"]
        )
        assert pytest.approx(0.097845, rel=1e-3) == value(
            adm.fs.unit.overflow.alkalinity[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, adm):
        assert (
            abs(
                value(
                    adm.fs.unit.inlet.flow_vol[0] * adm.fs.props.dens_mass
                    - adm.fs.unit.overflow.flow_vol[0] * adm.fs.props.dens_mass
                    - adm.fs.unit.underflow.flow_vol[0] * adm.fs.props.dens_mass
                )
            )
            <= 1e-6
        )

    @pytest.mark.unit
    def test_report(self, adm):
        adm.fs.unit.report()
