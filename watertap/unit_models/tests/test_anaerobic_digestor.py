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
Tests for anaerobic digestor example.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

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

from idaes.core.util.scaling import (
    unscaled_variables_generator,
)

from idaes.core.util.testing import initialization_tester

from watertap.unit_models.anaerobic_digestor import AD
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

    m.fs.unit = AD(
        liquid_property_package=m.fs.props,
        vapor_property_package=m.fs.props_vap,
        reaction_package=m.fs.rxn_props,
    )

    assert len(m.fs.unit.config) == 16

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.liquid_property_package is m.fs.props
    assert m.fs.unit.config.vapor_property_package is m.fs.props_vap
    assert m.fs.unit.config.reaction_package is m.fs.rxn_props


# -----------------------------------------------------------------------------
class TestAdm(object):
    @pytest.fixture(scope="class")
    def adm(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.props_vap = ADM1_vaporParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.unit = AD(
            liquid_property_package=m.fs.props,
            vapor_property_package=m.fs.props_vap,
            reaction_package=m.fs.rxn_props,
            has_heat_transfer=True,
            has_pressure_change=False,
        )

        m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature.fix(308.15)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

        m.fs.unit.inlet.cations[0].fix(0.04)
        m.fs.unit.inlet.anions[0].fix(0.02)

        m.fs.unit.volume_liquid.fix(3400)
        m.fs.unit.volume_vapor.fix(300)
        m.fs.unit.liquid_outlet.temperature.fix(308.15)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, adm):

        assert hasattr(adm.fs.unit, "inlet")
        assert len(adm.fs.unit.inlet.vars) == 6
        assert hasattr(adm.fs.unit.inlet, "flow_vol")
        assert hasattr(adm.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(adm.fs.unit.inlet, "temperature")
        assert hasattr(adm.fs.unit.inlet, "pressure")
        assert hasattr(adm.fs.unit.inlet, "anions")
        assert hasattr(adm.fs.unit.inlet, "cations")

        assert hasattr(adm.fs.unit, "liquid_outlet")
        assert len(adm.fs.unit.liquid_outlet.vars) == 6
        assert hasattr(adm.fs.unit.liquid_outlet, "flow_vol")
        assert hasattr(adm.fs.unit.liquid_outlet, "conc_mass_comp")
        assert hasattr(adm.fs.unit.liquid_outlet, "temperature")
        assert hasattr(adm.fs.unit.liquid_outlet, "pressure")
        assert hasattr(adm.fs.unit.liquid_outlet, "anions")
        assert hasattr(adm.fs.unit.liquid_outlet, "cations")

        assert hasattr(adm.fs.unit, "vapor_outlet")
        assert len(adm.fs.unit.vapor_outlet.vars) == 4
        assert hasattr(adm.fs.unit.vapor_outlet, "flow_vol")
        assert hasattr(adm.fs.unit.vapor_outlet, "conc_mass_comp")
        assert hasattr(adm.fs.unit.vapor_outlet, "temperature")
        assert hasattr(adm.fs.unit.vapor_outlet, "pressure")

        assert hasattr(adm.fs.unit, "ad_performance_eqn")
        assert hasattr(adm.fs.unit, "volume_AD")
        assert hasattr(adm.fs.unit, "volume_liquid")
        assert hasattr(adm.fs.unit, "volume_vapor")
        assert hasattr(adm.fs.unit, "heat_duty")

        assert number_variables(adm.fs.unit) == 201
        assert number_total_constraints(adm.fs.unit) == 169
        assert number_unused_variables(adm.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, adm):
        assert_units_consistent(adm.fs.unit)
        assert_units_equivalent(adm.fs.unit.volume_AD[0], units.m**3)
        assert_units_equivalent(adm.fs.unit.heat_duty[0], units.W)

    @pytest.mark.unit
    def test_dof(self, adm):
        assert degrees_of_freedom(adm) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, adm):
        initialization_tester(adm)

    @pytest.mark.component
    def test_var_scaling(self, adm):

        unscaled_var_list = list(
            unscaled_variables_generator(adm.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

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
        assert pytest.approx(101325.0, abs=1e-0) == value(
            adm.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(308.15, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.temperature[0]
        )
        assert pytest.approx(0.328772, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_I"]
        )
        assert pytest.approx(0.00531408, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_aa"]
        )
        assert pytest.approx(0.197783, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ac"]
        )
        assert pytest.approx(0.0132484, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_bu"]
        )
        assert pytest.approx(0.0549707, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ch4"]
        )
        assert pytest.approx(0.0986058, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_fa"]
        )
        assert pytest.approx(2.35916e-07, rel=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_h2"]
        )
        assert pytest.approx(0.01578123, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_pro"]
        )
        assert pytest.approx(0.0119533, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_su"]
        )
        assert pytest.approx(0.0116230, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_va"]
        )
        assert pytest.approx(25.6217, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_I"]
        )
        assert pytest.approx(1.1793, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_aa"]
        )
        assert pytest.approx(0.760653, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ac"]
        )
        assert pytest.approx(0.308718, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c"]
        )
        assert pytest.approx(0.431974, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c4"]
        )
        assert pytest.approx(0.0279475, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ch"]
        )
        assert pytest.approx(0.243068, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_fa"]
        )
        assert pytest.approx(0.3170629, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_h2"]
        )
        assert pytest.approx(0.0294834, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_li"]
        )
        assert pytest.approx(0.102574, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pr"]
        )
        assert pytest.approx(0.137323, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pro"]
        )
        assert pytest.approx(0.420219, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "X_su"]
        )
        assert pytest.approx(1.8321, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IC"]
        )
        assert pytest.approx(1.8232, abs=1e-3) == value(
            adm.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IN"]
        )
        assert pytest.approx(0.02, abs=1e-2) == value(
            adm.fs.unit.liquid_outlet.anions[0]
        )
        assert pytest.approx(0.04, abs=1e-2) == value(
            adm.fs.unit.liquid_outlet.cations[0]
        )
        assert pytest.approx(106659, abs=1e-0) == value(
            adm.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(308.15, abs=1e-3) == value(
            adm.fs.unit.vapor_outlet.temperature[0]
        )
        assert pytest.approx(0.034, abs=1e-2) == value(
            adm.fs.unit.vapor_outlet.flow_vol[0]
        )
        assert pytest.approx(1.6216, abs=1e-3) == value(
            adm.fs.unit.vapor_outlet.conc_mass_comp[0, "S_ch4"]
        )
        assert pytest.approx(0.169417, abs=1e-3) == value(
            adm.fs.unit.vapor_outlet.conc_mass_comp[0, "S_co2"]
        )
        assert pytest.approx(0.0271, abs=1e-3) == value(adm.fs.unit.KH_co2[0])
        assert pytest.approx(0.00116, abs=1e-3) == value(adm.fs.unit.KH_ch4[0])
        assert pytest.approx(7.38e-4, rel=1e-2) == value(adm.fs.unit.KH_h2[0])
        assert pytest.approx(0.2054, abs=1e-3) == value(
            adm.fs.unit.electricity_consumption[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, adm):
        assert (
            abs(
                value(
                    adm.fs.unit.inlet.flow_vol[0] * adm.fs.props.dens_mass
                    - adm.fs.unit.liquid_outlet.flow_vol[0] * adm.fs.props.dens_mass
                    - adm.fs.unit.vapor_outlet.flow_vol[0] * adm.fs.props_vap.dens_mass
                )
            )
            <= 1e-6
        )

        assert pytest.approx(-13.5835, abs=1e-3) == value(
            adm.fs.unit.liquid_phase.enthalpy_transfer[0]
        )
        assert (
            abs(
                value(
                    (
                        adm.fs.unit.inlet.flow_vol[0]
                        * adm.fs.props.dens_mass
                        * adm.fs.props.cp_mass
                        * (
                            adm.fs.unit.inlet.temperature[0]
                            - adm.fs.props.temperature_ref
                        )
                    )
                    - (
                        adm.fs.unit.liquid_outlet.flow_vol[0]
                        * adm.fs.props.dens_mass
                        * adm.fs.props.cp_mass
                        * (
                            adm.fs.unit.liquid_outlet.temperature[0]
                            - adm.fs.props.temperature_ref
                        )
                    )
                    + adm.fs.unit.liquid_phase.enthalpy_transfer[0]
                )
            )
            <= 1e-2
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, adm):
        m = adm

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.liquid_phase.properties_out[0].flow_vol)
        solver = get_solver()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1083290.8, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(5.04295, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.unit
    def test_get_performance_contents(self, adm):
        perf_dict = adm.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Volume": adm.fs.unit.volume_AD[0],
                "Heat Duty": adm.fs.unit.heat_duty[0],
            }
        }

    @pytest.mark.unit
    def test_report(self, adm):
        adm.fs.unit.report()
