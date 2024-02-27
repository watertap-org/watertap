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

import pytest
from pyomo.environ import (
    Block,
    ConcreteModel,
    assert_optimal_termination,
    Var,
    value,
    units,
)
from idaes.core import FlowsheetBlock
from watertap.unit_models.anaerobic_digester import AD
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)
from watertap.unit_models.tests.test_anaerobic_digester import build as AD_frame
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import calculate_scaling_factors
from pyomo.util.check_units import assert_units_consistent
from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def test_costing():
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
    # m.dummy_unit.area = Var()
    # m.dummy_unit.area.fix(10)

    # Add unit model costing
    m.fs.costing = WaterTAPCosting()

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()

    iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-6)

    assert value(m.fs.unit.costing.capital_cost) == pytest.approx(100000, rel=1e-2)


# class TestElectroNP:
#     @pytest.fixture(scope="class")
#     def build(self):
#         m = AD_frame
#
#         return m
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_initialize(self, ElectroNP_frame):
#         initialization_tester(ElectroNP_frame)
#
#         # check that all variables have scaling factors
#         unscaled_var_list = list(iscale.unscaled_variables_generator(ElectroNP_frame))
#         assert len(unscaled_var_list) == 0
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_solve(self, ElectroNP_frame):
#         m = ElectroNP_frame
#         results = solver.solve(m)
#
#         # Check for optimal solution
#         assert_optimal_termination(results)
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_conservation(self, ElectroNP_frame):
#         m = ElectroNP_frame
#         assert (
#             abs(
#                 value(
#                     m.fs.unit.inlet.flow_vol[0] * m.fs.properties.dens_mass
#                     - m.fs.unit.treated.flow_vol[0] * m.fs.properties.dens_mass
#                     - m.fs.unit.byproduct.flow_vol[0] * m.fs.properties.dens_mass
#                 )
#             )
#             <= 1e-6
#         )
#         for j in m.fs.properties.solute_set:
#             assert 1e-6 >= abs(
#                 value(
#                     m.fs.unit.inlet.flow_vol[0] * m.fs.unit.inlet.conc_mass_comp[0, j]
#                     - m.fs.unit.treated.flow_vol[0]
#                     * m.fs.unit.treated.conc_mass_comp[0, j]
#                     - m.fs.unit.byproduct.flow_vol[0]
#                     * m.fs.unit.byproduct.conc_mass_comp[0, j]
#                 )
#             )
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_solution(self, ElectroNP_frame):
#         m = ElectroNP_frame
#         assert value(m.fs.unit.treated.flow_vol[0]) == pytest.approx(0.213495, rel=1e-3)
#
#         assert value(m.fs.unit.treated.temperature[0]) == pytest.approx(
#             298.15, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_A"]) == pytest.approx(
#             0.02, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_F"]) == pytest.approx(
#             0.03, rel=1e-2
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
#             0.03, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
#             0.0112, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_O2"]) == pytest.approx(
#             0.01, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
#             7.2e-5, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
#             0.03, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
#             0.025, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_PP"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
#             0.125, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.treated.conc_mass_comp[0, "S_K"]) == pytest.approx(
#             0, abs=1e-4
#         )
#         assert value(m.fs.unit.byproduct.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
#             480, rel=1e-4
#         )
#         assert value(m.fs.unit.byproduct.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
#             352.8, rel=1e-4
#         )
#         assert value(m.fs.unit.treated.alkalinity[0]) == pytest.approx(
#             0.005083, rel=1e-4
#         )
#         assert value(m.fs.unit.energy_electric_flow_mass) == pytest.approx(
#             0.044, rel=1e-4
#         )
#         assert value(m.fs.unit.electricity[0]) == pytest.approx(0.1193, rel=1e-4)
#         assert value(m.fs.unit.MgCl2_flowrate[0]) == pytest.approx(1.0521, rel=1e-4)
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_costing(self, ElectroNP_frame):
#         m = ElectroNP_frame
#
#         m.fs.costing = WaterTAPCosting()
#
#         m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
#         m.fs.costing.cost_process()
#         m.fs.costing.add_LCOW(m.fs.unit.properties_treated[0].flow_vol)
#         results = solver.solve(m)
#
#         assert_optimal_termination(results)
#
#         # Check solutions
#         assert pytest.approx(2.0 * 1036611.9, rel=1e-5) == value(
#             m.fs.unit.costing.capital_cost
#         )
#         assert pytest.approx(0.04431857, rel=1e-5) == value(m.fs.costing.LCOW)
#
#     @pytest.mark.unit
#     def test_report(self, ElectroNP_frame):
#         m = ElectroNP_frame
#         m.fs.unit.report()
