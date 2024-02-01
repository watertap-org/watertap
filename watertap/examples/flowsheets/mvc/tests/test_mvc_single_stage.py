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
    ConcreteModel,
    Block,
    Var,
    Constraint,
    Expression,
    value,
    Objective,
    assert_optimal_termination,
)
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Feed, Product, Mixer, Separator
from idaes.models.unit_models.heat_exchanger import HeatExchanger
from idaes.models.unit_models.translator import Translator
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.model_statistics import degrees_of_freedom, number_total_objectives

from watertap.examples.flowsheets.mvc.mvc_single_stage import (
    build,
    set_operating_conditions,
    add_Q_ext,
    initialize_system,
    scale_costs,
    fix_outlet_pressures,
    solve,
    set_up_optimization,
    main,
)

import watertap.property_models.seawater_prop_pack as props_seawater
import watertap.property_models.water_prop_pack as props_water
from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser


solver = get_solver()


class TestMVC:
    @pytest.fixture(scope="class")
    def mvc_single_stage(self):
        m = build()
        add_Q_ext(m, time_point=m.fs.config.time)
        return m

    @pytest.mark.unit
    def test_build_model(self, mvc_single_stage):
        m = mvc_single_stage

        # test model set up
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties_feed, props_seawater.SeawaterParameterBlock)
        assert isinstance(m.fs.properties_vapor, props_water.WaterParameterBlock)
        assert isinstance(m.fs.costing, Block)

        # test unit models
        assert isinstance(m.fs.feed, Feed)
        assert isinstance(m.fs.distillate, Product)
        assert isinstance(m.fs.brine, Product)
        assert isinstance(m.fs.pump_feed, Pump)
        assert isinstance(m.fs.pump_distillate, Pump)
        assert isinstance(m.fs.pump_brine, Pump)
        assert isinstance(m.fs.separator_feed, Separator)
        assert isinstance(m.fs.hx_distillate, HeatExchanger)
        assert isinstance(m.fs.hx_brine, HeatExchanger)
        assert isinstance(m.fs.mixer_feed, Mixer)
        assert isinstance(m.fs.evaporator, Evaporator)
        assert isinstance(m.fs.compressor, Compressor)
        assert isinstance(m.fs.condenser, Condenser)
        assert isinstance(m.fs.tb_distillate, Translator)

        # unit model options
        assert isinstance(m.fs.hx_distillate.cold.deltaP, Var)
        assert isinstance(m.fs.hx_distillate.hot.deltaP, Var)
        assert isinstance(m.fs.hx_brine.cold.deltaP, Var)
        assert isinstance(m.fs.hx_brine.hot.deltaP, Var)

        # additional constraints, variables, and expressions
        assert isinstance(m.fs.recovery, Var)
        assert isinstance(m.fs.recovery_equation, Constraint)
        assert isinstance(m.fs.split_ratio_recovery_equality, Constraint)
        assert isinstance(m.fs.Q_ext, Var)
        assert isinstance(m.fs.costing.annual_water_production, Expression)
        assert isinstance(m.fs.costing.specific_energy_consumption, Expression)

        # costing blocks
        for blk_str in ("evaporator", "compressor", "hx_distillate", "hx_brine"):
            blk = getattr(m.fs, blk_str)
            c_blk = getattr(blk, "costing")
            assert isinstance(c_blk, Block)
            assert isinstance(getattr(c_blk, "capital_cost"), Var)

        var_str_list = [
            "total_capital_cost",
            "total_operating_cost",
        ]
        for var_str in var_str_list:
            var = getattr(m.fs.costing, var_str)
            assert isinstance(var, Var)

        # arcs
        arc_dict = {
            m.fs.s01: (m.fs.feed.outlet, m.fs.pump_feed.inlet),
            m.fs.s02: (m.fs.pump_feed.outlet, m.fs.separator_feed.inlet),
            m.fs.s03: (
                m.fs.separator_feed.hx_distillate_cold,
                m.fs.hx_distillate.cold_inlet,
            ),
            m.fs.s04: (m.fs.separator_feed.hx_brine_cold, m.fs.hx_brine.cold_inlet),
            m.fs.s05: (
                m.fs.hx_distillate.cold_outlet,
                m.fs.mixer_feed.hx_distillate_cold,
            ),
            m.fs.s06: (m.fs.hx_brine.cold_outlet, m.fs.mixer_feed.hx_brine_cold),
            m.fs.s07: (m.fs.mixer_feed.outlet, m.fs.evaporator.inlet_feed),
            m.fs.s08: (m.fs.evaporator.outlet_vapor, m.fs.compressor.inlet),
            m.fs.s09: (m.fs.compressor.outlet, m.fs.condenser.inlet),
            m.fs.s10: (m.fs.evaporator.outlet_brine, m.fs.pump_brine.inlet),
            m.fs.s11: (m.fs.pump_brine.outlet, m.fs.hx_brine.hot_inlet),
            m.fs.s12: (m.fs.hx_brine.hot_outlet, m.fs.brine.inlet),
            m.fs.s13: (m.fs.condenser.outlet, m.fs.tb_distillate.inlet),
            m.fs.s14: (m.fs.tb_distillate.outlet, m.fs.pump_distillate.inlet),
            m.fs.s15: (m.fs.pump_distillate.outlet, m.fs.hx_distillate.hot_inlet),
            m.fs.s16: (m.fs.hx_distillate.hot_outlet, m.fs.distillate.inlet),
        }
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

        # units
        assert_units_consistent(m.fs)

    @pytest.mark.component
    def test_set_operating_conditions(self, mvc_single_stage):
        m = mvc_single_stage
        set_operating_conditions(m)

        # check fixed variables
        # feed
        assert m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].is_fixed()
        assert value(m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]) == 0.1
        assert m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].is_fixed()
        assert value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]) == 40
        assert m.fs.feed.temperature[0].is_fixed()
        assert value(m.fs.feed.temperature[0]) == 298.15
        assert m.fs.feed.pressure[0].is_fixed()
        assert value(m.fs.feed.pressure[0]) == 101325

        # recovery
        assert m.fs.recovery[0].is_fixed()
        assert value(m.fs.recovery[0]) == 0.5

        # pumps
        assert m.fs.pump_feed.efficiency_pump[0].is_fixed()
        assert value(m.fs.pump_feed.efficiency_pump[0]) == 0.8
        assert m.fs.pump_feed.control_volume.deltaP[0].is_fixed()
        assert value(m.fs.pump_feed.control_volume.deltaP[0]) == 7e3

        assert m.fs.pump_distillate.efficiency_pump[0].is_fixed()
        assert value(m.fs.pump_distillate.efficiency_pump[0]) == 0.8
        assert m.fs.pump_distillate.control_volume.deltaP[0].is_fixed()
        assert value(m.fs.pump_distillate.control_volume.deltaP[0]) == 4e4

        assert m.fs.pump_brine.efficiency_pump[0].is_fixed()
        assert value(m.fs.pump_brine.efficiency_pump[0]) == 0.8
        assert m.fs.pump_brine.control_volume.deltaP[0].is_fixed()
        assert value(m.fs.pump_brine.control_volume.deltaP[0]) == 4e4

        # heat exchangers
        assert m.fs.hx_distillate.overall_heat_transfer_coefficient[0].is_fixed()
        assert value(m.fs.hx_distillate.overall_heat_transfer_coefficient[0]) == 2e3
        assert m.fs.hx_distillate.area.is_fixed()
        assert value(m.fs.hx_distillate.area) == 125
        assert m.fs.hx_distillate.cold.deltaP[0].is_fixed()
        assert value(m.fs.hx_distillate.cold.deltaP[0]) == 7e3
        assert m.fs.hx_distillate.hot.deltaP[0].is_fixed()
        assert value(m.fs.hx_distillate.hot.deltaP[0]) == 7e3

        assert m.fs.hx_brine.overall_heat_transfer_coefficient[0].is_fixed()
        assert value(m.fs.hx_brine.overall_heat_transfer_coefficient[0]) == 2e3
        assert m.fs.hx_brine.area.is_fixed()
        assert value(m.fs.hx_brine.area) == 115
        assert m.fs.hx_brine.cold.deltaP[0].is_fixed()
        assert value(m.fs.hx_brine.cold.deltaP[0]) == 7e3
        assert m.fs.hx_brine.hot.deltaP[0].is_fixed()
        assert value(m.fs.hx_brine.hot.deltaP[0]) == 7e3

        # evaporator
        assert m.fs.evaporator.outlet_brine.temperature[0].is_fixed()
        assert value(m.fs.evaporator.outlet_brine.temperature[0]) == 343.15
        assert m.fs.evaporator.U.is_fixed()
        assert value(m.fs.evaporator.U) == 3e3

        # compressor
        assert m.fs.compressor.pressure_ratio.is_fixed()
        assert value(m.fs.compressor.pressure_ratio) == 1.6
        assert m.fs.compressor.efficiency.is_fixed()
        assert value(m.fs.compressor.efficiency) == 0.8

        # 0 TDS in distillate
        assert (
            m.fs.tb_distillate.properties_out[0]
            .flow_mass_phase_comp["Liq", "TDS"]
            .is_fixed()
        )
        assert (
            value(
                m.fs.tb_distillate.properties_out[0].flow_mass_phase_comp["Liq", "TDS"]
            )
            == 1e-5
        )

        # Costing
        assert m.fs.costing.TIC.is_fixed()
        assert value(m.fs.costing.TIC) == 2
        assert m.fs.costing.heat_exchanger.material_factor_cost.is_fixed()
        assert value(m.fs.costing.heat_exchanger.material_factor_cost) == 5
        assert m.fs.costing.evaporator.material_factor_cost.is_fixed()
        assert value(m.fs.costing.evaporator.material_factor_cost) == 5
        assert m.fs.costing.compressor.unit_cost.is_fixed()
        assert value(m.fs.costing.compressor.unit_cost) == 7364

        # Temperature upper bounds
        assert value(m.fs.evaporator.properties_vapor[0].temperature.ub) == 348.15
        assert (
            value(m.fs.compressor.control_volume.properties_out[0].temperature.ub)
            == 450
        )

        # check degrees of freedom
        assert degrees_of_freedom(m) == 1

    @pytest.mark.component
    @pytest.mark.requires_idaes_solver
    def test_initialize_system(self, mvc_single_stage):
        m = mvc_single_stage

        initialize_system(m)

        # mass flows in evaporator
        assert value(
            m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
        ) == pytest.approx(22.222, rel=1e-3)
        assert value(
            m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"]
        ) == pytest.approx(1e-8, abs=1e-3)
        assert value(
            m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
        ) == pytest.approx(17.777, rel=1e-3)
        assert value(
            m.fs.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(4.444, rel=1e-3)

        # evaporator pressure
        assert value(m.fs.evaporator.properties_vapor[0].pressure) == pytest.approx(
            26.231e3, rel=1e-3
        )

        # compressor pressure outlet
        assert value(m.fs.compressor.control_volume.properties_out[0].pressure) / value(
            m.fs.compressor.control_volume.properties_in[0].pressure
        ) == pytest.approx(1.6, rel=1e-3)

    @pytest.mark.component
    @pytest.mark.requires_idaes_solver
    def test_simulation_Q_ext(self, mvc_single_stage):
        m = mvc_single_stage

        scale_costs(m)
        fix_outlet_pressures(m)

        assert degrees_of_freedom(m) == 1

        m.fs.objective = Objective(expr=m.fs.Q_ext[0])
        results = solve(m, solver=solver, tee=False)
        assert_optimal_termination(results)

        # Check system metrics
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            22.02, rel=1e-2
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(23.47, rel=1e-2)

        # Check mass balance
        assert pytest.approx(
            value(m.fs.feed.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]), rel=1e-3
        ) == value(m.fs.distillate.inlet.flow_mass_phase_comp[0, "Liq", "H2O"]) + value(
            m.fs.brine.inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(
            value(m.fs.feed.outlet.flow_mass_phase_comp[0, "Liq", "TDS"]), rel=1e-3
        ) == value(m.fs.distillate.inlet.flow_mass_phase_comp[0, "Liq", "TDS"]) + value(
            m.fs.brine.inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_optimization(self, mvc_single_stage):
        m = mvc_single_stage
        m.fs.Q_ext[0].fix(0)  # no longer want external heating in evaporator
        del m.fs.objective
        set_up_optimization(m)
        assert number_total_objectives(m) == 1
        assert degrees_of_freedom(m) == 4
        results = solve(m, solver=solver, tee=False)
        assert_optimal_termination(results)
        # Check decision variables
        assert value(m.fs.evaporator.properties_brine[0].temperature) == pytest.approx(
            348.15, rel=1e-2
        )
        assert value(m.fs.evaporator.properties_brine[0].pressure) == pytest.approx(
            32448.24, rel=1e-2
        )
        assert value(m.fs.hx_brine.area) == pytest.approx(173.99, rel=1e-2)
        assert value(m.fs.hx_distillate.area) == pytest.approx(206.31, rel=1e-2)
        assert value(m.fs.compressor.pressure_ratio) == pytest.approx(1.61, rel=1e-2)
        assert value(m.fs.evaporator.area) == pytest.approx(777.37, rel=1e-2)
        assert value(m.fs.evaporator.lmtd) == pytest.approx(22.59, rel=1e-2)

        # Check system metrics
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            22.36, rel=1e-2
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(4.52, rel=1e-2)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.unit
    def test_main_fun(self):
        main()
