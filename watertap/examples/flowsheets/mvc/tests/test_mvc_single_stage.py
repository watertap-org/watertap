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
    Expression
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Feed, Product, Mixer, Separator
from idaes.models.unit_models.heat_exchanger import HeatExchanger
from idaes.models.unit_models.translator import Translator
from pyomo.util.check_units import assert_units_consistent

from watertap.examples.flowsheets.mvc.mvc_single_stage import (
    build,
    set_operating_conditions,
    add_Q_ext,
    solve,
    set_up_optimization,
    display_metrics,
    display_design
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
        assert isinstance(m.fs.Q_ext,Var)
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
            "maintenance_labor_chemical_operating_cost",
            "total_operating_cost",
        ]
        for var_str in var_str_list:
            var = getattr(m.fs.costing, var_str)
            assert isinstance(var, Var)

        # arcs
        arc_dict = {
            m.fs.s01: (m.fs.feed.outlet, m.fs.pump_feed.inlet),
            m.fs.s02: (m.fs.pump_feed.outlet, m.fs.separator_feed.inlet),
            m.fs.s03: (m.fs.separator_feed.hx_distillate_cold, m.fs.hx_distillate.cold_inlet),
            m.fs.s04: (m.fs.separator_feed.hx_brine_cold, m.fs.hx_brine.cold_inlet),
            m.fs.s05: (m.fs.hx_distillate.cold_outlet, m.fs.mixer_feed.hx_distillate_cold),
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
            m.fs.s16: (m.fs.hx_distillate.hot_outlet, m.fs.distillate.inlet)
        }
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

        # units
        assert_units_consistent(m.fs)