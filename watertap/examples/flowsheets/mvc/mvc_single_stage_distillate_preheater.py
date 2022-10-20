###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
from pyomo.environ import (
    ConcreteModel,
    value,
    Set,
    Constraint,
    Expression,
    Objective,
    Param,
    Var,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
    check_optimal_termination
)
from pyomo.network import Arc

from pyomo.util.check_units import assert_units_consistent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_diagnostics import DegeneracyHunter
from idaes.models.unit_models import Feed, Separator, Mixer, Product
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.generic_models.costing import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.costing import WaterTAPCosting, PumpType
import math
import pandas as pd

def main():
    # build
    m = build()
    add_Q_ext(m, time_point=m.fs.config.time)
    set_operating_conditions(m)
    initialize_system(m)
    scale_costs(m)
    fix_outlet_pressures(m)

    print('DOF after initialization: ', degrees_of_freedom(m))
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    #m.fs.objective = Objective(expr=1)
    solver = get_solver()
    results = solve(m,tee=False)
    print('First solve - simulation')
    print(results.solver.termination_condition)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)
    display_results(m)

    print('Second solve - optimize')
    # m.fs.recovery[0].fix(0.7)
    m.fs.Q_ext[0].fix(0)
    del m.fs.objective
    set_up_optimization(m)
    results = solve(m, tee=False)
    print(results.solver.termination_condition)
    display_results(m)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)

    return m

def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties_feed})
    m.fs.pump_feed = Pump(default={"property_package": m.fs.properties_feed})

    m.fs.hx_distillate = HeatExchanger(
        default={
            "hot_side_name": "hot",
            "cold_side_name": "cold",
            "hot": {"property_package": m.fs.properties_feed},
            "cold": {"property_package": m.fs.properties_feed},
            "delta_temperature_callback": delta_temperature_chen_callback,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
        }
    )
    # Set lower bound of approach temperatures
    m.fs.hx_distillate.delta_temperature_in.setlb(0)
    m.fs.hx_distillate.delta_temperature_out.setlb(0)
    m.fs.hx_distillate.area.setlb(10)
    add_pressure_drop_to_hx(m.fs.hx_distillate, m.fs.config.time)

    m.fs.evaporator = Evaporator(
        default={
            "property_package_feed": m.fs.properties_feed,
            "property_package_vapor": m.fs.properties_vapor,
        }
    )

    m.fs.compressor = Compressor(default={"property_package": m.fs.properties_vapor})

    m.fs.condenser = Condenser(default={"property_package": m.fs.properties_vapor})

    m.fs.tb_distillate = Translator(
        default={
            "inlet_property_package": m.fs.properties_vapor,
            "outlet_property_package": m.fs.properties_feed
        }
    )

    # Translator block to convert distillate exiting condenser from water to seawater prop pack
    @m.fs.tb_distillate.Constraint()
    def eq_flow_mass_comp(blk):
        return (
            blk.properties_in[0].flow_mass_phase_comp['Liq','H2O']
            == blk.properties_out[0].flow_mass_phase_comp["Liq", 'H2O']
        )

    @m.fs.tb_distillate.Constraint()
    def eq_temperature(blk):
        return (
                blk.properties_in[0].temperature
                == blk.properties_out[0].temperature
        )

    @m.fs.tb_distillate.Constraint()
    def eq_pressure(blk):
        return (
                blk.properties_in[0].pressure
                == blk.properties_out[0].pressure
        )

    m.fs.pump_brine = Pump(default={"property_package": m.fs.properties_feed})

    m.fs.pump_distillate = Pump(default={"property_package": m.fs.properties_feed})

    m.fs.distillate = Product(default={"property_package": m.fs.properties_feed})

    m.fs.brine = Product(default={"property_package": m.fs.properties_feed})

    # Connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_feed.inlet)
    # m.fs.s02 = Arc(source=m.fs.pump_feed.outlet, destination=m.fs.separator_feed.inlet)
    m.fs.s03 = Arc(
        source=m.fs.pump_feed.outlet,
        destination=m.fs.hx_distillate.cold_inlet,
    )
    # m.fs.s04 = Arc(
    #     source=m.fs.separator_feed.hx_brine_cold, destination=m.fs.hx_brine.cold_inlet
    # )
    # m.fs.s05 = Arc(
    #     source=m.fs.hx_distillate.cold_outlet,
    #     destination=m.fs.mixer_feed.hx_distillate_cold,
    # )
    # m.fs.s06 = Arc(
    #     source=m.fs.hx_brine.cold_outlet, destination=m.fs.mixer_feed.hx_brine_cold
    # )
    m.fs.s07 = Arc(
        source=m.fs.hx_distillate.cold_outlet, destination=m.fs.evaporator.inlet_feed
    )
    m.fs.s08 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s09 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    m.fs.s10 = Arc(
        source=m.fs.evaporator.outlet_brine, destination=m.fs.pump_brine.inlet
    )
    m.fs.s11 = Arc(source=m.fs.pump_brine.outlet, destination=m.fs.brine.inlet)
    # m.fs.s12 = Arc(source=m.fs.hx_brine.hot_outlet, destination=m.fs.brine.inlet)
    m.fs.s13 = Arc(source=m.fs.condenser.outlet, destination=m.fs.tb_distillate.inlet)
    m.fs.s14 = Arc(source=m.fs.tb_distillate.outlet, destination=m.fs.pump_distillate.inlet)
    m.fs.s15 = Arc(
        source=m.fs.pump_distillate.outlet, destination=m.fs.hx_distillate.hot_inlet
    )
    m.fs.s16 = Arc(
        source=m.fs.hx_distillate.hot_outlet, destination=m.fs.distillate.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.evaporator.connect_to_condenser(m.fs.condenser)

    # Add costing
    add_costing(m)

    # Add recovery ratio
    m.fs.recovery = Var(m.fs.config.time,initialize=0.5,bounds=(0,1))
    m.fs.recovery_equation = Constraint(expr= m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] ==
                                              m.fs.recovery[0]*(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"] +
                                              m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]))

    # Make split ratio equal to recovery
    # m.fs.split_ratio_recovery_equality = Constraint(expr=m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"] ==
    #                                                 m.fs.recovery[0])

    # Scaling
    # properties
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    # unit model values
    # pumps
    iscale.set_scaling_factor(m.fs.pump_feed.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_brine.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_distillate.control_volume.work, 1e-3)

    # distillate HX
    iscale.set_scaling_factor(m.fs.hx_distillate.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_distillate.cold.heat, 1e-3)
    iscale.set_scaling_factor(
        m.fs.hx_distillate.overall_heat_transfer_coefficient, 1e-3
    )
    iscale.set_scaling_factor(m.fs.hx_distillate.area, 1e-1)
    iscale.constraint_scaling_transform(m.fs.hx_distillate.cold.pressure_drop, 1e-5)
    iscale.constraint_scaling_transform(m.fs.hx_distillate.hot.pressure_drop, 1e-5)

    # evaporator
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)

    # compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)

    # condenser
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-6)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def add_Q_ext(m,time_point=None):
    if time_point is None:
        time_point = m.fs.config.time
    m.fs.Q_ext = Var(time_point, initialize=0, units=pyunits.J/pyunits.s)#, bounds=(0,1e7))
    m.fs.Q_ext[0].setlb(0)
    m.fs.evaporator.eq_energy_balance.deactivate()
    m.fs.evaporator.eq_energy_balance_with_additional_Q = Constraint(expr=
        m.fs.evaporator.heat_transfer + m.fs.Q_ext[0] + m.fs.evaporator.properties_feed[0].enth_flow == m.fs.evaporator.properties_brine[0].enth_flow
        + m.fs.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
    )
    iscale.set_scaling_factor(m.fs.Q_ext, 1e-6)

def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.pump_feed.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.pump_distillate.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.pump_brine.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.hx_distillate.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.evaporator.costing = UnitModelCostingBlock(
        default={'flowsheet_costing_block': m.fs.costing}
    )
    m.fs.compressor.costing = UnitModelCostingBlock(
        default={'flowsheet_costing_block': m.fs.costing}
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.distillate.properties[0].flow_vol)

    # Add costing expressions
    m.fs.costing.MVC_comp = Set(initialize=["feed_pump",
                                            "distillate_pump",
                                            "brine_pump",
                                            "hx_distillate",
                                            "evaporator",
                                            "compressor"])

    # Percentage of capital costs
    m.fs.costing.MVC_captial_cost_percentage = Expression(m.fs.costing.MVC_comp)
    m.fs.costing.MVC_captial_cost_percentage["feed_pump"] = (m.fs.pump_feed.costing.capital_cost/m.fs.costing.aggregate_capital_cost)
    m.fs.costing.MVC_captial_cost_percentage["distillate_pump"] = (m.fs.pump_distillate.costing.capital_cost/m.fs.costing.aggregate_capital_cost)
    m.fs.costing.MVC_captial_cost_percentage["brine_pump"] = (m.fs.pump_brine.costing.capital_cost/m.fs.costing.aggregate_capital_cost)
    m.fs.costing.MVC_captial_cost_percentage["hx_distillate"] = (m.fs.hx_distillate.costing.capital_cost/m.fs.costing.aggregate_capital_cost)
    m.fs.costing.MVC_captial_cost_percentage["evaporator"] = (m.fs.evaporator.costing.capital_cost/m.fs.costing.aggregate_capital_cost)
    m.fs.costing.MVC_captial_cost_percentage["compressor"] = (m.fs.compressor.costing.capital_cost/m.fs.costing.aggregate_capital_cost)

    # Percentage of costs normalized to LCOW
    m.fs.costing.annual_operating_costs = Expression(expr=m.fs.costing.total_investment_cost*m.fs.costing.factor_capital_annualization+m.fs.costing.total_operating_cost)
    m.fs.costing.MVC_LCOW_comp = Set(initialize=["feed_pump",
                                            "distillate_pump",
                                            "brine_pump",
                                            "hx_distillate",
                                            "evaporator",
                                            "compressor",
                                            "electricity",
                                            "MLC",
                                            "capital_costs",
                                            "operating_costs"])
    m.fs.costing.LCOW_percentage = Expression(m.fs.costing.MVC_LCOW_comp)
    m.fs.costing.LCOW_percentage["feed_pump"] = (m.fs.pump_feed.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["distillate_pump"] = (m.fs.pump_distillate.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["brine_pump"] = (m.fs.pump_brine.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["hx_distillate"] = (m.fs.hx_distillate.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["evaporator"] = (m.fs.evaporator.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["compressor"] = (m.fs.compressor.costing.capital_cost*m.fs.costing.factor_total_investment*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage['electricity'] = (m.fs.costing.aggregate_flow_costs['electricity']*m.fs.costing.utilization_factor/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage['MLC'] = (m.fs.costing.maintenance_labor_chemical_operating_cost/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["capital_costs"] = (m.fs.costing.total_investment_cost*m.fs.costing.factor_capital_annualization/m.fs.costing.annual_operating_costs)
    m.fs.costing.LCOW_percentage["operating_costs"] = (m.fs.costing.total_operating_cost/m.fs.costing.annual_operating_costs)

def add_pressure_drop_to_hx(hx_blk, time_point):
    # input: hx_blk - heat exchanger block
    # output: deactivates control volume pressure balance and adds pressure drop to pressure balance equation
    hx_blk.cold.pressure_balance.deactivate()
    hx_blk.hot.pressure_balance.deactivate()
    hx_blk.cold.deltaP = Var(time_point, initialize=7e4, units=pyunits.Pa)
    hx_blk.hot.deltaP = Var(time_point, initialize=7e4, units=pyunits.Pa)
    hx_blk.cold.pressure_drop = Constraint(
        expr=hx_blk.cold.properties_in[0].pressure
        == hx_blk.cold.properties_out[0].pressure + hx_blk.cold.deltaP[0]
    )
    hx_blk.hot.pressure_drop = Constraint(
        expr=hx_blk.hot.properties_in[0].pressure
        == hx_blk.hot.properties_out[0].pressure + hx_blk.hot.deltaP[0]
    )

def set_operating_conditions(m):
    # Feed inlet
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].fix(0.1)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40)
    #m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(2)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.recovery[0].fix(0.5)

    # Feed pump
    m.fs.pump_feed.efficiency_pump.fix(0.8)
    m.fs.pump_feed.control_volume.deltaP[0].fix(7e3)

    # Separator
    # m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"] = m.fs.recovery[0].value

    # distillate HX
    m.fs.hx_distillate.overall_heat_transfer_coefficient.fix(1e3)
    m.fs.hx_distillate.area.fix(200)
    m.fs.hx_distillate.cold.deltaP[0].fix(7e3)
    m.fs.hx_distillate.hot.deltaP[0].fix(7e3)

    # evaporator specifications
    #m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(1000)  # m^2
    #m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] = 5
    #m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(20)

    # compressor
    m.fs.compressor.pressure_ratio.fix(2)
    #m.fs.compressor.control_volume.properties_out[0].temperature = 400
    # m.fs.compressor.pressure_ratio = 2
    m.fs.compressor.efficiency.fix(0.8)

    # Brine pump
    m.fs.pump_brine.efficiency_pump.fix(0.8)
    m.fs.pump_brine.control_volume.deltaP[0].fix(4e4)

    # Distillate pump
    m.fs.pump_distillate.efficiency_pump.fix(0.8)
    m.fs.pump_distillate.control_volume.deltaP[0].fix(4e4)

    # Fix 0 TDS
    m.fs.tb_distillate.properties_out[0].flow_mass_phase_comp['Liq','TDS'].fix(1e-5)

    # Costing
    m.fs.costing.factor_total_investment.fix(2)
    m.fs.costing.electricity_base_cost = 0.15
    # print(value(m.fs.costing.electricity_base_cost))
    m.fs.costing.heat_exchanger_unit_cost.fix(2000)
    m.fs.costing.evaporator_unit_cost.fix(3000)

    # Change upper bound of compressed vapor temperature
    # m.fs.compressor.control_volume.properties_out[0].temperature.setub(450)
    # check degrees of freedom
    print("DOF after setting operating conditions: ", degrees_of_freedom(m))

def set_operating_conditions_validation(m):
    set_operating_conditions(m)
    m_f = 40 # kg/s
    w_f = 50*1e-3 # kg/kg
    rr = 0.5
    # Feed inlet
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix((1-w_f)*m_f)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(w_f*m_f)

    # Feed pump
    m.fs.pump_feed.control_volume.deltaP[0].fix(2e3)

    # Separator
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].fix(0.5)

    # distillate HX
    m.fs.hx_distillate.overall_heat_transfer_coefficient.fix(1e3)
    m.fs.hx_distillate.area.fix(50)
    m.fs.hx_distillate.cold.deltaP[0].fix(2e3)
    m.fs.hx_distillate.hot.deltaP[0].fix(2e3)

    # brine HX
    m.fs.hx_brine.overall_heat_transfer_coefficient.fix(1e3)
    m.fs.hx_brine.area.fix(50)  # = m.fs.hx_distillate.area.value
    m.fs.hx_brine.cold.deltaP[0].fix(2e3)
    m.fs.hx_brine.hot.deltaP[0].fix(2e3)

    # evaporator specifications
    m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area = 2000  # m^2
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(rr*m_f)

    # compressor
    m.fs.compressor.pressure_ratio.fix(2)
    m.fs.compressor.efficiency.fix(0.8)

    # Fix 0 TDS
    m.fs.tb_distillate.properties_out[0].flow_mass_phase_comp['Liq','TDS'].fix(1e-5)

    # Adjust scaling for higher feed flow rate
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    # compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-8)

    # condenser
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-8)

    iscale.calculate_scaling_factors(m)

    # Costing
    m.fs.costing.factor_total_investment.fix(0.5)

    # check degrees of freedom
    print("DOF after setting operating conditions: ", degrees_of_freedom(m))

def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # Touch feed mass fraction property
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS']

    # Propagate vapor flow rate
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] = m.fs.recovery[0] * (m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"] + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"])

    # initialize feed pump
    propagate_state(m.fs.s01)
    # m.fs.pump_feed.control_volume.deltaP[0].fix(7e4)
    m.fs.pump_feed.initialize(optarg=optarg)

    # initialize separator
    # propagate_state(m.fs.s02)
    # Touch property for initialization
    # m.fs.separator_feed.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    # m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].fix(m.fs.recovery[0].value)

    # m.fs.separator_feed.mixed_state.initialize(optarg=optarg)
    # Touch properties for initialization
    # m.fs.separator_feed.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    # m.fs.separator_feed.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    # m.fs.separator_feed.initialize(optarg=optarg, outlvl=idaeslog.INFO_HIGH)
    # m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].unfix()

    # initialize distillate heat exchanger
    propagate_state(m.fs.s03)
    m.fs.hx_distillate.cold_outlet.temperature[
        0
    ] = m.fs.evaporator.inlet_feed.temperature[0].value
    m.fs.hx_distillate.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[
        0
    ].value
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    m.fs.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
    m.fs.hx_distillate.hot_inlet.temperature[
        0
    ] = m.fs.evaporator.outlet_brine.temperature[0].value
    m.fs.hx_distillate.hot_inlet.pressure[0] = 101325

    # fix for initialization
    m.fs.hx_distillate.initialize()

    # initialize brine heat exchanger
    # propagate_state(m.fs.s04)
    # m.fs.hx_brine.cold_outlet.temperature[0] = m.fs.evaporator.inlet_feed.temperature[
    #     0
    # ].value
    # m.fs.hx_brine.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[0].value
    # m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
    #     m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    #     - m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    # )
    # m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = (
    #     m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].value
    # )
    # m.fs.hx_brine.hot_inlet.temperature[0] = m.fs.evaporator.outlet_brine.temperature[
    #     0
    # ].value
    # m.fs.hx_brine.hot_inlet.pressure[0] = 101325
    # m.fs.hx_brine.initialize()

    # initialize mixer
    # propagate_state(m.fs.s05)
    # propagate_state(m.fs.s06)
    # m.fs.mixer_feed.initialize()
    # m.fs.mixer_feed.pressure_equality_constraints[0,2].deactivate()

    # initialize evaporator
    propagate_state(m.fs.s07)
    m.fs.Q_ext[0].fix()
    m.fs.evaporator.initialize(
        delta_temperature_in=70, delta_temperature_out=9
    )  # fixes and unfixes those values
    m.fs.Q_ext[0].unfix()

    # initialize compressor
    propagate_state(m.fs.s08)
    m.fs.compressor.initialize_build()

    # initialize condenser
    propagate_state(m.fs.s09)
    m.fs.condenser.initialize(heat=-m.fs.evaporator.heat_transfer.value)

    # initialize brine pump
    propagate_state(m.fs.s10)
    m.fs.pump_brine.initialize(optarg=optarg)

    # initialize distillate pump
    propagate_state(m.fs.s13) # to translator block
    propagate_state(m.fs.s14) # from translator block to pump
    m.fs.pump_distillate.initialize(optarg=optarg)

    m.fs.costing.initialize()

    print('Initialization done')

def fix_outlet_pressures(m):
    # unfix pump heads
    m.fs.pump_brine.control_volume.deltaP[0].unfix()
    # m.fs.pump_distillate.control_volume.deltaP[0].unfix()

    # Fix outlet pressures
    m.fs.brine.properties[0].pressure.fix(101325)
    # m.fs.distillate.properties[0].pressure.fix(101325)

    return

def calculate_cost_sf(cost):
    sf = 10**-(math.log10(abs(cost.value)))
    iscale.set_scaling_factor(cost, sf)

def scale_costs(m):
    calculate_cost_sf(m.fs.hx_distillate.costing.capital_cost)
    calculate_cost_sf(m.fs.evaporator.costing.capital_cost)
    calculate_cost_sf(m.fs.compressor.costing.capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_flow_costs['electricity'])
    calculate_cost_sf(m.fs.costing.total_investment_cost)
    calculate_cost_sf(m.fs.costing.maintenance_labor_chemical_operating_cost)
    calculate_cost_sf(m.fs.costing.total_operating_cost)

    iscale.calculate_scaling_factors(m)

    print('Scaled costs')

def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    results = solver.solve(model, tee=tee)
    if check_optimal_termination(results):
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        print(msg)
        return results

def sweep_solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    # First simulate minimizing Q_ext
    model.fs.objective = Objective(expr=model.fs.Q_ext[0])
    results = solver.solve(model, tee=tee)

    # Now optimize
    model.fs.Q_ext[0].fix(0)
    del model.fs.objective
    set_up_optimization(model)
    results = solve(model)

    if check_optimal_termination(results):
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        print(msg)
        return results

def sweep_solve_fixed_brine_temp(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    # First simulate minimizing Q_ext
    model.fs.objective = Objective(expr=model.fs.Q_ext[0])
    results = solver.solve(model, tee=tee)

    # Now optimize
    model.fs.Q_ext[0].fix(0)
    del model.fs.objective
    set_up_optimization_fixed_brine_temp(model)
    results = solve(model)

    if check_optimal_termination(results):
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        print(msg)
        return results

def debug_infeasible(m, solver):
    print("\n---infeasible constraints---")
    infeas.log_infeasible_constraints(m)
    print("\n---infeasible bounds---")
    infeas.log_infeasible_bounds(m)
    print("\n---close to bounds---")
    infeas.log_close_to_bounds(m)
    print("\n---poor scaling---")
    bsv_gen = iscale.badly_scaled_var_generator(m)
    # for output in bsv_gen:
    #     var = output[0]
    #     val = output[1]
    #     print(var.name, val)
    for (var, val) in bsv_gen:
        print(var.name, val)
    # Create Degeneracy Hunter object
    print("\n---degeneracy hunter---")
    dh = DegeneracyHunter(m, solver=solver)
    dh.check_residuals(tol=1e-8)
    # dh.find_candidate_equations(verbose=True,tee=True)
    # dh.check_rank_equality_constraints()

def set_up_optimization_fixed_brine_temp(m):
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    m.fs.Q_ext[0].fix(0)
    m.fs.evaporator.area.unfix()
    # m.fs.evaporator.outlet_brine.temperature[0].unfix()
    m.fs.compressor.pressure_ratio.unfix()
    m.fs.hx_distillate.area.unfix()
    m.fs.hx_brine.area.unfix()
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].unfix()

    print("DOF for optimization: ", degrees_of_freedom(m))

def set_up_optimization(m):
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    m.fs.Q_ext[0].fix(0)
    m.fs.evaporator.area.unfix()
    m.fs.evaporator.outlet_brine.temperature[0].unfix()
    m.fs.compressor.pressure_ratio.unfix()
    m.fs.hx_distillate.area.unfix()

    print("DOF for optimization: ", degrees_of_freedom(m))

def display_results(m):
    print("Feed flow rate:                          ", m.fs.feed.properties[0].flow_mass_phase_comp['Liq','H2O'].value+
          m.fs.feed.properties[0].flow_mass_phase_comp['Liq','TDS'].value)
    print("Feed mass fraction:                      ", m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].value)
    print("Brine salinity: ", m.fs.brine.properties[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3," g/kg")
    print("Vapor flow rate:                         ", m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value)
    print('Recovery:                                ', m.fs.recovery[0].value)
    print('Distillate preheater area:               ', m.fs.hx_distillate.area.value)
    print("Preheated feed temperature:              ", m.fs.evaporator.properties_feed[0].temperature.value)
    print("Evaporator temperature:                  ", m.fs.evaporator.properties_brine[0].temperature.value)
    print("Compressed vapor temperature:            ", m.fs.compressor.control_volume.properties_out[0].temperature.value)
    print("Condensed vapor temperature:             ", m.fs.condenser.control_volume.properties_out[0].temperature.value)
    print("Compressor work:                         ", m.fs.compressor.control_volume.work[0].value)
    print("Compressor pressure ratio:               ", m.fs.compressor.pressure_ratio.value)
    print("Evaporator area:                         ", m.fs.evaporator.area.value)
    print('Evaporator LMTD:                         ', m.fs.evaporator.lmtd.value)
    print('Specific energy consumption:             ', value(m.fs.costing.specific_energy_consumption))
    print('Total investment factor:                 ', m.fs.costing.factor_total_investment.value)
    print('LCOW:                                    ', m.fs.costing.LCOW.value)
    print('External Q:                              ', m.fs.Q_ext[0].value)

def save_results(m,filename=None):

    feed_concentration = round(m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].value,3)*1000
    recovery = round(m.fs.recovery[0].value,2)*100

    if filename is None:
        filename = "C:/Users/carso/Documents/MVC/watertap_results/wf_" + str(feed_concentration) + "_rr_" + str(recovery) + ".csv"

    outputs = {}
    outputs['Feed concentration'] = [feed_concentration]
    outputs['Recovery'] = [recovery]
    outputs['Feed mass flow water'] = [m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value]
    outputs['Feed mass flow salt'] = [m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'TDS'].value]
    outputs['Feed mass fraction'] = [m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].value]
    outputs['Feed temperature'] = [m.fs.feed.properties[0].temperature.value]
    outputs['Feed pressure'] = [m.fs.feed.properties[0].pressure.value]

    # Brine from evaporator
    outputs['Brine mass flow water'] = [m.fs.brine.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value]
    outputs['Brine mass flow salt'] = [m.fs.brine.properties[0].flow_mass_phase_comp['Liq', 'TDS'].value]
    outputs['Brine temperature'] = [m.fs.evaporator.properties_brine[0].temperature.value]
    outputs['Brine pressure'] = [m.fs.evaporator.properties_brine[0].pressure.value]

    # Vapor
    outputs['Vapor mass flow'] = [m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O'].value]
    outputs['Vapor temperature'] = [m.fs.evaporator.properties_vapor[0].temperature.value]
    outputs['Vapor pressure'] = [m.fs.evaporator.properties_vapor[0].pressure.value]

    # Compressed vapor
    outputs['Compressed vapor temperature'] = [m.fs.compressor.control_volume.properties_out[0].temperature.value]
    outputs['Compressed vapor pressure'] = [m.fs.compressor.control_volume.properties_out[0].pressure.value]

    # Condensed vapor/distillate
    outputs['Distillate temperature'] = [m.fs.condenser.control_volume.properties_out[0].temperature.value]
    outputs['Distillate pressure'] = [m.fs.condenser.control_volume.properties_out[0].pressure.value]

    # Exiting distillate
    outputs['Exiting distillate pressure'] = [m.fs.distillate.properties[0].pressure.value]

    # Exiting brine
    outputs['Exiting brine pressure'] = [m.fs.brine.properties[0].pressure.value]

    # Evaporator performance
    outputs['Evaporator area'] = [m.fs.evaporator.area.value]
    outputs['Evaporator LMTD'] = [m.fs.evaporator.lmtd.value]
    outputs['Evaporator heat transfer'] = [m.fs.evaporator.heat_transfer.value]
    outputs['Evaporator overall heat transfer coefficient'] = [m.fs.evaporator.U.value]
    outputs['Evaporator approach temperature in'] = [m.fs.evaporator.delta_temperature_in.value]
    outputs['Evaporator approach temperature out'] = [m.fs.evaporator.delta_temperature_out.value]

    # Compressor performance
    outputs['Compressor pressure ratio'] = [m.fs.compressor.pressure_ratio.value]
    outputs['Compressor work'] = [m.fs.compressor.control_volume.work[0].value]
    outputs['Compressor efficiency'] = [m.fs.compressor.efficiency.value]

    # Preheater performance
    outputs['Preheater split ratio'] = [m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].value]

    # Feed exiting distillate heat exchanger
    outputs['Feed exiting distillate hx temperature'] = [m.fs.hx_distillate.cold.properties_out[0].temperature.value]
    outputs['Feed exiting distillate hx pressure'] = [m.fs.hx_distillate.cold.properties_out[0].pressure.value]

    # Feed exiting brine heat exchanger
    outputs['Feed exiting brine hx temperature'] = [m.fs.hx_brine.cold.properties_out[0].temperature.value]
    outputs['Feed exiting brine hx pressure'] = [m.fs.hx_brine.cold.properties_out[0].pressure.value]

    # Preheated feed
    outputs['Preheated feed temperature'] = [m.fs.evaporator.properties_feed[0].temperature.value]
    outputs['Preheated feed pressure'] = [m.fs.evaporator.properties_feed[0].pressure.value]

    # Distillate heat exchanger performance
    outputs['Distillate hx area'] = [m.fs.hx_distillate.area.value]
    outputs['Distillate hx delta temp in'] = [m.fs.hx_distillate.delta_temperature_in[0].value]
    outputs['Distillate hx delta temp out'] = [m.fs.hx_distillate.delta_temperature_out[0].value]
    outputs['Distillate hx heat transfer'] = [m.fs.hx_distillate.heat_duty[0].value]
    outputs['Distillate hx overall heat transfer coefficient'] = [m.fs.hx_distillate.overall_heat_transfer_coefficient[0].value]

    # Brine heat exchanger performance
    outputs['Brine hx area'] = [m.fs.hx_brine.area.value]
    outputs['Brine hx delta temp in'] = [m.fs.hx_brine.delta_temperature_in[0].value]
    outputs['Brine hx delta temp out'] = [m.fs.hx_brine.delta_temperature_out[0].value]
    outputs['Brine hx heat transfer'] = [m.fs.hx_brine.heat_duty[0].value]
    outputs['Brine hx overall heat transfer coefficient'] = [m.fs.hx_brine.overall_heat_transfer_coefficient[0].value]

    # External Q
    outputs['Q external'] = [m.fs.Q_ext[0].value]

    # Cost outcome metrics
    outputs['Feed pump capital cost'] = [m.fs.pump_feed.costing.capital_cost.value]
    outputs['Distillate pump captial cost'] = [m.fs.pump_distillate.costing.capital_cost.value]
    outputs['Brine pump captial cost'] = [m.fs.pump_distillate.costing.capital_cost.value]
    outputs['Distillate hx capital cost'] = [m.fs.hx_distillate.costing.capital_cost.value]
    outputs['Brine hx capital cost'] = [m.fs.hx_brine.costing.capital_cost.value]
    outputs['Mixer capital cost'] = [m.fs.mixer_feed.costing.capital_cost.value]
    outputs['Evaporator capital cost'] = [m.fs.evaporator.costing.capital_cost.value]
    outputs['Compressor capital cost'] = [m.fs.compressor.costing.capital_cost.value]
    outputs['Aggregate capital cost'] = [m.fs.costing.aggregate_capital_cost.value]
    outputs['Aggregate electricity flow cost'] = [value(m.fs.costing.aggregate_flow_costs['electricity'])]
    outputs['Total investment cost'] = [m.fs.costing.total_investment_cost.value]
    outputs['Total investment factor'] = [m.fs.costing.factor_total_investment.value]
    outputs['Evaporator material factor'] = [m.fs.costing.evaporator_material_factor_cost.value]
    outputs['Distillate hx material factor'] = [m.fs.costing.heat_exchanger_material_factor_cost.value]
    outputs['Brine hx material factor'] = [m.fs.costing.heat_exchanger_material_factor_cost.value]
    outputs['Electricity cost'] = [m.fs.costing.electricity_base_cost.value]
    outputs['Capital annualization factor'] = [m.fs.costing.factor_capital_annualization.value]
    outputs['Utilization factor'] = [m.fs.costing.utilization_factor.value]
    outputs['MLC cost'] = [value(m.fs.costing.maintenance_labor_chemical_operating_cost)]
    outputs['Total operating cost'] = [value(m.fs.costing.total_operating_cost)]
    outputs['LCOW'] = [m.fs.costing.LCOW.value]
    outputs['SEC'] = [value(m.fs.costing.specific_energy_consumption)]

    # Normalized capital costs
    outputs['CC normalized feed pump'] = [value(m.fs.costing.MVC_captial_cost_percentage['feed_pump'])]
    outputs['CC normalized distillate pump'] = [value(m.fs.costing.MVC_captial_cost_percentage["distillate_pump"])]
    outputs['CC normalized brine pump'] = [value(m.fs.costing.MVC_captial_cost_percentage["brine_pump"])]
    outputs['CC normalized distiallte hx'] = [value(m.fs.costing.MVC_captial_cost_percentage["hx_distillate"])]
    outputs['CC normalized brine hx'] = [value(m.fs.costing.MVC_captial_cost_percentage["hx_brine"])]
    outputs['CC normalized mixer'] =[value(m.fs.costing.MVC_captial_cost_percentage["mixer"])]
    outputs['CC normalized evaportor'] =[value(m.fs.costing.MVC_captial_cost_percentage["evaporator"])]
    outputs['CC normalized compressor'] =[value(m.fs.costing.MVC_captial_cost_percentage["compressor"])]

    # Normalized LCOW costs
    outputs['LCOW normalized feed pump'] = [value(m.fs.costing.LCOW_percentage["feed_pump"])]
    outputs['LCOW normalized distillate pump'] = [value(m.fs.costing.LCOW_percentage["distillate_pump"])]
    outputs['LCOW normalized brine pump'] = [value(m.fs.costing.LCOW_percentage["brine_pump"])]
    outputs['LCOW normalized distillate hx'] = [value(m.fs.costing.LCOW_percentage["hx_distillate"])]
    outputs['LCOW normalized brine hx'] = [value(m.fs.costing.LCOW_percentage["hx_brine"])]
    outputs['LCOW normalized mixer'] = [value(m.fs.costing.LCOW_percentage["mixer"])]
    outputs['LCOW normalized evaporator'] = [value(m.fs.costing.LCOW_percentage["evaporator"])]
    outputs['LCOW normalized compressor'] = [value(m.fs.costing.LCOW_percentage["compressor"])]
    outputs['LCOW normalized electricity'] = [value(m.fs.costing.LCOW_percentage['electricity'])]
    outputs['LCOW normalized MLC'] = [value(m.fs.costing.LCOW_percentage['MLC'])]
    outputs['LCOW normalized capex'] = [value(m.fs.costing.LCOW_percentage["capital_costs"])]
    outputs['LCOW normalized opex'] = [value(m.fs.costing.LCOW_percentage["operating_costs"])]

    df = pd.DataFrame(outputs)
    df = df.set_index('Feed concentration').T

    df.to_csv(filename,index=False)

def display_system(m):
    recovery = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ].value / (
        m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
        + m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    print(
        "Feed salinity: ",
        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3,
        " g/kg",
    )
    print(
        "Brine salinity: ",
        m.fs.brine.properties[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3,
        " g/kg",
    )
    print("Recovery: ", recovery)
    # print('\nSplitter')
    print("\nDistillate heat exchanger")
    print("Area: ", m.fs.hx_distillate.area.value, " m^2")
    print(
        "U: ", m.fs.hx_distillate.overall_heat_transfer_coefficient[0].value, " W/m^2-K"
    )
    print("Heat transfer: ", m.fs.hx_distillate.heat_duty[0].value, " W")
    print("\nBrine heat exchanger")
    print("\nMixed feed")
    print("\nEvaporator")
    print(
        "Temperature: ", m.fs.evaporator.properties_brine[0].temperature.value, " K"
    )  # , ', Fixed? ',  m.fs.evaporator.outlet_brine.temperature[0].fixed())
    print("Pressure: ", m.fs.evaporator.properties_brine[0].pressure.value, " Pa")
    print(
        "Area: ", m.fs.evaporator.area.value, " m^2"
    )  # , ', Fixed? ', m.fs.evaporator.area.isfixed())
    print(
        "U: ", m.fs.evaporator.U.value, " W/m^2-K"
    )  # , ', Fixed? ', m.fs.evaporator.U.isfixed())
    print("heat transfer: ", m.fs.evaporator.heat_transfer.value, " W")
    print("\nCompressor")
    print("Work: ", m.fs.compressor.control_volume.work[0].value, " W")
    print("Pressure ratio: ", m.fs.compressor.pressure_ratio.value)
    print("Efficiency: ", m.fs.compressor.efficiency.value)
    print("\nCondenser")
    print("Heat transfer: ", m.fs.condenser.control_volume.heat[0].value, " W")

def display_seawater_states(state_blk):
    print("water mass flow ", state_blk.flow_mass_phase_comp["Liq", "H2O"].value)
    print("TDS mass flow   ", state_blk.flow_mass_phase_comp["Liq", "TDS"].value)
    print("temperature     ", state_blk.temperature.value)
    print("pressure        ", state_blk.pressure.value)
    print("specific enthalpy", value(state_blk.enth_flow))

def display_water_states(state_blk):
    print("Liquid mass flow ", state_blk.flow_mass_phase_comp["Liq", "H2O"].value)
    print("Vapor mass flow  ", state_blk.flow_mass_phase_comp["Vap", "H2O"].value)
    print("temperature      ", state_blk.temperature.value)
    print("pressure         ", state_blk.pressure.value)
    print("specific enthalpy vapor", state_blk.enth_flow_phase['Vap'].value)
    print("specific enthalpy liquid", state_blk.enth_flow_phase['Liq'].value)


if __name__ == "__main__":
    m = main()