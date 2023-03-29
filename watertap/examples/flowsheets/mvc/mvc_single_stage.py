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
    Constraint,
    Objective,
    Var,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
    check_optimal_termination
)
from pyomo.network import Arc

import pyomo.environ as pyo
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
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.costing import WaterTAPCosting
import math
import pandas as pd

def main():
    # build, set operating conditions, initialize for simulation
    m = build()
    set_operating_conditions(m)
    add_Q_ext(m, time_point=m.fs.config.time)
    initialize_system(m)
    scale_costs(m) # rescale costs after initialization because scaling depends on flow rates
    fix_outlet_pressures(m) # outlet pressure are initially unfixed for initialization

    # set up for minimizing Q_ext in first solve
    print('DOF after initialization: ', degrees_of_freedom(m)) # should be 1 because Q_ext is unfixed
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])

    print('\nFirst solve - simulation results')
    solver = get_solver()
    results = solve(m, solver=solver, tee=False)
    print(results.solver.termination_condition)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)
    display_results(m)

    print('\nSecond solve - optimization results')
    m.fs.Q_ext[0].fix(0) # no longer want external heating in evaporator
    del m.fs.objective
    set_up_optimization(m)
    results = solve(m, tee=False)
    print(results.solver.termination_condition)
    display_results(m)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)

    assert_units_consistent(m)


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Unit models
    m.fs.feed = Feed(property_package=m.fs.properties_feed)

    m.fs.pump_feed = Pump(property_package=m.fs.properties_feed)

    m.fs.separator_feed = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["hx_distillate_cold", "hx_brine_cold"],
        split_basis=SplittingType.totalFlow
    )

    m.fs.hx_distillate = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed},
        cold={"property_package": m.fs.properties_feed},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent
    )

    # Set lower bound of approach temperatures
    m.fs.hx_distillate.delta_temperature_in.setlb(0)
    m.fs.hx_distillate.delta_temperature_out.setlb(0)
    m.fs.hx_distillate.area.setlb(10)
    add_pressure_drop_to_hx(m.fs.hx_distillate, m.fs.config.time)

    m.fs.hx_brine = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed},
        cold={"property_package": m.fs.properties_feed},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent
    )
    # Set lower bound of approach temperatures
    m.fs.hx_brine.delta_temperature_in.setlb(0)
    m.fs.hx_brine.delta_temperature_out.setlb(0)
    m.fs.hx_brine.area.setlb(10)
    add_pressure_drop_to_hx(m.fs.hx_brine, m.fs.config.time)

    m.fs.mixer_feed = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["hx_distillate_cold", "hx_brine_cold"]
    )
    m.fs.mixer_feed.pressure_equality_constraints[0,2].deactivate()

    m.fs.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )

    m.fs.compressor = Compressor(property_package=m.fs.properties_vapor)

    m.fs.condenser = Condenser(property_package=m.fs.properties_vapor)

    m.fs.tb_distillate = Translator(
        inlet_property_package=m.fs.properties_vapor,
        outlet_property_package=m.fs.properties_feed
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

    m.fs.pump_brine = Pump(property_package=m.fs.properties_feed)

    m.fs.pump_distillate = Pump(property_package=m.fs.properties_feed)

    m.fs.distillate = Product(property_package=m.fs.properties_feed)

    m.fs.brine = Product(property_package=m.fs.properties_feed)

    # Connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_feed.inlet)
    m.fs.s02 = Arc(source=m.fs.pump_feed.outlet, destination=m.fs.separator_feed.inlet)
    m.fs.s03 = Arc(
        source=m.fs.separator_feed.hx_distillate_cold,
        destination=m.fs.hx_distillate.cold_inlet,
    )
    m.fs.s04 = Arc(
        source=m.fs.separator_feed.hx_brine_cold, destination=m.fs.hx_brine.cold_inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.hx_distillate.cold_outlet,
        destination=m.fs.mixer_feed.hx_distillate_cold,
    )
    m.fs.s06 = Arc(
        source=m.fs.hx_brine.cold_outlet, destination=m.fs.mixer_feed.hx_brine_cold
    )
    m.fs.s07 = Arc(
        source=m.fs.mixer_feed.outlet, destination=m.fs.evaporator.inlet_feed
    )
    m.fs.s08 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s09 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    m.fs.s10 = Arc(
        source=m.fs.evaporator.outlet_brine, destination=m.fs.pump_brine.inlet
    )
    m.fs.s11 = Arc(source=m.fs.pump_brine.outlet, destination=m.fs.hx_brine.hot_inlet)
    m.fs.s12 = Arc(source=m.fs.hx_brine.hot_outlet, destination=m.fs.brine.inlet)
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
    m.fs.split_ratio_recovery_equality = Constraint(expr=m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"] ==
                                                    m.fs.recovery[0])

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

    # brine HX
    iscale.set_scaling_factor(m.fs.hx_brine.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.cold.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_brine.area, 1e-1)
    iscale.constraint_scaling_transform(m.fs.hx_brine.cold.pressure_drop, 1e-5)
    iscale.constraint_scaling_transform(m.fs.hx_brine.hot.pressure_drop, 1e-5)

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
    m.fs.Q_ext = Var(time_point, initialize=0, units=pyunits.J/pyunits.s)
    m.fs.Q_ext[0].setlb(0)
    m.fs.evaporator.eq_energy_balance.deactivate()
    m.fs.evaporator.eq_energy_balance_with_additional_Q = Constraint(expr=
        m.fs.evaporator.heat_transfer + m.fs.Q_ext[0] + m.fs.evaporator.properties_feed[0].enth_flow == m.fs.evaporator.properties_brine[0].enth_flow
        + m.fs.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
    )
    iscale.set_scaling_factor(m.fs.Q_ext, 1e-6)


def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.pump_feed.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.pump_distillate.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.pump_brine.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.hx_distillate.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.hx_brine.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.mixer_feed.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.evaporator.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.compressor.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.distillate.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.distillate.properties[0].flow_vol)


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
    m.fs.costing.base_currency = pyo.units.USD_2020

    # Feed inlet
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].fix(0.1)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(4)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.recovery[0].fix(0.5)

    # Feed pump
    m.fs.pump_feed.efficiency_pump.fix(0.8)
    m.fs.pump_feed.control_volume.deltaP[0].fix(7e3)

    # Separator
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"] = m.fs.recovery[0].value

    # Distillate HX
    m.fs.hx_distillate.overall_heat_transfer_coefficient.fix(2e3)
    m.fs.hx_distillate.area.fix(125)
    m.fs.hx_distillate.cold.deltaP[0].fix(7e3)
    m.fs.hx_distillate.hot.deltaP[0].fix(7e3)

    # Brine HX
    m.fs.hx_brine.overall_heat_transfer_coefficient.fix(2e3)
    m.fs.hx_brine.area.fix(115)
    m.fs.hx_brine.cold.deltaP[0].fix(7e3)
    m.fs.hx_brine.hot.deltaP[0].fix(7e3)

    # Evaporator
    m.fs.evaporator.inlet_feed.temperature[0] = 50+273.15  # provide guess
    m.fs.evaporator.outlet_brine.temperature[0].fix(70+273.15)
    m.fs.evaporator.U.fix(3e3)  # W/K-m^2
    m.fs.evaporator.area.setub(1e4) # m^2

    # Compressor
    m.fs.compressor.pressure_ratio.fix(1.6)
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
    m.fs.costing.electricity_cost = 0.1 # 0.15
    m.fs.costing.heat_exchanger.material_factor_cost.fix(5)
    m.fs.costing.evaporator.material_factor_cost.fix(5)
    m.fs.costing.compressor.unit_cost.fix(1*7364)

    # Temperature bounds
    m.fs.evaporator.properties_vapor[0].temperature.setub(75+273.15)
    # m.fs.evaporator.properties_vapor[0].temperature.setlb(50+273.15)
    m.fs.compressor.control_volume.properties_out[0].temperature.setub(450)

    # check degrees of freedom
    print("DOF after setting operating conditions: ", degrees_of_freedom(m))


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # Touch feed mass fraction property
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq','TDS']
    solver.solve(m.fs.feed)

    # Propagate vapor flow rate
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] = m.fs.recovery[0] * (m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"] + m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"])
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Liq','H2O'] = 0
    # m.fs.evaporator.properties_vapor[0].display()

    # Propagate brine salinity and flow rate
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp['Liq','TDS'] = m.fs.feed.properties[0].mass_frac_phase_comp['Liq','TDS']/(1-m.fs.recovery[0])
    m.fs.evaporator.properties_brine[0].mass_frac_phase_comp['Liq','H2O'] = 1-m.fs.evaporator.properties_brine[0].mass_frac_phase_comp['Liq','TDS'].value
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp['Liq','TDS'] = m.fs.feed.properties[0].flow_mass_phase_comp['Liq','TDS']
    m.fs.evaporator.properties_brine[0].flow_mass_phase_comp['Liq','H2O'] = m.fs.feed.properties[0].flow_mass_phase_comp['Liq','H2O'] - m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap','H2O']
    # m.fs.evaporator.properties_brine[0].display()

    # initialize feed pump
    propagate_state(m.fs.s01)
    m.fs.pump_feed.initialize_build(optarg=optarg)

    # initialize separator
    propagate_state(m.fs.s02)
    # Touch property for initialization
    m.fs.separator_feed.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].fix(m.fs.recovery[0].value)
    m.fs.separator_feed.mixed_state.initialize(optarg=optarg)
    # Touch properties for initialization
    m.fs.separator_feed.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.separator_feed.initialize_build(optarg=optarg, outlvl=idaeslog.INFO_HIGH)
    m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"].unfix()
    # m.fs.separator_feed.report()

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
    m.fs.hx_distillate.initialize_build()

    # initialize brine heat exchanger
    propagate_state(m.fs.s04)
    m.fs.hx_brine.cold_outlet.temperature[0] = m.fs.evaporator.inlet_feed.temperature[
        0
    ].value
    m.fs.hx_brine.cold_outlet.pressure[0] = m.fs.evaporator.inlet_feed.pressure[0].value
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp['Liq','H2O']
    )
    m.fs.hx_brine.hot_inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'] = (
        m.fs.evaporator.properties_brine[0].flow_mass_phase_comp['Liq','TDS']
    )
    m.fs.hx_brine.hot_inlet.temperature[0] = m.fs.evaporator.outlet_brine.temperature[0].value
    m.fs.hx_brine.hot_inlet.pressure[0] = 101325
    m.fs.hx_brine.initialize_build()

    # initialize mixer
    propagate_state(m.fs.s05)
    propagate_state(m.fs.s06)
    m.fs.mixer_feed.initialize_build()
    m.fs.mixer_feed.pressure_equality_constraints[0,2].deactivate()

    # initialize evaporator
    propagate_state(m.fs.s07)
    m.fs.Q_ext[0].fix()
    # m.fs.evaporator.initialize_build(
    #     delta_temperature_in=60, delta_temperature_out=10
    # )  # fixes and unfixes those values
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    m.fs.evaporator.initialize_build(
        delta_temperature_in=60
    )  # fixes and unfixes those values
    m.fs.Q_ext[0].unfix()
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    m.fs.evaporator.display()
    # initialize compressor
    propagate_state(m.fs.s08)
    m.fs.compressor.initialize_build()
    m.fs.compressor.report()

    # initialize condenser
    propagate_state(m.fs.s09)
    m.fs.condenser.initialize_build(heat=-m.fs.evaporator.heat_transfer.value)

    # initialize brine pump
    propagate_state(m.fs.s10)
    m.fs.pump_brine.initialize_build(optarg=optarg)

    # initialize distillate pump
    propagate_state(m.fs.s13) # to translator block
    propagate_state(m.fs.s14) # from translator block to pump
    m.fs.pump_distillate.control_volume.properties_in[0].temperature = m.fs.condenser.control_volume.properties_out[0].temperature.value
    m.fs.pump_distillate.control_volume.properties_in[0].pressure = m.fs.condenser.control_volume.properties_out[0].pressure.value
    m.fs.pump_distillate.initialize_build(optarg=optarg)

    # propagate brine state
    propagate_state(m.fs.s12)
    propagate_state(m.fs.s16)
    #
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
    calculate_cost_sf(m.fs.hx_brine.costing.capital_cost)
    calculate_cost_sf(m.fs.mixer_feed.costing.capital_cost)
    calculate_cost_sf(m.fs.evaporator.costing.capital_cost)
    calculate_cost_sf(m.fs.compressor.costing.capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_capital_cost)
    calculate_cost_sf(m.fs.costing.aggregate_flow_costs['electricity'])
    calculate_cost_sf(m.fs.costing.total_capital_cost)
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

def set_up_optimization(m):
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    m.fs.Q_ext[0].fix(0)
    m.fs.evaporator.area.unfix()
    m.fs.evaporator.outlet_brine.temperature[0].unfix()
    m.fs.compressor.pressure_ratio.unfix()
    m.fs.hx_distillate.area.unfix()
    m.fs.hx_brine.area.unfix()

    print("DOF for optimization: ", degrees_of_freedom(m))

def display_results(m): # TODO: break up display results into operating, design, metrics
    print("Feed flow rate:                          ", m.fs.feed.properties[0].flow_mass_phase_comp['Liq','H2O'].value+
          m.fs.feed.properties[0].flow_mass_phase_comp['Liq','TDS'].value)
    print("Feed salinity:                           ", m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].value *1e3, " g/kg")
    print("Brine salinity:                          ", m.fs.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3," g/kg")
    print("Vapor flow rate:                         ", m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value)
    print('Recovery:                                ', m.fs.recovery[0].value)
    print("Preheated feed temperature:              ", m.fs.evaporator.properties_feed[0].temperature.value)
    print("Evaporator temperature:                  ", m.fs.evaporator.properties_brine[0].temperature.value)
    print('Evaporator pressure:                     ', m.fs.evaporator.properties_vapor[0].pressure.value)
    print('Vapor temperature:                       ', m.fs.evaporator.properties_vapor[0].temperature.value)
    print("Compressed vapor temperature:            ", m.fs.compressor.control_volume.properties_out[0].temperature.value)
    print('Compressed vapor pressure:               ', m.fs.compressor.control_volume.properties_out[0].pressure.value)
    print("Condensed vapor temperature:             ", m.fs.condenser.control_volume.properties_out[0].temperature.value)
    print("Compressor work:                         ", m.fs.compressor.control_volume.work[0].value)
    print("Compressor pressure ratio:               ", m.fs.compressor.pressure_ratio.value)
    print("Evaporator area:                         ", m.fs.evaporator.area.value)
    print("Evaporator U:                            ", m.fs.evaporator.U.value)
    print('Evaporator LMTD:                         ', m.fs.evaporator.lmtd.value)
    print('Evaporator approach temperature in:      ', m.fs.evaporator.delta_temperature_in.value)
    print('Evaporator approach temperature out:     ', m.fs.evaporator.delta_temperature_out.value)
    print("Brine HX area:                           ", m.fs.hx_brine.area.value)
    print("Distillate HX area:                      ", m.fs.hx_distillate.area.value)
    print('Specific energy consumption:             ', value(m.fs.costing.specific_energy_consumption))
    print('Electricity cost:                        ', m.fs.costing.electricity_cost.value)
    print('HX unit cost:                            ', m.fs.costing.heat_exchanger.unit_cost.value)
    print('HX material factor:                      ', m.fs.costing.heat_exchanger.material_factor_cost.value)
    print('Evaporator unit cost:                    ', m.fs.costing.evaporator.unit_cost.value)
    print('Evaporator material factor:              ', m.fs.costing.evaporator.material_factor_cost.value)
    print('Total investment factor:                 ', m.fs.costing.factor_total_investment.value)
    print('LCOW:                                    ', m.fs.costing.LCOW.value)
    print('External Q:                              ', m.fs.Q_ext[0].value)


if __name__ == "__main__":
    m = main()
