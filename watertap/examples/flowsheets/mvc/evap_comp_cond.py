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


def main():
    # build
    m = build()
    set_operating_conditions(m)
    #set_operating_conditions_validation(m)
    initialize_system(m)
    scale_costs(m)
    print('DOF after initialization: ', degrees_of_freedom(m))
    m.fs.objective = Objective(expr=1)
    solver = get_solver()
    results = solve(m)

    print('First solve - simulation')
    print(results.solver.termination_condition)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)
    display_results(m)
    assert False

    print('Second solve')
    results = solve(m)
    print(results.solver.termination_condition)
    display_results(m)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)

    print('Third solve')
    results = solve(m)
    print(results.solver.termination_condition)
    display_results(m)
    if results.solver.termination_condition == "infeasible":
        debug_infeasible(m.fs, solver)


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties_feed})

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

    m.fs.distillate = Product(default={"property_package": m.fs.properties_feed})

    m.fs.brine = Product(default={"property_package": m.fs.properties_feed})

    # Connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.evaporator.inlet_feed)
    m.fs.s02 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s03 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    m.fs.s04 = Arc(
        source=m.fs.evaporator.outlet_brine, destination=m.fs.brine.inlet)
    m.fs.s05 = Arc(source=m.fs.condenser.outlet, destination=m.fs.tb_distillate.inlet)
    m.fs.s06 = Arc(source=m.fs.tb_distillate.outlet, destination=m.fs.distillate.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.evaporator.connect_to_condenser(m.fs.condenser)

    # Add costing
    add_costing(m)

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

def add_costing(m):
    m.fs.costing = WaterTAPCosting()
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

def set_operating_conditions(m):
    m_f = 40  # kg/s
    w_f = 50 * 1e-3  # kg/kg
    rr = 0.25
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix((1 - w_f) * m_f)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(w_f * m_f)
    m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(rr*m_f)

    # Feed inlet
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40) # 10
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(2) # 0.5
    m.fs.feed.properties[0].temperature.fix(273.15 + 50)
    m.fs.feed.properties[0].pressure.fix(101325)

    # evaporator specifications
    #m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(1000)  # m^2
    #m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix(20) # 5

    # compressor
    m.fs.compressor.pressure_ratio.fix(1.5)
    #m.fs.compressor.control_volume.properties_out[0].temperature = 400
    m.fs.compressor.efficiency.fix(0.8)

    # Fix 0 TDS
    m.fs.tb_distillate.properties_out[0].flow_mass_phase_comp['Liq','TDS'].fix(1e-5)

    # rescale based on mass flow rate
    sf = 10 ** -(math.log10(m.fs.feed.properties[0].flow_mass_phase_comp['Liq','H2O'].value))
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", sf, index=("Liq", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", sf, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", sf, index=("Liq", "H2O")
    )

    sf = 10 ** -(math.log10(m.fs.feed.properties[0].flow_mass_phase_comp['Liq','TDS'].value))
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", sf, index=("Liq", "TDS")
    )
    iscale.calculate_scaling_factors(m)

    # check degrees of freedom
    print("DOF after setting operating conditions: ", degrees_of_freedom(m))


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # Touch feed mass fraction property
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS']

    # initialize evaporator
    propagate_state(m.fs.s01)
    m.fs.evaporator.initialize(
        delta_temperature_in=70, delta_temperature_out=9
    )  # fixes and unfixes those values

    # initialize compressor
    propagate_state(m.fs.s02)
    m.fs.compressor.initialize()

    # initialize condenser
    propagate_state(m.fs.s03)
    m.fs.condenser.initialize(heat=-m.fs.evaporator.heat_transfer.value)

    # propagate translator block
    propagate_state(m.fs.s05) # to translator block
    #m.fs.pump_distillate.initialize(optarg=optarg)

    m.fs.costing.initialize()

    print('Initialization done')

def calculate_cost_sf(cost):
    sf = 10**-(math.log10(cost.value))
    iscale.set_scaling_factor(cost, sf)

def scale_costs(m):
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


def optimize(m):
    m.fs.objective = Objective(
        expr=-m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )
    print("Set objective")
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    recovery = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ].value / (
        m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
        + m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    print("Recovery after optimization: ", recovery)


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

    m.fs.evaporator.area.unfix()
    m.fs.evaporator.inlet_feed.temperature[0].unfix()
    m.fs.evaporator.outlet_brine.temperature[0].unfix()
    m.fs.compressor.control_volume.work[0].unfix()
    m.fs.compressor.pressure_ratio.unfix()
    m.fs.hx_distillate.area.unfix()
    m.fs.hx_brine.area.unfix()


def optimize(m):
    print("\nSet objective to minimize cost")
    set_up_optimization(m)
    #results = solver.solve(m, tee=False)
    #assert_optimal_termination(results)


def display_results(m):
    print("Feed flow rate:                          ", m.fs.feed.properties[0].flow_mass_phase_comp['Liq','H2O'].value+
          m.fs.feed.properties[0].flow_mass_phase_comp['Liq','TDS'].value)
    print("Feed mass fraction:                      ", m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].value)
    print("Vapor flow rate:                         ", m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value)
    print("Preheated feed temperature:              ", m.fs.evaporator.properties_feed[0].temperature.value)
    print("Evaporator temperature:                  ", m.fs.evaporator.properties_brine[0].temperature.value)
    print("Evaporator pressure:                     ", m.fs.evaporator.properties_vapor[0].pressure.value)
    print("Compressed vapor temperature:            ", m.fs.compressor.control_volume.properties_out[0].temperature.value)
    print("Compressed vapor pressure:               ", m.fs.compressor.control_volume.properties_out[0].pressure.value)
    print("Compressed vapor specific enthalpy:      ", m.fs.compressor.control_volume.properties_out[0].enth_mass_phase['Vap'].value)
    print("Condensed vapor temperature:             ", m.fs.condenser.control_volume.properties_out[0].temperature.value)
    print("Condensed vapor pressure:                ", m.fs.condenser.control_volume.properties_out[0].pressure.value)
    print("Compressor work:                         ", m.fs.compressor.control_volume.work[0].value)
    print("Specific work:                           ", m.fs.compressor.control_volume.work[0].value/m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap','H2O'].value)
    print("Compressor pressure ratio:               ", m.fs.compressor.pressure_ratio.value)
    print("Evaporator area:                         ", m.fs.evaporator.area.value)
    print('Evaporator LMTD:                         ', m.fs.evaporator.lmtd.value)
    print('Specific energy consumption:             ', value(m.fs.costing.specific_energy_consumption))
    print('LCOW:                                    ', m.fs.costing.LCOW.value)

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


def display_water_states(state_blk):
    print("Liquid mass flow ", state_blk.flow_mass_phase_comp["Liq", "H2O"].value)
    print("Vapor mass flow  ", state_blk.flow_mass_phase_comp["Vap", "H2O"].value)
    print("temperature      ", state_blk.temperature.value)
    print("pressure         ", state_blk.pressure.value)


if __name__ == "__main__":
    main()
