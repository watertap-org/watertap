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

from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Objective,
    Var,
    NonNegativeReals,
    TransformationFactory,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (
    propagate_state,
)
from idaes.models.unit_models import Heater, Separator, Mixer, Product, Feed
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
from idaes.core import UnitModelCostingBlock, FlowDirection
import idaes.core.util.scaling as iscale
from watertap.unit_models.MD.membrane_distillation_0D import MembraneDistillation0D
from watertap.unit_models.MD.MD_channel_base import (
    ConcentrationPolarizationType,
    TemperaturePolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.pump import (
    cost_pump,
)
from watertap.costing.unit_models.heater_chiller import (
    cost_heater_chiller,
)
from watertap.core.util.initialization import assert_degrees_of_freedom

__author__ = "Elmira Shamlou"


def main():
    solver = get_solver()
    m = build()
    set_operating_conditions(m)
    initialize_system(m, solver=solver)

    optimize_set_up(m)
    solve(m, solver=solver)

    print("\n***---optimization results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    return m


def build():

    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # properties
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.properties_hot_ch)
    m.fs.permeate = Product(property_package=m.fs.properties_cold_ch)
    m.fs.reject = Product(property_package=m.fs.properties_hot_ch)

    # overall process water recovery
    m.fs.overall_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        * m.fs.overall_recovery
        == m.fs.permeate.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )

    # unit models

    m.fs.MD = MembraneDistillation0D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.calculated,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.calculated,
            "flow_direction": FlowDirection.backward,
        },
    )

    m.fs.hx = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_hot_ch, "has_pressure_change": True},
        cold={"property_package": m.fs.properties_cold_ch, "has_pressure_change": True},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    # Set lower bound of approach temperatures
    m.fs.hx.delta_temperature_in.setlb(0)
    m.fs.hx.delta_temperature_out.setlb(0)
    m.fs.hx.area.setlb(10)

    m.fs.heater = Heater(
        property_package=m.fs.properties_hot_ch, has_pressure_change=True
    )
    m.fs.chiller = Heater(
        property_package=m.fs.properties_cold_ch, has_pressure_change=True
    )

    m.fs.eq_equal_flow = Constraint(
        expr=m.fs.chiller.control_volume.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ]
        == m.fs.heater.control_volume.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ]
    )

    m.fs.separator_permeate = Separator(
        property_package=m.fs.properties_cold_ch,
        outlet_list=["permeate", "cold_loop_stream"],
        split_basis=SplittingType.totalFlow,
    )
    m.fs.separator_concentrate = Separator(
        property_package=m.fs.properties_cold_ch,
        outlet_list=["reject", "recycle"],
        split_basis=SplittingType.totalFlow,
    )

    m.fs.recycle_ratio = Var(
        m.fs.config.time,
        initialize=10,
        bounds=(0, 100),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_recycle_ratio = Constraint(
        expr=sum(m.fs.feed.flow_mass_phase_comp[0, "Liq", j] for j in ["H2O", "TDS"])
        * m.fs.recycle_ratio[0]
        == sum(
            m.fs.separator_concentrate.recycle.flow_mass_phase_comp[0, "Liq", j]
            for j in ["H2O", "TDS"]
        )
    )

    m.fs.mixer = Mixer(
        property_package=m.fs.properties_hot_ch,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["feed", "recycle"],
    )
    m.fs.mixer.pressure_equality_constraints[0, 2].deactivate()

    m.fs.pump_feed = Pump(property_package=m.fs.properties_hot_ch)
    m.fs.pump_permeate = Pump(property_package=m.fs.properties_cold_ch)
    m.fs.pump_brine = Pump(property_package=m.fs.properties_hot_ch)

    # connections
    # brine (MD hot side) loop
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_feed.inlet)
    m.fs.s02 = Arc(source=m.fs.pump_feed.outlet, destination=m.fs.mixer.feed)
    m.fs.s03 = Arc(source=m.fs.mixer.outlet, destination=m.fs.hx.cold_inlet)
    m.fs.s04 = Arc(source=m.fs.hx.cold_outlet, destination=m.fs.pump_brine.inlet)
    m.fs.s05 = Arc(source=m.fs.pump_brine.outlet, destination=m.fs.heater.inlet)
    m.fs.s06 = Arc(source=m.fs.heater.outlet, destination=m.fs.MD.hot_ch_inlet)
    m.fs.s07 = Arc(
        source=m.fs.MD.hot_ch_outlet, destination=m.fs.separator_concentrate.inlet
    )
    m.fs.s08 = Arc(
        source=m.fs.separator_concentrate.reject, destination=m.fs.reject.inlet
    )
    m.fs.s09 = Arc(
        source=m.fs.separator_concentrate.recycle, destination=m.fs.mixer.recycle
    )

    # pemeate (MD cold side) loop
    m.fs.s10 = Arc(source=m.fs.chiller.outlet, destination=m.fs.MD.cold_ch_inlet)
    m.fs.s11 = Arc(source=m.fs.MD.cold_ch_outlet, destination=m.fs.hx.hot_inlet)
    m.fs.s12 = Arc(source=m.fs.hx.hot_outlet, destination=m.fs.separator_permeate.inlet)
    m.fs.s13 = Arc(
        source=m.fs.separator_permeate.permeate, destination=m.fs.permeate.inlet
    )
    m.fs.s14 = Arc(
        source=m.fs.separator_permeate.cold_loop_stream,
        destination=m.fs.pump_permeate.inlet,
    )
    m.fs.s15 = Arc(source=m.fs.pump_permeate.outlet, destination=m.fs.chiller.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)
    add_costs(m)

    # set default property values
    m.fs.properties_hot_ch.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_hot_ch.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_cold_ch.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_cold_ch.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    # set unit model values
    iscale.set_scaling_factor(m.fs.pump_permeate.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_brine.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.pump_feed.control_volume.work, 1e-3)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]

    iscale.set_scaling_factor(m.fs.hx.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx.cold.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(m.fs.hx.area, 1e-1)
    iscale.set_scaling_factor(m.fs.MD.hot_ch.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.MD.cold_ch.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.heater.control_volume.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.chiller.control_volume.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.MD.area, 1e-1)

    iscale.calculate_scaling_factors(m)

    return m


def add_costs(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.MD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.heater.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_heater_chiller,
        costing_method_arguments={"HC_type": "electric_heater"},
    )
    m.fs.chiller.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_heater_chiller,
        costing_method_arguments={"HC_type": "chiller"},
    )
    m.fs.pump_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_pump,
        costing_method_arguments={"pump_type": "low_pressure"},
    )
    m.fs.pump_permeate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_pump,
        costing_method_arguments={"pump_type": "low_pressure"},
    )
    m.fs.pump_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_pump,
        costing_method_arguments={"pump_type": "low_pressure"},
    )
    m.fs.hx.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.mixer.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.permeate.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.permeate.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.permeate.properties[0].flow_vol)


def set_operating_conditions(m):
    # overall recovery
    m.fs.overall_recovery.fix(0.5)
    # feed
    feed_flow_mass = 1
    feed_mass_frac_TDS = 0.035
    feed_pressure = 101325  # atmospheric
    feed_temperature = 273.15 + 25
    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.feed.properties[0].pressure.fix(feed_pressure)  # Pa
    m.fs.feed.properties[0].temperature.fix(feed_temperature)  # K

    # concentrate separator
    # either the overall recovery or the split fraction will be fixed, but not both
    recycle_fraction = 0.9
    recycle_flow_mass = 10 * feed_flow_mass
    # MD
    membrane_pressure_drop = -5e5
    membrane_area = 10 * recycle_flow_mass
    m.fs.MD.area.fix(membrane_area)  # m^2
    m.fs.MD.permeability_coef.fix(1e-10)  # kg/m-s-Pa
    m.fs.MD.membrane_thickness.fix(1e-4)  # m
    m.fs.MD.membrane_thermal_conductivity.fix(0.2)  # W/K-m
    length = 3
    m.fs.MD.length.fix(length)  # m

    m.fs.MD.hot_ch.channel_height.fix(0.0019)  # m
    m.fs.MD.hot_ch.spacer_porosity.fix(0.77)
    m.fs.MD.cold_ch.channel_height.fix(0.0019)  # m
    m.fs.MD.cold_ch.spacer_porosity.fix(0.77)

    # heat exchanger
    heat_exchanger_area = 4 * membrane_area
    heat_exchanger_pressure_drop = -1e5
    m.fs.hx.overall_heat_transfer_coefficient.fix(2e3)  # W/K-m^2
    m.fs.hx.area.fix(heat_exchanger_area)  # m^2
    m.fs.hx.cold.deltaP[0].fix(heat_exchanger_pressure_drop)  # Pa
    m.fs.hx.hot.deltaP[0].fix(heat_exchanger_pressure_drop)  # Pa

    # feed pump
    m.fs.pump_feed.efficiency_pump.fix(0.8)
    m.fs.pump_feed.control_volume.deltaP[0].fix(7e5)  # Pa

    # brine pump
    LEP = 7e5
    heater_pressure_drop = -1e5
    m.fs.pump_brine.efficiency_pump.fix(0.8)
    m.fs.pump_brine.control_volume.properties_out[0].pressure.fix(
        feed_pressure + LEP + heater_pressure_drop + heat_exchanger_pressure_drop
    )

    # permeate pump
    chiller_pressure_drop = -1e5
    m.fs.pump_permeate.efficiency_pump.fix(0.8)
    m.fs.pump_permeate.control_volume.properties_out[0].pressure.fix(
        feed_pressure + LEP + chiller_pressure_drop + heat_exchanger_pressure_drop
    )

    # heater
    heater_top_temperature = 273.15 + 90
    mixed_brine_mass_frac_TDS = 0.07
    mixed_brine_mass_frac_H2O = 1 - mixed_brine_mass_frac_TDS
    mixed_brine_flow_mass = feed_flow_mass + recycle_flow_mass
    m.fs.heater.control_volume.properties_out[0].temperature.fix(
        heater_top_temperature
    )  # K
    m.fs.heater.control_volume.deltaP[0].fix(heater_pressure_drop)
    m.fs.heater.control_volume.properties_out[0].pressure.value = (
        m.fs.pump_brine.control_volume.properties_out[0].pressure.value
        - heater_pressure_drop
    )
    m.fs.heater.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "H2O"
    ].value = (mixed_brine_flow_mass * mixed_brine_mass_frac_H2O)
    m.fs.heater.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "TDS"
    ].value = (mixed_brine_flow_mass * mixed_brine_mass_frac_TDS)

    # chiller
    chiller_outlet_temperature = 273.15 + 10
    # The flow rate in MD cold and hot channels is assumed to be equal.
    cold_loop_stream_flow_mass = recycle_flow_mass + feed_flow_mass
    m.fs.chiller.control_volume.properties_out[0].temperature.fix(
        chiller_outlet_temperature
    )  # K
    m.fs.chiller.control_volume.deltaP[0].fix(chiller_pressure_drop)
    m.fs.chiller.control_volume.properties_out[0].pressure.value = (
        m.fs.pump_permeate.control_volume.properties_out[0].pressure.value
        - chiller_pressure_drop
    )

    m.fs.chiller.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "TDS"
    ].value = 0
    m.fs.chiller.control_volume.properties_out[0].flow_mass_phase_comp[
        "Liq", "H2O"
    ] = m.fs.heater.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )


def solve(blk, solver=None, tee=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if not check_optimal_termination(results):
        results = solver.solve(blk, tee=tee)
    return results


def initialize_system(m, solver=None, verbose=True):
    if solver is None:
        solver = get_solver()

    # initialize feed block
    m.fs.feed.initialize()

    # with fixed chiller and heater outlet temperature, and initial guess for the flowrates, MD is initialized
    propagate_state(m.fs.s06)
    propagate_state(m.fs.s10)
    m.fs.MD.initialize()

    # hot loop
    propagate_state(m.fs.s07)
    m.fs.separator_concentrate.initialize()

    propagate_state(m.fs.s08)
    m.fs.reject.initialize()

    propagate_state(m.fs.s01)
    m.fs.pump_feed.initialize()

    propagate_state(m.fs.s02)
    propagate_state(m.fs.s09)
    m.fs.mixer.initialize()
    m.fs.mixer.pressure_equality_constraints[0, 2].deactivate()

    propagate_state(m.fs.s03)
    propagate_state(m.fs.s11)
    m.fs.hx.initialize()

    propagate_state(m.fs.s04)
    m.fs.pump_brine.initialize()

    propagate_state(m.fs.s05)
    m.fs.heater.initialize()

    # cold loop
    propagate_state(m.fs.s12)
    m.fs.separator_permeate.initialize()

    propagate_state(m.fs.s13)
    m.fs.permeate.initialize()

    propagate_state(m.fs.s14)
    m.fs.pump_permeate.initialize()

    propagate_state(m.fs.s15)
    m.fs.chiller.initialize()

    propagate_state(m.fs.s06)
    propagate_state(m.fs.s10)
    m.fs.MD.initialize()

    print(f"DOF: {degrees_of_freedom(m)}")

    m.fs.costing.initialize()


def optimize_set_up(m):
    # add objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    m.fs.MD.area.unfix()
    m.fs.MD.area.setlb(1)
    m.fs.MD.area.setub(150)

    m.fs.MD.length.unfix()
    m.fs.MD.length.setlb(0)
    m.fs.MD.length.setub(150)

    m.fs.hx.area.unfix()
    m.fs.hx.area.setlb(1)
    m.fs.hx.area.setub(750)

    m.fs.heater.control_volume.properties_out[0].temperature.unfix()
    m.fs.heater.control_volume.properties_out[0].temperature.setlb(273.15 + 10)
    m.fs.heater.control_volume.properties_out[0].temperature.setub(273.15 + 90)

    m.fs.chiller.control_volume.properties_out[0].temperature.unfix()
    m.fs.chiller.control_volume.properties_out[0].temperature.setlb(273.15 + 10)
    m.fs.chiller.control_volume.properties_out[0].temperature.setub(273.15 + 25)

    m.fs.pump_brine.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_brine.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.pump_brine.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.pump_brine.deltaP.setlb(0)

    # additional constraints
    LEP = 10e5
    m.fs.eq_liquid_entry_pressure = Constraint(
        expr=m.fs.MD.hot_ch_inlet.pressure[0] <= LEP
    )

    assert_degrees_of_freedom(m, 6)


def optimize(m, solver=None):
    # --solve---
    return solve(m, solver=solver)


def display_system(m):
    print("---system metrics---")
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    feed_mass_frac_TDS = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].value / feed_flow_mass
    )
    print("Feed: %.2f kg/s, %.0f ppm" % (feed_flow_mass, feed_mass_frac_TDS * 1e6))

    permeate_flow_mass = sum(
        m.fs.permeate.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    permeate_mass_frac_TDS = (
        m.fs.permeate.flow_mass_phase_comp[0, "Liq", "TDS"].value / permeate_flow_mass
    )
    print(
        "Permeate: %.3f kg/s, %.0f ppm"
        % (permeate_flow_mass, permeate_mass_frac_TDS * 1e6)
    )

    reject_flow_mass = sum(
        m.fs.reject.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    reject_mass_frac_TDS = (
        m.fs.reject.flow_mass_phase_comp[0, "Liq", "TDS"].value / reject_flow_mass
    )
    print(
        "Reject: %.3f kg/s, %.0f ppm" % (reject_flow_mass, reject_mass_frac_TDS * 1e6)
    )

    print("Overall recovery: %.1f%%" % (value(m.fs.overall_recovery) * 100))

    print("Recycle ratio: %.1f" % (value(m.fs.recycle_ratio[0])))

    print("Thermal efficiency: %.1f%%" % (value(m.fs.MD.thermal_efficiency[0]) * 100))

    print("effectiveness: %.1f%%" % (value(m.fs.MD.effectiveness[0]) * 100))

    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))


def display_design(m):
    print("---decision variables---")
    print("Membrane area %.1f m2" % (m.fs.MD.area.value))
    print("Heat exchanger area %.1f m2" % (m.fs.hx.area.value))
    print(
        "Heater\noutlet temperature: %.1f Celcius\npower %.2f kW"
        % (
            m.fs.heater.outlet.temperature[0].value - 273.15,
            m.fs.heater.heat_duty[0].value / 1e3,
        )
    )
    print(
        "Chiller\noutlet temperature: %.1f Celcius\npower %.2f kW"
        % (
            m.fs.chiller.outlet.temperature[0].value - 273.15,
            m.fs.chiller.heat_duty[0].value / 1e3,
        )
    )
    print(
        "brine pump\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.pump_brine.outlet.pressure[0].value / 1e5,
            m.fs.pump_brine.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "permeate pump\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.pump_permeate.outlet.pressure[0].value / 1e5,
            m.fs.pump_permeate.work_mechanical[0].value / 1e3,
        )
    )


def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "TDS"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        Temperature_Celsius = b.temperature[0].value - 273.15
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.1f bar, %.1f °C"
            % (flow_mass, mass_frac_ppm, pressure_bar, Temperature_Celsius)
        )

    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Feed pump out    ", m.fs.pump_feed.outlet)
    print_state("Mixer out    ", m.fs.mixer.outlet)
    print_state("Heat exchanger hot channel in", m.fs.hx.hot_side_inlet)
    print_state("Heat exchanger hot channel out", m.fs.hx.hot_side_outlet)
    print_state("Heat exchanger cold channel in", m.fs.hx.cold_side_inlet)
    print_state("Heat exchanger cold channel out", m.fs.hx.cold_side_outlet)
    print_state("brine pump out    ", m.fs.pump_brine.outlet)
    print_state("Heater out    ", m.fs.heater.outlet)
    print_state("MD hot channel in", m.fs.MD.hot_ch_inlet)
    print_state("MD hot channel out", m.fs.MD.hot_ch_outlet)
    print_state("MD cold channel in", m.fs.MD.cold_ch_inlet)
    print_state("MD cold channel out", m.fs.MD.cold_ch_outlet)
    print_state("Concentrate separator recycle ", m.fs.separator_concentrate.recycle)
    print_state("Concentrate separator reject ", m.fs.separator_concentrate.reject)
    print_state("permeate pump out    ", m.fs.pump_permeate.outlet)
    print_state("Concentrate permeate permeate ", m.fs.separator_permeate.permeate)
    print_state(
        "Concentrate permeate cold loop stream ",
        m.fs.separator_permeate.cold_loop_stream,
    )
    print_state("chiller out    ", m.fs.chiller.outlet)


if __name__ == "__main__":
    m = main()
