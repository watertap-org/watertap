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
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc

from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.generic_models.unit_models import Feed, Separator, Mixer, Product
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.generic_models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w


def main():
    # build
    m = build()
    set_operating_conditions(m)

    initialize_system(m)
    # m.fs.evaporator.connect_to_condenser(m.fs.condenser)
    print(degrees_of_freedom(m))
    print("Initialization")

    # assert False
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    recovery = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ].value / (
        m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
        + m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    print("First solve")
    print("Recovery: ", recovery)
    print(
        "Evaporator temperature: ",
        m.fs.evaporator.properties_brine[0].temperature.value,
    )
    print("Evaporator pressure: ", m.fs.evaporator.properties_brine[0].pressure.value)
    print("Evaporator area: ", m.fs.evaporator.area.value)
    print(
        "Compressor outlet temperature: ", m.fs.compressor.outlet.temperature[0].value
    )
    print("Compressor outlet pressure: ", m.fs.compressor.outlet.pressure[0].value)
    print("Compressor work: ", m.fs.compressor.control_volume.work[0].value)
    print("Condenser outlet temperature: ", m.fs.condenser.outlet.temperature[0].value)
    print("Condenser outlet pressure: ", m.fs.condenser.outlet.pressure[0].value)

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
    # m.fs.separator_feed = Separator(
    #     default={
    #         "property_package": m.fs.properties_feed,
    #         "outlet_list": [""]
    #     }
    # )
    m.fs.hx_distillate = HeatExchanger(
        default={
            "hot_side_name": "hot",
            "cold_side_name": "cold",
            "hot": {"property_package": m.fs.properties_vapor},
            "cold": {"property_package": m.fs.properties_feed},
            "delta_temperature_callback": delta_temperature_chen_callback,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
        }
    )
    # m.fs.mixer_feed = Mixer(
    #     default={
    #         "property_package": m.fs.properties_feed,
    #         "momentum_mixing_type": MomentumMixingType.equality,
    #         "inlet_list": []
    #     }
    # )

    m.fs.evaporator = Evaporator(
        default={
            "property_package_feed": m.fs.properties_feed,
            "property_package_vapor": m.fs.properties_vapor,
        }
    )

    m.fs.compressor = Compressor(default={"property_package": m.fs.properties_vapor})

    m.fs.condenser = Condenser(default={"property_package": m.fs.properties_vapor})

    m.fs.pump_brine = Pump(default={"property_package": m.fs.properties_feed})

    m.fs.pump_distillate = Pump(default={"property_package": m.fs.properties_vapor})

    m.fs.distillate = Product(default={"property_package": m.fs.properties_vapor})

    m.fs.brine = Product(default={"property_package": m.fs.properties_feed})

    # Connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_feed.inlet)
    m.fs.s02 = Arc(
        source=m.fs.pump_feed.outlet, destination=m.fs.hx_distillate.cold_inlet
    )
    m.fs.s03 = Arc(
        source=m.fs.hx_distillate.cold_outlet, destination=m.fs.evaporator.inlet_feed
    )
    m.fs.s04 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s05 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    m.fs.s06 = Arc(
        source=m.fs.evaporator.outlet_brine, destination=m.fs.pump_brine.inlet
    )
    m.fs.s07 = Arc(source=m.fs.pump_brine.outlet, destination=m.fs.brine.inlet)
    m.fs.s08 = Arc(source=m.fs.condenser.outlet, destination=m.fs.pump_distillate.inlet)
    m.fs.s09 = Arc(
        source=m.fs.pump_distillate.outlet, destination=m.fs.hx_distillate.hot_inlet
    )
    m.fs.s10 = Arc(
        source=m.fs.hx_distillate.hot_outlet, destination=m.fs.distillate.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.evaporator.connect_to_condenser(m.fs.condenser)

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
    # heat exchangers
    iscale.set_scaling_factor(m.fs.hx_distillate.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.hx_distillate.cold.heat, 1e-3)
    iscale.set_scaling_factor(
        m.fs.hx_distillate.overall_heat_transfer_coefficient, 1e-3
    )
    iscale.set_scaling_factor(m.fs.hx_distillate.area, 1)
    # evaporator
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)
    # iscale.set_scaling_factor(m.fs.evaporator.heat_transfer, 1e-6)

    # compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)

    # condenser
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-6)

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # Feed inlet
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)
    # evaporator feed inlet
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(10)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.05)

    # evaporator specifications
    m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(400)  # m^2
    # m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap','H2O'].fix(5)

    # compressor
    # m.fs.compressor.pressure_ratio.fix(2)
    m.fs.compressor.pressure_ratio = 2
    m.fs.compressor.control_volume.work.fix(5.8521e05)
    m.fs.compressor.efficiency.fix(0.8)

    # check degrees of freedom
    print(degrees_of_freedom(m))


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # initialize evaporator
    m.fs.evaporator.initialize_build(
        delta_temperature_in=30, delta_temperature_out=5
    )  # fixes and unfixes those values

    # initialize compressor
    propagate_state(m.fs.s01)
    m.fs.compressor.initialize_build()

    # initialize condenser
    propagate_state(m.fs.s02)
    m.fs.condenser.initialize_build(heat=-m.fs.evaporator.heat_transfer.value)
    # m.fs.condenser.outlet.temperature[0].fix(m.fs.evaporator.properties_brine[0].temperature.value+1)
    # m.fs.condenser.initialize_build()
    # m.fs.condenser.outlet.temperature[0].unfix()


def display_system(m):
    recovery = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ].value / (
        m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "TDS"].value
        + m.fs.evaporator.properties_feed[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    print("Recovery: ", recovery)


def display_stream_states(m):
    print("Feed")
    # print('Temperature: %.3f')
    print("Evaporator heat transfer: ", m.fs.evaporator.heat_transfer.value)
    print("Condenser heat transfer: ", m.fs.condenser.control_volume.heat[0].value)
    print("Feed inlet enth_flow: ", value(m.fs.evaporator.properties_feed[0].enth_flow))
    print(
        "Brine inlet enth_flow: ", value(m.fs.evaporator.properties_brine[0].enth_flow)
    )
    print(
        "Vapor inlet enth_flow: ",
        m.fs.evaporator.properties_vapor[0].enth_flow_phase["Vap"].value,
    )
    print(
        "Condenser inlet enth_flow: ",
        m.fs.condenser.control_volume.properties_in[0].enth_flow_phase["Vap"].value,
    )
    print(
        "Condenser outlet enth_flow: ",
        m.fs.condenser.control_volume.properties_out[0].enth_flow_phase["Liq"].value,
    )


if __name__ == "__main__":
    main()
