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
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.NDMA_prop_pack as props
from watertap.unit_models.uv_aop import Ultraviolet0D
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting


def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)
    # display_results(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m)
    display_costing(m)

    return m, results

def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NDMAParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.unit = Ultraviolet0D(default={"property_package": m.fs.properties})
    m.fs.product = Product(default={"property_package": m.fs.properties})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.unit.inlet)
    m.fs.s02 = Arc(source=m.fs.unit.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e4, index=("Liq", "NDMA")
    )
    # set unit model values (no value)

    # unused scaling factors needed by IDAES base costing module
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def set_operating_conditions(m):
    # ---specifications---
    feed_flow_mass = 2053 * pyunits.kg / pyunits.s
    feed_mass_frac_NDMA = 74e-9
    feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
    feed_pressure = 101325 * pyunits.Pa
    feed_temperature = (273.15 + 25) * pyunits.K
    uv_intensity = 1 * pyunits.mW / pyunits.cm ** 2
    exporure_time = 500 * pyunits.s
    inactivation_rate = 2.3 * pyunits.cm ** 2 / pyunits.J
    reaction_rate_constant = 0 * pyunits.min ** -1
    EEO = 0.25 * pyunits.kWh / pyunits.m ** 3
    lamp_efficiency = 0.8
    UVT = 0.9

    # feed
    # state variables
    m.fs.unit.inlet.pressure[0].fix(feed_pressure) # feed pressure [Pa]
    m.fs.unit.inlet.temperature[0].fix(feed_temperature) # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
        feed_flow_mass * feed_mass_frac_NDMA
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    # uv aop unit
    m.fs.unit.uv_intensity.fix(uv_intensity)
    m.fs.unit.exposure_time.fix(exporure_time)
    m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
    m.fs.unit.reaction_rate_constant["Liq", "NDMA"].fix(reaction_rate_constant)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure)
    m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
    m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
    m.fs.unit.UVT.fix(UVT)

def initialize_system(m):
    # ---initialize UV-AOP---
    m.fs.unit.initialize(outlvl=idaeslog.DEBUG)

    # ---initialize feed block---
    m.fs.feed.initialize(outlvl=idaeslog.DEBUG)

def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results

def display_results(m):
    unit_list = ["feed", "unit"]
    for u in unit_list:
        m.fs.component(u).report()

def add_costing(m):
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.unit.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
        },
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)


def initialize_costing(m):
    m.fs.costing.initialize()

def display_costing(m):
    m.fs.costing.total_capital_cost.display()
    m.fs.costing.total_operating_cost.display()
    m.fs.costing.LCOW.display()

    print("\nUnit Capital Costs\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020)),
        )

    # print("\nUtility Costs\n")
    # for f in m.fs.costing.flow_types:
    #     print(
    #         f,
    #         " :   ",
    #         value(
    #             pyunits.convert(
    #                 m.fs.costing.aggregate_flow_costs[f],
    #                 to_units=pyunits.USD_2020 / pyunits.year,
    #             )
    #         ),
    #     )

    print("")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.2f} M$")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.2f} M$/year")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m^3")

if __name__ == "__main__":
    m, results = main()