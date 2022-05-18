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
import os
from pyomo.environ import (
    ConcreteModel,
    Set,
    Expression,
    value,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.generic_models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    PumpElectricityZO,
    NanofiltrationZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting


def main():
    m = build()
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    # display_results(m)

    add_costing(m)
    m.fs.costing.initialize()
    #
    # adjust_default_parameters(m)
    results = solve(m)
    assert_optimal_termination(results)
    display_costing(m)
    return m


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(default={"solute_list": ["dye", "tds"]})

    # unit model
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})

    m.fs.P1 = PumpElectricityZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "default",
        }
    )

    m.fs.nanofiltration = NanofiltrationZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "rHGO_dye_rejection",
        }
    )

    m.fs.permeate1 = Product(default={"property_package": m.fs.prop})
    m.fs.retentate1 = Product(default={"property_package": m.fs.prop})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.nanofiltration.inlet)
    m.fs.s03 = Arc(source=m.fs.nanofiltration.treated, destination=m.fs.permeate1.inlet)
    m.fs.s04 = Arc(
        source=m.fs.nanofiltration.byproduct, destination=m.fs.retentate1.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # feed
    flow_vol = 120 / 3600 * pyunits.m**3 / pyunits.s
    conc_mass_dye = 2.5 * pyunits.kg / pyunits.m**3
    conc_mass_tds = 50.0 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "dye"].fix(conc_mass_dye)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    solve(m.fs.feed)

    # nanofiltration
    m.fs.nanofiltration.load_parameters_from_database(use_default_removal=True)

    # pump
    m.fs.P1.load_parameters_from_database(use_default_removal=True)
    m.fs.P1.applied_pressure.fix(m.fs.nanofiltration.applied_pressure.get_values()[0])
    m.fs.P1.lift_height.unfix()


def initialize_system(m):
    seq = SequentialDecomposition()
    seq.options.tear_set = []
    seq.options.iterLim = 1
    seq.run(m, lambda u: u.initialize())


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dye_desalination_global_costing.yaml",
    )

    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})

    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.nanofiltration.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.P1.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.costing.cost_process()

    # # electricity per cost of influent

    m.fs.costing.annual_energy_cost = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * m.fs.costing.electricity_cost
            * pyunits.convert(
                m.fs.P1.electricity[0],
                to_units=pyunits.kWh / m.fs.costing.base_period,
            )
        ),
        doc="Energy cost associated with pumping",
    )

    # amount of dye in retentate - permeate = dye removed [kg/time]
    m.fs.costing.annual_dye_removal = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.retentate1.flow_mass_comp[0, "dye"]
                - m.fs.permeate1.flow_mass_comp[0, "dye"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        ),
        doc="Annual dye removal",
    )

    # total cost per year - as defined in metab.py
    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            + m.fs.costing.annual_energy_cost
        ),
        doc="Annual cost of treatment",
    )

    # levelized cost of dye removal
    m.fs.costing.LCODS = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_dye_removal),
        doc="Levelized Cost of Dye Separation",
    )

    # # annual cost of dye disposal
    m.fs.costing.annual_dye_disposal_cost = Expression(
        expr=(m.fs.costing.annual_dye_removal * m.fs.costing.dye_disposal_cost),
        doc="Cost of dye disposal",
    )

    m.fs.costing.LCODS = Expression(
        expr=m.fs.costing.utilization_factor
        * m.fs.costing.total_annualized_cost
        / pyunits.convert(
            m.fs.feed.flow_vol[0], to_units=pyunits.m**3 / m.fs.costing.base_period
        ),
        doc="Levelized cost of dye separation, on a feed flow volume basis",
    )


def display_results(m):
    unit_list = ["feed", "permeate1", "retentate1", "nanofiltration"]
    for u in unit_list:
        m.fs.component(u).report()


def display_costing(m):
    m.fs.costing.annual_dye_disposal_cost.display()
    m.fs.costing.annual_dye_removal.display()
    m.fs.costing.annual_energy_cost.display()

    print("\nUnit Capital Costs\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018)),
        )

    print("")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2018)
    )
    print(f"Total Capital Costs: {total_capital_cost:.4f} M$")

    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.4f} M$/year")

    levelized_cost_dye_separation = value(
        pyunits.convert(m.fs.costing.LCODS, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(
        f"\nLevelized Cost of Dye Separation: {levelized_cost_dye_separation:.2f} $/m3"
    )


if __name__ == "__main__":
    model = main()
