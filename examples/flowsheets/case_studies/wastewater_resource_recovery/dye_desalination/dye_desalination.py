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
    Block,
    Expression,
    value,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.core import UnitModelCostingBlock

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

    display_results(m)

    add_costing(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m)
    assert_optimal_termination(results)

    display_costing(m)
    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(solute_list=["dye", "tds"])

    # unit model
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # define block to integrate with dye_desalination_withRO
    dye_sep = m.fs.dye_separation = Block()

    dye_sep.P1 = PumpElectricityZO(
        property_package=m.fs.prop, database=m.db, process_subtype="default"
    )

    dye_sep.nanofiltration = NanofiltrationZO(
        property_package=m.fs.prop, database=m.db, process_subtype="rHGO_dye_rejection"
    )

    m.fs.permeate = Product(property_package=m.fs.prop)
    m.fs.dye_retentate = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=dye_sep.P1.inlet)
    m.fs.s02 = Arc(source=dye_sep.P1.outlet, destination=dye_sep.nanofiltration.inlet)
    m.fs.s03 = Arc(
        source=dye_sep.nanofiltration.treated, destination=m.fs.permeate.inlet
    )
    m.fs.s04 = Arc(
        source=dye_sep.nanofiltration.byproduct, destination=m.fs.dye_retentate.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    dye_sep = m.fs.dye_separation
    # feed
    flow_vol = 120 / 3600 * pyunits.m**3 / pyunits.s
    conc_mass_dye = 2.5 * pyunits.kg / pyunits.m**3
    conc_mass_tds = 50.0 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "dye"].fix(conc_mass_dye)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    solve(m.fs.feed)

    # nanofiltration
    dye_sep.nanofiltration.load_parameters_from_database(use_default_removal=True)

    # pump
    dye_sep.P1.load_parameters_from_database(use_default_removal=True)
    dye_sep.P1.applied_pressure.fix(
        dye_sep.nanofiltration.applied_pressure.get_values()[0]
    )
    dye_sep.P1.lift_height.unfix()

    return


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
    # initialize block
    dye_sep = m.fs.dye_separation

    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dye_desalination_global_costing.yaml",
    )

    # zero order costing
    m.fs.zo_costing = ZeroOrderCosting(case_study_definition=source_file)

    costing_kwargs = {"flowsheet_costing_block": m.fs.zo_costing}

    # create costing blocks
    dye_sep.nanofiltration.costing = UnitModelCostingBlock(**costing_kwargs)
    dye_sep.P1.costing = UnitModelCostingBlock(**costing_kwargs)

    # aggregate unit level costs
    m.fs.zo_costing.cost_process()

    # create system level cost metrics
    m.fs.brine_disposal_cost = Expression(
        expr=(
            m.fs.zo_costing.utilization_factor
            * (
                m.fs.zo_costing.waste_disposal_cost
                * pyunits.convert(
                    m.fs.permeate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            )
        ),
        doc="Cost of disposing of saline brine/ NF permeate",
    )

    m.fs.dye_recovery_revenue = Expression(
        expr=(
            m.fs.zo_costing.utilization_factor
            * m.fs.zo_costing.dye_mass_cost
            * pyunits.convert(
                m.fs.dye_retentate.flow_mass_comp[0, "dye"],
                to_units=pyunits.kg / m.fs.zo_costing.base_period,
            )
        ),
        doc="Savings from dye-retentate recovered back to the plant",
    )

    # combine results for system level costs - to be the same syntax as dye_desalination_withRO
    @m.fs.Expression(doc="Total capital cost")
    def total_capital_cost(b):
        return pyunits.convert(
            m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2020
        )

    @m.fs.Expression(doc="Total operating cost")
    def total_operating_cost(b):
        return pyunits.convert(
            b.zo_costing.total_fixed_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        ) + pyunits.convert(
            b.zo_costing.total_variable_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )

    @m.fs.Expression(doc="Total cost of dye recovered and brine disposed")
    def total_externalities(b):
        return pyunits.convert(
            b.dye_recovery_revenue - b.brine_disposal_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )

    @m.fs.Expression(
        doc="Levelized cost of treatment with respect to volumetric feed flowrate"
    )
    def LCOT(b):
        return (
            b.total_capital_cost * b.zo_costing.capital_recovery_factor
            + b.total_operating_cost
            - b.total_externalities
        ) / (
            pyunits.convert(
                b.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.zo_costing.utilization_factor
        )

    assert_units_consistent(m)
    m.fs.zo_costing.initialize()
    return


def display_results(m):
    unit_list = ["P1", "nanofiltration"]
    for u in unit_list:
        m.fs.dye_separation.component(u).report()


def display_costing(m):
    print("\nUnit Capital Costs")
    for u in m.fs.zo_costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020)),
            "$",
        )

    print("\nSystem Costs")
    total_capital_cost = value(
        pyunits.convert(m.fs.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.4f} M$")

    total_operating_cost = value(
        pyunits.convert(
            m.fs.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.4f} M$/year")

    total_externalities = value(
        pyunits.convert(
            m.fs.total_externalities, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Externalities: {total_externalities:.4f} M$/year")

    levelized_cost_treatment = value(
        pyunits.convert(m.fs.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"Levelized Cost of Treatment (LCOT): {levelized_cost_treatment:.2f} $/m3 feed"
    )


if __name__ == "__main__":
    m, results = main()
