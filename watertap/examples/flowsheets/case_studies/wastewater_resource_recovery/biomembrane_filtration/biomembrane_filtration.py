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
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.core import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    MABRZO,
    DMBRZO,
    PumpElectricityZO,
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
    # display_reports(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    results = solve(m)
    assert_optimal_termination(results)

    display_metrics_results(m)
    display_additional_results(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=["bod", "tss", "ammonium_as_nitrogen", "nitrate", "nitrogen"]
    )

    # unit models
    m.fs.feed = FeedZO(property_package=m.fs.prop)
    m.fs.mabr = MABRZO(property_package=m.fs.prop, database=m.db)
    m.fs.dmbr = DMBRZO(property_package=m.fs.prop, database=m.db)
    m.fs.P1 = PumpElectricityZO(property_package=m.fs.prop, database=m.db)

    m.fs.product_H2O = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.mabr.inlet)
    m.fs.s03 = Arc(source=m.fs.mabr.treated, destination=m.fs.dmbr.inlet)
    m.fs.s04 = Arc(source=m.fs.dmbr.treated, destination=m.fs.product_H2O.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 417 * pyunits.m**3 / pyunits.hr
    conc_mass_bod = 103 * pyunits.mg / pyunits.liter
    conc_mass_tss = 100 * pyunits.mg / pyunits.liter
    conc_mass_ammonium_as_nitrogen = 26 * pyunits.mg / pyunits.liter
    conc_mass_nitrate = 1.3 * pyunits.mg / pyunits.liter

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "bod"].fix(conc_mass_bod)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"].fix(
        conc_mass_ammonium_as_nitrogen
    )
    m.fs.feed.conc_mass_comp[0, "nitrate"].fix(conc_mass_nitrate)
    m.fs.feed.conc_mass_comp[0, "nitrogen"].fix(1e-5)
    solve(m.fs.feed)

    # pump
    m.fs.P1.load_parameters_from_database(use_default_removal=True)
    m.fs.P1.lift_height.fix(1)

    # mabr
    m.fs.mabr.load_parameters_from_database(use_default_removal=True)

    # dmbr
    m.fs.dmbr.load_parameters_from_database(use_default_removal=True)


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


def display_reports(m):
    unit_list = ["feed", "P1", "mabr", "dmbr"]
    for u in unit_list:
        m.fs.component(u).report()


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "biomembrane_filtration_global_costing.yaml",
    )
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    # typing aid
    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.mabr.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.dmbr.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.P1.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product_H2O.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product_H2O.properties[0].flow_vol)

    # other levelized costs
    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.annual_water_production = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.product_H2O.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
    )

    m.fs.costing.LCOT = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_water_inlet),
        doc="Levelized Cost of Treatment",
    )


def initialize_costing(m):
    m.fs.costing.initialize()


def display_metrics_results(m):
    print("----------Levelized costs----------")
    LCOT = value(
        pyunits.convert(
            m.fs.costing.LCOT, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Treatment: {LCOT:.4f} $/m3 of feed")
    LCOW = value(
        pyunits.convert(
            m.fs.costing.LCOW, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m3 of product")

    print("----------Capital costs----------")
    DCC_normalized = value(
        pyunits.convert(
            (
                m.fs.mabr.costing.direct_capital_cost
                + m.fs.dmbr.costing.direct_capital_cost
                + m.fs.P1.costing.direct_capital_cost
            )
            / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.4f} $/(m3/day)")
    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.4f} $/(m3/day)")

    print("----------Operating costs----------")
    FMC_normalized = value(
        pyunits.convert(
            m.fs.costing.maintenance_cost / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.a,
        )
    )
    print(f"Normalized maintenance costs: {FMC_normalized:.4f} 1/year")
    BRC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_fixed_operating_cost
            / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.a,
        )
    )
    print(f"Normalized bead replacement cost: {BRC_normalized:.4f} 1/year")
    EC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.4f} $/m3 of feed")

    print("----------Performance metrics----------")
    volumetric_recovery = value(
        m.fs.product_H2O.properties[0].flow_vol / m.fs.feed.properties[0].flow_vol
    )
    print(f"Water recovery: {volumetric_recovery:.4f} m3 of product/m3 of feed")
    BODR_normalized = value(
        pyunits.convert(
            1
            - m.fs.product_H2O.properties[0].flow_mass_comp["bod"]
            / m.fs.feed.properties[0].flow_mass_comp["bod"],
            to_units=pyunits.dimensionless,
        )
    )
    print(f"BOD removal: {BODR_normalized:.4f} dimensionless")
    TSSR_normalized = value(
        pyunits.convert(
            1
            - m.fs.product_H2O.properties[0].flow_mass_comp["tss"]
            / m.fs.feed.properties[0].flow_mass_comp["tss"],
            to_units=pyunits.dimensionless,
        )
    )
    print(f"TSS removal: {TSSR_normalized:.4f} dimensionless")
    NR_normalized = value(
        pyunits.convert(
            1
            - (
                m.fs.product_H2O.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
                + m.fs.product_H2O.properties[0].flow_mass_comp["nitrate"]
            )
            / (
                m.fs.feed.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
                + m.fs.feed.properties[0].flow_mass_comp["nitrate"]
            ),
            to_units=pyunits.dimensionless,
        )
    )
    print(f"Total N removal: {NR_normalized:.4f} dimensionless")

    print("----------Energy intensity----------")
    SEC = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Specific electricity consumption: {SEC:.4f} kWh/m3 of feed")


def display_additional_results(m):
    print("----------Outlets----------")
    product_H2O_flow = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.hr,
        )
    )
    print(f"H2O outlet flow: {product_H2O_flow:.4f} m3/h")
    product_H2O_BOD = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["bod"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet BOD conc: {product_H2O_BOD:.4f} g/L")
    product_H2O_TSS = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["tss"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet TSS conc: {product_H2O_TSS:.4f} g/L")
    product_H2O_total_nitrogen = value(
        pyunits.convert(
            (
                m.fs.product_H2O.properties[0].conc_mass_comp["ammonium_as_nitrogen"]
                + m.fs.product_H2O.properties[0].conc_mass_comp["nitrate"]
            ),
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet total nitrogen conc: {product_H2O_total_nitrogen:.4f} g/L")

    print("----------Capital costs----------")
    total_capital_costs = value(m.fs.costing.total_capital_cost) / 1e6
    print(f"Total capital costs: {total_capital_costs:.4f} $M")
    mabr_capital_costs = value(m.fs.mabr.costing.capital_cost) / 1e6
    print(f"MABR capital costs: {mabr_capital_costs:.4f} $M")
    dmbr_capital_costs = value(m.fs.dmbr.costing.capital_cost) / 1e6
    print(f"DMBR capital costs: {dmbr_capital_costs:.4f} $M")
    pump_capital_costs = value(m.fs.P1.costing.capital_cost) / 1e6
    print(f"Pump capital costs: {pump_capital_costs:.4f} $M")

    print("----------Operating costs----------")
    total_operating_costs = value(m.fs.costing.total_operating_cost) / 1e6
    print(f"Total operating costs: {total_operating_costs:.4f} $M/year")
    fixed_operating_costs = value(m.fs.costing.total_fixed_operating_cost) / 1e6
    print(f"Fixed operating costs: {fixed_operating_costs:.4f} $M/year")
    electricity_operating_costs = (
        value(m.fs.costing.aggregate_flow_costs["electricity"]) / 1e3
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.4f} $k/year")


if __name__ == "__main__":
    m, results = main()
