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
    units as pyunits,
    assert_optimal_termination,
    Expression,
    value,
    TransformationFactory,
)

from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

# from idaes.generic_models.unit_models import Product
from idaes.models.unit_models import Product
import idaes.core.util.scaling as iscale

# from idaes.generic_models.costing import UnitModelCostingBlock
from idaes.core import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    ClothMediaFiltrationZO,
    AnaerobicDigestionReactiveZO,
    MembraneEvaporatorZO,
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
    m.fs.costing.initialize()

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)
    results = solve(m)
    assert_optimal_termination(results)
    display_costing(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=[
            "tss",
            "cod",
            "tkn",
            "acetic_acid",
            "ammonium_as_nitrogen",
        ]
    )

    # unit models
    m.fs.cmf = ClothMediaFiltrationZO(
        property_package=m.fs.prop,
        database=m.db,
    )
    m.fs.ad = AnaerobicDigestionReactiveZO(
        property_package=m.fs.prop,
        database=m.db,
    )
    m.fs.me = MembraneEvaporatorZO(
        property_package=m.fs.prop,
        database=m.db,
    )

    # feed and product streams
    m.fs.feed = FeedZO(property_package=m.fs.prop)
    m.fs.filtered_water = Product(property_package=m.fs.prop)
    m.fs.ad_byproduct = Product(property_package=m.fs.prop)
    m.fs.me_treated = Product(property_package=m.fs.prop)
    m.fs.me_byproduct = Product(property_package=m.fs.prop)

    # connections
    m.fs.s1 = Arc(source=m.fs.feed.outlet, destination=m.fs.cmf.inlet)
    m.fs.s2 = Arc(source=m.fs.cmf.treated, destination=m.fs.filtered_water.inlet)
    m.fs.s3 = Arc(source=m.fs.cmf.byproduct, destination=m.fs.ad.inlet)
    m.fs.s4 = Arc(source=m.fs.ad.byproduct, destination=m.fs.ad_byproduct.inlet)
    m.fs.s5 = Arc(source=m.fs.ad.treated, destination=m.fs.me.inlet)
    m.fs.s6 = Arc(source=m.fs.me.treated, destination=m.fs.me_treated.inlet)
    m.fs.s7 = Arc(source=m.fs.me.byproduct, destination=m.fs.me_byproduct.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # feed conditions
    flow_vol = 25000000.0 * pyunits.gal / pyunits.day
    conc_mass_tss = 500.0 * pyunits.mg / pyunits.L
    conc_mass_cod = 700.0 * pyunits.mg / pyunits.L
    conc_mass_tkn = 65 * pyunits.mg / pyunits.L
    conc_mass_aa = 1e-8 * pyunits.mg / pyunits.L
    conc_mass_nh4 = 1e-8 * pyunits.mg / pyunits.L

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "cod"].fix(conc_mass_cod)
    m.fs.feed.conc_mass_comp[0, "tkn"].fix(conc_mass_tkn)
    m.fs.feed.conc_mass_comp[0, "acetic_acid"].fix(conc_mass_aa)
    m.fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"].fix(conc_mass_nh4)

    solve(m.fs.feed)

    m.fs.cmf.load_parameters_from_database(use_default_removal=True)
    m.fs.ad.load_parameters_from_database(use_default_removal=True)
    m.fs.me.load_parameters_from_database(use_default_removal=True)


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


def display_results(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY RESULTS ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    unit_list = ["cmf", "ad", "me"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n---------- Feed volumetric flowrate ----------")
    feed_vol_flow = value(
        pyunits.convert(m.fs.feed.flow_vol[0], to_units=pyunits.m**3 / pyunits.hr)
    )
    print(f"{feed_vol_flow:.0f} m^3/hr")

    print("\n---------- Feed water concentrations ----------")
    tss_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "tss"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"TSS: {tss_feed_conc:.1f} mg/L")

    cod_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "cod"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"COD: {cod_feed_conc:.1f} mg/L")

    tkn_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "tkn"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"TKN: {tkn_feed_conc:.1f} mg/L")

    aa_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "acetic_acid"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Acetic acid: {aa_feed_conc:.1f} mg/L")

    nh4_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"TSS: {nh4_feed_conc:.1f} mg/L")

    print("\n---------- Anaerobic digester results ----------")
    biogas_production = value(
        pyunits.convert(
            m.fs.ad.biogas_production[0], to_units=pyunits.m**3 / pyunits.hr
        )
    )
    print(f"Biogas production: {biogas_production:.2f} m^3/hr")

    byproduct_tss_rate = value(
        pyunits.convert(
            m.fs.ad_byproduct.flow_mass_comp[0, "tss"], to_units=pyunits.kg / pyunits.hr
        )
    )
    print(f"TSS byproduct production: {byproduct_tss_rate:.2f} kg/hr")

    byproduct_cod_rate = value(
        pyunits.convert(
            m.fs.ad_byproduct.flow_mass_comp[0, "cod"], to_units=pyunits.kg / pyunits.hr
        )
    )
    print(f"COD byproduct production: {byproduct_cod_rate:.2f} kg/hr")

    byproduct_tkn_rate = value(
        pyunits.convert(
            m.fs.ad_byproduct.flow_mass_comp[0, "tkn"], to_units=pyunits.kg / pyunits.hr
        )
    )
    print(f"TKN byproduct production: {byproduct_tkn_rate:.2f} kg/hr")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "amo_1690_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.cmf.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.ad.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.me.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.biogas_recovery_volume = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.ad.biogas_production[0],
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        ),
        doc="Volume of biogas recovered per year",
    )

    m.fs.costing.value_biogas_recovery = Expression(
        expr=m.fs.costing.biogas_recovery_volume * m.fs.costing.biogas_cost,
        doc="Value of biogas recovered per year",
    )

    m.fs.costing.fertilizer_recovery_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.me_byproduct.flow_mass_comp[0, "ammonium_as_nitrogen"]
            + m.fs.me_byproduct.flow_mass_comp[0, "acetic_acid"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Volume of fertilizer (NH4 + acetic acid) recovered per year",
    )

    m.fs.costing.value_fertilizer_recovery = Expression(
        expr=m.fs.costing.fertilizer_recovery_mass * m.fs.costing.fertilizer_cost,
        doc="Value of fertilizer recovered per year",
    )

    m.fs.costing.LCOT = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate",
    )

    m.fs.costing.LCOT_with_revenue = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            - m.fs.costing.value_biogas_recovery
            - m.fs.costing.value_fertilizer_recovery
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate including revenue from biogas and fertilizer",
    )

    m.fs.costing.LC_biogas = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / m.fs.costing.biogas_recovery_volume,
        doc="Levelized Cost of Biogas",
    )

    m.fs.costing.LC_biogas_with_revenue = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            - m.fs.costing.value_fertilizer_recovery
        )
        / m.fs.costing.biogas_recovery_volume,
        doc="Levelized Cost of Biogas including revenue from fertilizer",
    )

    m.fs.costing.LC_fertilizer = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / m.fs.costing.fertilizer_recovery_mass,
        doc="Levelized Cost of Fertilizer",
    )

    m.fs.costing.LC_fertilizer_with_revenue = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            - m.fs.costing.value_biogas_recovery
        )
        / m.fs.costing.fertilizer_recovery_mass,
        doc="Levelized Cost of Fertilize including revenue from biogas",
    )


def display_costing(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY COSTING ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    print("\n---------- Levelized costs ----------")
    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Treatment: {LCOT:.3f} $/m^3 feed water")

    LCOT_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOT_with_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(
        f"Levelized Cost of Treatment including revenue: {LCOT_with_revenue:.3f} $/m^3 feed water"
    )

    LC_biogas = value(
        pyunits.convert(
            m.fs.costing.LC_biogas, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Biogas: {LC_biogas:.3f} $/m^3")

    LC_biogas_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LC_biogas_with_revenue,
            to_units=pyunits.USD_2020 / pyunits.m**3,
        )
    )
    print(
        f"Levelized Cost of Biogas including fertilizer revenue: {LC_biogas_with_revenue:.3f} $/m^3"
    )

    LC_fertilizer = value(
        pyunits.convert(
            m.fs.costing.LC_fertilizer, to_units=pyunits.USD_2020 / pyunits.kg
        )
    )
    print(f"Levelized Cost of Fertilizer: {LC_fertilizer:.3f} $/kg")

    LC_fertilizer_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LC_fertilizer_with_revenue,
            to_units=pyunits.USD_2020 / pyunits.kg,
        )
    )
    print(
        f"Levelized Cost of Fertilizer including biogas revenue: {LC_fertilizer_with_revenue:.3f} $/kg"
    )

    print("\n------------- Capital costs -------------")
    DCC_normalized = value(
        pyunits.convert(
            (
                m.fs.cmf.costing.direct_capital_cost
                + m.fs.ad.costing.direct_capital_cost
                + m.fs.me.costing.direct_capital_cost
            )
            / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kUSD_2020 / (pyunits.m**3 / pyunits.hr),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.2f} k$/(m^3 feed/hr)")

    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kUSD_2020 / (pyunits.m**3 / pyunits.hr),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.2f} k$/(m^3/hr)")

    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.2f} M$")

    print(
        f"Cloth media filtration capital cost: {value(pyunits.convert(m.fs.cmf.costing.capital_cost, to_units=pyunits.MUSD_2020)):.2f} M$"
    )

    print(
        f"Anaerobic digester capital cost: {value(pyunits.convert(m.fs.ad.costing.capital_cost, to_units=pyunits.MUSD_2020)):.2f} M$"
    )

    print(
        f"Membrane evaporator capital cost: {value(pyunits.convert(m.fs.me.costing.capital_cost, to_units=pyunits.MUSD_2020)):.2f} M$"
    )

    print("\n------------- Operating costs -------------")
    FMC_normalized = value(
        pyunits.convert(
            m.fs.costing.maintenance_cost / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.year,
        )
    )
    print(f"Normalized maintenance costs: {FMC_normalized:.2f} 1/year")

    OFOC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_fixed_operating_cost
            / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.year,
        )
    )
    print(f"Normalized other fixed operating cost: {OFOC_normalized:.2f} 1/year")

    EC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            * m.fs.costing.utilization_factor
            / m.fs.costing.annual_water_inlet,
            to_units=pyunits.USD_2020 / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.3f} $/m^3 of feed")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(
        f"Electricity Intensity with respect to influent flowrate: {electricity_intensity:.3f} kWh/m^3"
    )

    total_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.kUSD_2020 / pyunits.year
        )
    )
    print(f"Total operating costs: {total_operating_costs:.1f} k$/year")

    fixed_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_fixed_operating_cost,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Fixed operating costs: {fixed_operating_costs:.1f} k$/year")

    variable_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_variable_operating_cost,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Variable operating costs: {variable_operating_costs:.1f} k$/year")

    electricity_operating_costs = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            * m.fs.costing.utilization_factor,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.1f} k$/year")

    print("\n------------- Product recovery and revenue -------------")
    biogas_recovery = value(
        pyunits.convert(
            m.fs.costing.biogas_recovery_volume, to_units=pyunits.m**3 / pyunits.year
        )
    )
    print(f"Biogas recovery: {biogas_recovery:.0f} m^3/year")

    biogas_recovery_normalized = value(
        pyunits.convert(
            m.fs.costing.biogas_recovery_volume / m.fs.costing.annual_water_inlet,
            to_units=pyunits.dimensionless,
        )
    )
    print(
        f"Normalized biogas recovery: {biogas_recovery_normalized:.3f} m^3 biogas/m^3 feed water"
    )

    biogas_revenue = value(
        pyunits.convert(
            m.fs.costing.value_biogas_recovery,
            to_units=pyunits.MUSD_2020 / pyunits.year,
        )
    )
    print(f"Biogas revenue: {biogas_revenue:.2f} M$/year")

    fertilizer_recovery = value(
        pyunits.convert(
            m.fs.costing.fertilizer_recovery_mass, to_units=pyunits.kg / pyunits.year
        )
    )
    print(f"Fertilizer recovery: {fertilizer_recovery:.0f} kg/year")

    fertilizer_recovery_normalized = value(
        pyunits.convert(
            m.fs.costing.fertilizer_recovery_mass / m.fs.costing.annual_water_inlet,
            to_units=pyunits.kg / pyunits.m**3,
        )
    )
    print(
        f"Normalized fertilizer recovery: {fertilizer_recovery_normalized:.3f} kg fertilizer/m^3 feed water"
    )

    fertilizer_revenue = value(
        pyunits.convert(
            m.fs.costing.value_fertilizer_recovery,
            to_units=pyunits.MUSD_2020 / pyunits.year,
        )
    )
    print(f"Fertilizer revenue: {fertilizer_revenue:.2f} M$/year")

    total_revenue = value(
        pyunits.convert(
            m.fs.costing.value_biogas_recovery + m.fs.costing.value_fertilizer_recovery,
            to_units=pyunits.MUSD_2020 / pyunits.year,
        )
    )
    print(f"Total revenue: {total_revenue:.2f} M$/year")

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
