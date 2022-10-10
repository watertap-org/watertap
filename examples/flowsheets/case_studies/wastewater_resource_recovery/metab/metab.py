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
    MetabZO,
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

    add_costing(m)
    assert_degrees_of_freedom(m, 0)
    m.fs.costing.initialize()

    # adjust_default_parameters(m)

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
    m.fs.prop = prop_ZO.WaterParameterBlock(solute_list=["cod", "hydrogen", "methane"])

    # unit models
    m.fs.feed = FeedZO(property_package=m.fs.prop)
    m.fs.metab_hydrogen = MetabZO(
        property_package=m.fs.prop, database=m.db, process_subtype="hydrogen"
    )
    m.fs.metab_methane = MetabZO(
        property_package=m.fs.prop, database=m.db, process_subtype="methane"
    )
    m.fs.product_hydrogen = Product(property_package=m.fs.prop)
    m.fs.product_methane = Product(property_package=m.fs.prop)
    m.fs.product_H2O = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.metab_hydrogen.inlet)
    m.fs.s02 = Arc(
        source=m.fs.metab_hydrogen.treated, destination=m.fs.metab_methane.inlet
    )
    m.fs.s03 = Arc(
        source=m.fs.metab_hydrogen.byproduct, destination=m.fs.product_hydrogen.inlet
    )
    m.fs.s04 = Arc(
        source=m.fs.metab_methane.byproduct, destination=m.fs.product_methane.inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.metab_methane.treated, destination=m.fs.product_H2O.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 3.286e-4 * pyunits.m**3 / pyunits.s
    conc_mass_cod = 6.76 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "cod"].fix(conc_mass_cod)
    m.fs.feed.properties[0].flow_mass_comp["hydrogen"].fix(1e-8)
    m.fs.feed.properties[0].flow_mass_comp["methane"].fix(1e-8)
    solve(m.fs.feed)

    # metab_hydrogen
    m.fs.metab_hydrogen.load_parameters_from_database(use_default_removal=True)

    # metab_methane
    m.fs.metab_methane.load_parameters_from_database(use_default_removal=True)


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


def display_reports(fs):
    unit_list = ["feed", "metab_hydrogen", "metab_methane"]
    for u in unit_list:
        fs.component(u).report()


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "metab_global_costing.yaml",
    )
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    # typing aid
    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.metab_hydrogen.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.metab_methane.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product_H2O.properties[0].flow_vol)

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
    m.fs.costing.annual_cod_removal = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.feed.outlet.flow_mass_comp[0, "cod"]
                - m.fs.product_H2O.inlet.flow_mass_comp[0, "cod"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    m.fs.costing.annual_hydrogen_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.metab_hydrogen.byproduct.flow_mass_comp[0, "hydrogen"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    m.fs.costing.annual_methane_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.metab_methane.byproduct.flow_mass_comp[0, "methane"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
    )

    m.fs.costing.LCOW = Expression(
        expr=(
            m.fs.costing.total_annualized_cost / m.fs.costing.annual_water_production
        ),
        doc="Levelized Cost of Water",
    )

    m.fs.costing.LCOT = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_water_inlet),
        doc="Levelized Cost of Treatment",
    )

    m.fs.costing.LCOCR = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_cod_removal),
        doc="Levelized Cost of COD Removal",
    )

    m.fs.costing.LCOH = Expression(
        expr=(
            (
                m.fs.metab_hydrogen.costing.capital_cost
                * m.fs.costing.capital_recovery_factor
                + m.fs.metab_hydrogen.costing.capital_cost
                * m.fs.costing.maintenance_costs_percent_FCI
                + m.fs.metab_hydrogen.costing.fixed_operating_cost
                + (
                    pyunits.convert(
                        m.fs.metab_hydrogen.heat[0] * m.fs.costing.heat_cost,
                        to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
                    )
                    + pyunits.convert(
                        m.fs.metab_hydrogen.electricity[0]
                        * m.fs.costing.electricity_cost,
                        to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
                    )
                )
                * m.fs.costing.utilization_factor
            )
            / m.fs.costing.annual_hydrogen_production
        ),
        doc="Levelized Cost of Hydrogen",
    )

    m.fs.costing.LCOM = Expression(
        expr=(
            (
                m.fs.metab_methane.costing.capital_cost
                * m.fs.costing.capital_recovery_factor
                + m.fs.metab_methane.costing.capital_cost
                * m.fs.costing.maintenance_costs_percent_FCI
                + m.fs.metab_methane.costing.fixed_operating_cost
                + (
                    pyunits.convert(
                        m.fs.metab_methane.heat[0] * m.fs.costing.heat_cost,
                        to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
                    )
                    + pyunits.convert(
                        m.fs.metab_methane.electricity[0]
                        * m.fs.costing.electricity_cost,
                        to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
                    )
                )
                * m.fs.costing.utilization_factor
            )
            / m.fs.costing.annual_methane_production
        ),
        doc="Levelized Cost of Methane",
    )

    m.fs.costing.LC_comp = Set(
        initialize=[
            "bead",
            "reactor",
            "mixer",
            "membrane",
            "vacuum",
            "heat",
            "electricity_mixer",
            "electricity_vacuum",
            "hydrogen_product",
            "methane_product",
        ]
    )
    m.fs.costing.LCOH_comp = Expression(m.fs.costing.LC_comp)

    m.fs.costing.LCOH_comp["bead"] = (
        m.fs.metab_hydrogen.costing.DCC_bead
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
        + m.fs.metab_hydrogen.costing.fixed_operating_cost
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["reactor"] = (
        m.fs.metab_hydrogen.costing.DCC_reactor
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["mixer"] = (
        m.fs.metab_hydrogen.costing.DCC_mixer
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["vacuum"] = (
        m.fs.metab_hydrogen.costing.DCC_vacuum
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["membrane"] = (
        m.fs.metab_hydrogen.costing.DCC_membrane
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["electricity_vacuum"] = (
        pyunits.convert(
            m.fs.metab_hydrogen.energy_electric_vacuum_flow_vol_byproduct
            * m.fs.metab_hydrogen.properties_byproduct[0].flow_mass_comp["hydrogen"]
            * m.fs.costing.electricity_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["electricity_mixer"] = (
        pyunits.convert(
            m.fs.metab_hydrogen.energy_electric_mixer_vol
            * m.fs.metab_hydrogen.volume
            * m.fs.costing.electricity_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOH_comp["heat"] = (
        pyunits.convert(
            m.fs.metab_hydrogen.heat[0] * m.fs.costing.heat_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_hydrogen_production

    m.fs.costing.LCOM_comp = Expression(m.fs.costing.LC_comp)

    m.fs.costing.LCOM_comp["bead"] = (
        m.fs.metab_methane.costing.DCC_bead
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
        + m.fs.metab_methane.costing.fixed_operating_cost
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["reactor"] = (
        m.fs.metab_methane.costing.DCC_reactor
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["mixer"] = (
        m.fs.metab_methane.costing.DCC_mixer
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["vacuum"] = (
        m.fs.metab_methane.costing.DCC_vacuum
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["membrane"] = (
        m.fs.metab_methane.costing.DCC_membrane
        * m.fs.costing.TIC
        * (
            m.fs.costing.capital_recovery_factor
            + m.fs.costing.maintenance_costs_percent_FCI
        )
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["electricity_vacuum"] = (
        pyunits.convert(
            m.fs.metab_methane.energy_electric_vacuum_flow_vol_byproduct
            * m.fs.metab_methane.properties_byproduct[0].flow_mass_comp["methane"]
            * m.fs.costing.electricity_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["electricity_mixer"] = (
        pyunits.convert(
            m.fs.metab_methane.energy_electric_mixer_vol
            * m.fs.metab_methane.volume
            * m.fs.costing.electricity_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_methane_production

    m.fs.costing.LCOM_comp["heat"] = (
        pyunits.convert(
            m.fs.metab_methane.heat[0] * m.fs.costing.heat_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
        * m.fs.costing.utilization_factor
    ) / m.fs.costing.annual_methane_production

    def rule_LCOW_comp(b, c):
        if c in ["hydrogen_product", "methane_product"]:
            return (
                m.fs.costing.aggregate_flow_costs[c]
                * m.fs.costing.utilization_factor
                / m.fs.costing.annual_water_production
            )
        else:
            return (
                m.fs.costing.LCOH_comp[c] * m.fs.costing.annual_hydrogen_production
                + m.fs.costing.LCOM_comp[c] * m.fs.costing.annual_methane_production
            ) / m.fs.costing.annual_water_production

    m.fs.costing.LCOW_comp = Expression(m.fs.costing.LC_comp, rule=rule_LCOW_comp)

    def rule_LCOCR_comp(b, c):
        if c in ["hydrogen_product", "methane_product"]:
            return (
                m.fs.costing.aggregate_flow_costs[c]
                * m.fs.costing.utilization_factor
                / m.fs.costing.annual_cod_removal
            )
        else:
            return (
                m.fs.costing.LCOH_comp[c] * m.fs.costing.annual_hydrogen_production
                + m.fs.costing.LCOM_comp[c] * m.fs.costing.annual_methane_production
            ) / m.fs.costing.annual_cod_removal

    m.fs.costing.LCOCR_comp = Expression(m.fs.costing.LC_comp, rule=rule_LCOCR_comp)


def adjust_default_parameters(m):
    m.fs.metab_hydrogen.hydraulic_retention_time.fix(6)  # default - 12 hours, 0.5x
    m.fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"].set_value(
        0.05
    )  # default - 0.005, 10x
    m.fs.costing.metab.bead_bulk_density["hydrogen"].fix(7.17)  # default 23.9, 0.3x
    m.fs.costing.metab.bead_replacement_factor["hydrogen"].fix(1)  # default 3.376, 0.3x
    m.fs.metab_hydrogen.energy_electric_mixer_vol.fix(0.049875)  # default 0.049875
    m.fs.metab_hydrogen.energy_electric_vacuum_flow_vol_byproduct.fix(
        9.190
    )  # default 9190, 0.001x
    m.fs.metab_hydrogen.energy_thermal_flow_vol_inlet.fix(7875)  # default 78750, 0.1x
    m.fs.costing.metab.bead_cost["hydrogen"].fix(14.40)  # default 1440, 0.01x
    m.fs.costing.metab.reactor_cost["hydrogen"].fix(78.9)  # default 789, 0.1x
    m.fs.costing.metab.vacuum_cost["hydrogen"].fix(5930)  # default 59300, 0.1x
    m.fs.costing.metab.mixer_cost["hydrogen"].fix(27.40)  # default 2740, 0.01x
    m.fs.costing.metab.membrane_cost["hydrogen"].fix(498)  # default 498

    m.fs.metab_methane.hydraulic_retention_time.fix(15)  # default 150, 0.1x
    m.fs.metab_methane.generation_ratio["cod_to_methane", "methane"].set_value(
        0.101
    )  # default 0.101, no change
    m.fs.costing.metab.bead_bulk_density["methane"].fix(7.17)  # default 23.9, 0.3x
    m.fs.costing.metab.bead_replacement_factor["methane"].fix(1)  # default 3.376, 0.3x
    m.fs.metab_methane.energy_electric_mixer_vol.fix(0.049875)  # default 0.049875
    m.fs.metab_methane.energy_electric_vacuum_flow_vol_byproduct.fix(
        1.53
    )  # default 15.3, 0.1x
    m.fs.metab_methane.energy_thermal_flow_vol_inlet.fix(0)  # default 0
    m.fs.costing.metab.bead_cost["methane"].fix(14.40)  # default 1440, 0.01x
    m.fs.costing.metab.reactor_cost["methane"].fix(78.9)  # default 789, 0.1x
    m.fs.costing.metab.vacuum_cost["methane"].fix(136.0)  # default 1360, 0.1x
    m.fs.costing.metab.mixer_cost["methane"].fix(27.40)  # default 2740, 0.01x
    m.fs.costing.metab.membrane_cost["methane"].fix(498)  # default 498


def display_metrics_results(m):
    print("----------Levelized costs----------")
    LCOT = value(
        pyunits.convert(
            m.fs.costing.LCOT, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Treatment: {LCOT:.2f} $/m3 of feed")
    LCOW = value(
        pyunits.convert(
            m.fs.costing.LCOW, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Water: {LCOW:.2f} $/m3 of product")
    LCOH = value(
        pyunits.convert(
            m.fs.costing.LCOH, to_units=m.fs.costing.base_currency / pyunits.kg
        )
    )
    print(f"Levelized Cost of Hydrogen: {LCOH:.2f} $/kg")
    LCOM = value(
        pyunits.convert(
            m.fs.costing.LCOM, to_units=m.fs.costing.base_currency / pyunits.kg
        )
    )
    print(f"Levelized Cost of Methane: {LCOM:.2f} $/kg")
    LCOCR = value(
        pyunits.convert(
            m.fs.costing.LCOCR, to_units=m.fs.costing.base_currency / pyunits.kg
        )
    )
    print(f"Levelized Cost of COD Removal: {LCOCR:.2f} $/kg")

    print("----------Capital costs----------")
    DCC_normalized = value(
        pyunits.convert(
            (
                m.fs.metab_hydrogen.costing.direct_capital_cost
                + m.fs.metab_methane.costing.direct_capital_cost
            )
            / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.2f} $/(m3/day)")
    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.2f} $/(m3/day)")

    print("----------Operating costs----------")
    FMC_normalized = value(
        pyunits.convert(
            m.fs.costing.maintenance_cost / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.a,
        )
    )
    print(f"Normalized maintenance costs: {FMC_normalized:.3f} 1/year")
    BRC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_fixed_operating_cost
            / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.a,
        )
    )
    print(f"Normalized bead replacement cost: {BRC_normalized:.3f} 1/year")
    EC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.2f} $/m3 of feed")
    HC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["heat"] / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized heating cost: {HC_normalized:.2f} $/m3 of feed")

    print("----------Revenue----------")
    H2R_normalized = value(
        pyunits.convert(
            -m.fs.costing.aggregate_flow_costs["hydrogen_product"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized hydrogen revenue: {H2R_normalized:.2f} $/m3 of feed")
    CH4R_normalized = value(
        pyunits.convert(
            -m.fs.costing.aggregate_flow_costs["methane_product"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized methane revenue: {CH4R_normalized:.2f} $/m3 of feed")

    print("----------Performance metrics----------")
    volumetric_recovery = value(
        m.fs.product_H2O.properties[0].flow_vol / m.fs.feed.properties[0].flow_vol
    )
    print(f"Water recovery: {volumetric_recovery:.3f} m3 of product/m3 of feed")
    CODR_normalized = value(
        pyunits.convert(
            1
            - m.fs.product_H2O.properties[0].flow_mass_comp["cod"]
            / m.fs.feed.properties[0].flow_mass_comp["cod"],
            to_units=pyunits.dimensionless,
        )
    )
    print(f"COD removal: {CODR_normalized:.4f} dimensionless")
    H2P_normalized = value(
        pyunits.convert(
            m.fs.product_hydrogen.properties[0].flow_mass_comp["hydrogen"]
            / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kg / pyunits.m**3,
        )
    )
    print(f"Normalized hydrogen production: {H2P_normalized:.4f} kg/m3")
    CH4P_normalized = value(
        pyunits.convert(
            m.fs.product_methane.properties[0].flow_mass_comp["methane"]
            / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kg / pyunits.m**3,
        )
    )
    print(f"Normalized methane production: {CH4P_normalized:.4f} kg/m3")

    print("----------Energy intensity----------")
    SEC = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Specific electricity consumption: {SEC:.3f} kWh/m3 of feed")
    STC = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_heat / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.MJ / pyunits.m**3,
        )
    )
    print(f"Specific thermal consumption: {STC:.3f} MJ/m3 of feed")


def display_additional_results(m):
    print("----------Outlets----------")
    product_H2O_flow = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.hr,
        )
    )
    print(f"H2O outlet flow: {product_H2O_flow:.2f} m3/h")
    product_H2O_COD = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["cod"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet COD conc: {product_H2O_COD:.2f} g/L")
    product_H2_flow = value(
        pyunits.convert(
            m.fs.product_hydrogen.properties[0].flow_mass_comp["hydrogen"],
            to_units=pyunits.kg / pyunits.hr,
        )
    )
    print(f"H2 outlet flow: {product_H2_flow:.4f} kg/h")
    product_CH4_flow = value(
        pyunits.convert(
            m.fs.product_methane.properties[0].flow_mass_comp["methane"],
            to_units=pyunits.kg / pyunits.hr,
        )
    )
    print(f"CH4 outlet flow: {product_CH4_flow:.3f} kg/h")

    print("----------Capital costs----------")
    total_capital_costs = value(m.fs.costing.total_capital_cost) / 1e6
    print(f"Total capital costs: {total_capital_costs:.1f} $M")
    hydrogen_capital_costs = value(m.fs.metab_hydrogen.costing.capital_cost) / 1e6
    print(f"Hydrogen capital costs: {hydrogen_capital_costs:.2f} $M")
    methane_capital_costs = value(m.fs.metab_methane.costing.capital_cost) / 1e6
    print(f"Methane capital costs: {methane_capital_costs:.1f} $M")

    print("----------Operating costs----------")
    total_operating_costs = value(m.fs.costing.total_operating_cost) / 1e6
    print(f"Total operating costs: {total_operating_costs:.1f} $M/year")
    fixed_operating_costs = value(m.fs.costing.total_fixed_operating_cost) / 1e6
    print(f"Fixed operating costs: {fixed_operating_costs:.1f} $M/year")
    electricity_operating_costs = (
        value(m.fs.costing.aggregate_flow_costs["electricity"]) / 1e3
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.1f} $k/year")
    heating_operating_costs = value(m.fs.costing.aggregate_flow_costs["heat"]) / 1e3
    print(f"Heating operating costs: {heating_operating_costs:.1f} $k/year")

    print("----------Revenue----------")
    total_revenue = value(
        -(
            m.fs.costing.aggregate_flow_costs["hydrogen_product"]
            + m.fs.costing.aggregate_flow_costs["methane_product"]
        )
    )
    print(f"Total revenue: {total_revenue:.1f} $/year")
    hydrogen_revenue = value(-(m.fs.costing.aggregate_flow_costs["hydrogen_product"]))
    print(f"Hydrogen revenue: {hydrogen_revenue:.1f} $/year")
    methane_revenue = value(-(m.fs.costing.aggregate_flow_costs["methane_product"]))
    print(f"Methane revenue: {methane_revenue:.1f} $/year")


if __name__ == "__main__":
    m, results = main()
