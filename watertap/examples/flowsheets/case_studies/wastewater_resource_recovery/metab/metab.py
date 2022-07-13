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
    display_results(m.fs)

    add_costing(m)
    m.fs.costing.initialize()

    adjust_default_parameters(m)

    assert_degrees_of_freedom(m, 0)
    results = solve(m)
    assert_optimal_termination(results)
    display_costing(m.fs)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={"solute_list": ["cod", "hydrogen", "methane"]}
    )

    # unit models
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})
    m.fs.metab_hydrogen = MetabZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "hydrogen",
        },
    )
    m.fs.metab_methane = MetabZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "methane",
        },
    )
    m.fs.product_hydrogen = Product(default={"property_package": m.fs.prop})
    m.fs.product_methane = Product(default={"property_package": m.fs.prop})
    m.fs.product_H2O = Product(default={"property_package": m.fs.prop})

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


def display_results(fs):
    unit_list = ["feed", "metab_hydrogen", "metab_methane"]
    for u in unit_list:
        fs.component(u).report()


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "metab_global_costing.yaml",
    )
    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})
    # typing aid
    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.metab_hydrogen.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.metab_methane.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product_H2O.properties[0].flow_vol)

    # m.fs.costing.add_LCOW(m.fs.product_H2O.properties[0].flow_vol)

    # other levelized costs
    # TODO -resolve discrepancy between sum of cost components and total levelized costs
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


def display_costing(fs):
    fs.costing.total_capital_cost.display()
    fs.costing.total_operating_cost.display()
    fs.costing.LCOW.display()

    print("\nUnit Capital Costs\n")
    for u in fs.costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018)),
        )

    print("\nUtility Costs\n")
    for f in fs.costing.flow_types:
        print(
            f,
            " :   ",
            value(
                pyunits.convert(
                    fs.costing.aggregate_flow_costs[f],
                    to_units=pyunits.USD_2018 / pyunits.year,
                )
            ),
        )

    print("")
    total_capital_cost = value(
        pyunits.convert(fs.costing.total_capital_cost, to_units=pyunits.MUSD_2018)
    )
    print(f"Total Capital Costs: {total_capital_cost:.4f} M$")

    total_operating_cost = value(
        pyunits.convert(
            fs.costing.total_operating_cost, to_units=pyunits.MUSD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.4f} M$/year")

    electricity_intensity = value(
        pyunits.convert(
            fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m^3")
    LCOW = value(
        pyunits.convert(
            fs.costing.LCOW, to_units=fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m^3")
    LCOCR = value(
        pyunits.convert(
            fs.costing.LCOCR, to_units=fs.costing.base_currency / pyunits.kg
        )
    )
    print(f"Levelized Cost of COD Removal: {LCOCR:.4f} $/kg")
    LCOH = value(
        pyunits.convert(fs.costing.LCOH, to_units=fs.costing.base_currency / pyunits.kg)
    )
    print(f"Levelized Cost of Hydrogen: {LCOH:.4f} $/kg")
    LCOM = value(
        pyunits.convert(fs.costing.LCOM, to_units=fs.costing.base_currency / pyunits.kg)
    )
    print(f"Levelized Cost of Methane: {LCOM:.4f} $/kg")


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


if __name__ == "__main__":
    m, results = main()
