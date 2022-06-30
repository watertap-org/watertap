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
    AnaerobicMBRMECZO,
    CofermentationZO,
    ConstructedWetlandsZO,
    GasSpargedMembraneZO,
    IonExchangeZO,
    SedimentationZO,
    VFARecoveryZO,
)
from idaes.models.unit_models import Mixer, MomentumMixingType, MixingType
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)
    display_results(m.fs)

    # add_costing(m)
    # m.fs.costing.initialize()
    #
    # adjust_default_parameters(m)
    #
    # assert_degrees_of_freedom(m, 0)
    # results = solve(m)
    # assert_optimal_termination(results)
    # display_costing(m.fs)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={
            "solute_list": [
                "cod",
                "nonbiodegradable_cod",
                "ammonium_as_nitrogen",
                "phosphates",
            ]
        }
    )

    # unit models
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})
    m.fs.mbr_mec = AnaerobicMBRMECZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.vfa_recovery = VFARecoveryZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.gas_sparged_membrane = GasSpargedMembraneZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.ion_exchange = IonExchangeZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "clinoptilolite",
        },
    )
    m.fs.sedimentation = SedimentationZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "phosphorus_capture",
        },
    )
    m.fs.food_waste = FeedZO(default={"property_package": m.fs.prop})

    m.fs.cofermentation = CofermentationZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.constructed_wetlands = ConstructedWetlandsZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.mixer = Mixer(
        default={
            "property_package": m.fs.prop,
            "inlet_list": ["inlet1", "inlet2"],
            "momentum_mixing_type": MomentumMixingType.none,
            "energy_mixing_type": MixingType.none,
        },
    )

    m.fs.product_water = Product(default={"property_package": m.fs.prop})
    m.fs.product_phosphate = Product(default={"property_package": m.fs.prop})
    m.fs.product_ammonia = Product(default={"property_package": m.fs.prop})
    m.fs.product_vfa = Product(default={"property_package": m.fs.prop})
    m.fs.waste_vfa = Product(default={"property_package": m.fs.prop})
    # TODO: because of gas-sparged membrane formulation, H2 gas flow is a unit variable
    #  instead of a mass flow via property model; hence, gas flow exiting the unit is
    #  not connected to a port or state block
    # m.fs.product_hydrogen = Product(default={"property_package": m.fs.prop})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.mbr_mec.inlet)
    m.fs.s02 = Arc(
        source=m.fs.mbr_mec.treated, destination=m.fs.gas_sparged_membrane.inlet
    )
    m.fs.s03 = Arc(
        source=m.fs.mbr_mec.byproduct,
        destination=m.fs.mixer.inlet1,  # Todo: need to add mixer inlet instead of straight to vfa_rec
    )
    m.fs.s04 = Arc(
        source=m.fs.gas_sparged_membrane.treated, destination=m.fs.ion_exchange.inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.gas_sparged_membrane.byproduct,
        destination=m.fs.cofermentation.inlet1,
    )
    m.fs.s06 = Arc(
        source=m.fs.food_waste.outlet, destination=m.fs.cofermentation.inlet2
    )
    m.fs.s07 = Arc(source=m.fs.cofermentation.treated, destination=m.fs.mixer.inlet2)
    m.fs.s08 = Arc(source=m.fs.mixer.outlet, destination=m.fs.vfa_recovery.inlet)
    m.fs.s09 = Arc(source=m.fs.vfa_recovery.treated, destination=m.fs.product_vfa.inlet)
    m.fs.s10 = Arc(source=m.fs.vfa_recovery.byproduct, destination=m.fs.waste_vfa.inlet)
    m.fs.s11 = Arc(
        source=m.fs.ion_exchange.byproduct, destination=m.fs.product_ammonia.inlet
    )
    m.fs.s12 = Arc(
        source=m.fs.ion_exchange.treated, destination=m.fs.sedimentation.inlet
    )
    m.fs.s13 = Arc(
        source=m.fs.sedimentation.byproduct, destination=m.fs.product_phosphate.inlet
    )
    m.fs.s14 = Arc(
        source=m.fs.sedimentation.treated, destination=m.fs.constructed_wetlands.inlet
    )
    m.fs.s15 = Arc(
        source=m.fs.constructed_wetlands.treated, destination=m.fs.product_water.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = pyunits.convert(
        3780 * pyunits.L / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    conc_mass_cod = pyunits.convert(
        2300 * pyunits.mg / pyunits.L, to_units=pyunits.kg / pyunits.m**3
    )
    conc_mass_nh4 = pyunits.convert(
        105 * pyunits.mg / pyunits.L, to_units=pyunits.kg / pyunits.m**3
    )
    conc_mass_po4 = pyunits.convert(
        50 * pyunits.mg / pyunits.L, to_units=pyunits.kg / pyunits.m**3
    )

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "cod"].fix(conc_mass_cod)
    m.fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"].fix(conc_mass_nh4)
    m.fs.feed.conc_mass_comp[0, "phosphates"].fix(conc_mass_po4)
    m.fs.feed.conc_mass_comp[0, "nonbiodegradable_cod"].fix(1e-8)
    solve(m.fs.feed)

    # anaerobic fermentation and MEC
    m.fs.mbr_mec.load_parameters_from_database(use_default_removal=True)
    assert_degrees_of_freedom(m.fs.mbr_mec, expected_dof=5)

    # gas-sparged membrane
    m.fs.gas_sparged_membrane.load_parameters_from_database(use_default_removal=True)
    assert_degrees_of_freedom(m.fs.gas_sparged_membrane, expected_dof=5)

    m.fs.food_waste.flow_vol[0].fix(25e3 * pyunits.L / pyunits.day)
    m.fs.food_waste.conc_mass_comp[0, "cod"].fix(29.23 * pyunits.mg / pyunits.L)
    m.fs.food_waste.conc_mass_comp[0, "ammonium_as_nitrogen"].fix(1e-8)
    m.fs.food_waste.conc_mass_comp[0, "phosphates"].fix(1e-8)
    m.fs.food_waste.conc_mass_comp[0, "nonbiodegradable_cod"].fix(1e-8)
    solve(m.fs.food_waste)

    # cofermentation
    m.fs.cofermentation.load_parameters_from_database(use_default_removal=True)
    assert_degrees_of_freedom(m.fs.cofermentation, expected_dof=10)

    # mixer

    # VFA recovery unit
    m.fs.vfa_recovery.load_parameters_from_database(use_default_removal=True)

    # ion exchange
    m.fs.ion_exchange.load_parameters_from_database(use_default_removal=True)

    # sedimentation
    m.fs.sedimentation.load_parameters_from_database(use_default_removal=True)

    # constructed wetlands
    m.fs.constructed_wetlands.load_parameters_from_database(use_default_removal=True)


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
    unit_list = [
        "feed",
        "mbr_mec",
        "gas_sparged_membrane",
        "ion_exchange",
        "sedimentation",
        "food_waste",
        "cofermentation",
        "mixer",
        "vfa_recovery",
        "constructed_wetlands",
    ]
    for u in unit_list:
        fs.component(u).report()


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "swine_wwt_global_costing.yaml",
    )
    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})
    m.fs.mixer_costing = WaterTAPCosting()
    # typing aid
    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}

    m.fs.mbr_mec.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.vfa_recovery.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.ion_exchange.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.sedimentation.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.cofermentation.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.constructed_wetlands.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.mixer.costing = UnitModelCostingBlock(
        {"default": {"flowsheet_costing_block": m.fs.mixer_costing}}
    )
    m.fs.costing.cost_process()
    m.fs.mixer_costing.cost_process()

    # Resource recovery in terms of annual production
    # Hydrogen gas production
    m.fs.costing.annual_hydrogen_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.gas_sparged_membrane.flow_mass_gas_extraction[0],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    # Phosphorus production
    m.fs.costing.annual_phosphate_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.product_phosphate.flow_mass_comp[0, "phosphates"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    # Nitrogen production
    m.fs.costing.annual_nitrogen_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.product_ammonia.flow_mass_comp[0, "ammonium_as_nitrogen"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    # VFA production
    m.fs.costing.vfa_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.product_vfa.flow_mass_comp[0, "nonbiodegradable_cod"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )
    # Treated effluent production
    m.fs.costing.annual_water_production = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.product_water.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    # Annual COD removed
    m.fs.costing.annual_cod_removed = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.flow_mass_comp[0, "cod"]
            + m.fs.food_waste.flow_mass_comp[0, "cod"]
            - m.fs.product_water.flow_mass_comp[0, "cod"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        )
    )
    # Annual waste sludge released
    m.fs.costing.annual_sludge_waste = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.waste_vfa.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    # Annual costs and revenues
    # Annual water revenue
    m.fs.costing.annual_water_revenue = Expression(
        expr=(
            m.fs.zo_costing.water_product_cost * m.fs.costing.annual_water_production
        ),
        doc="Value of annual water savings, considering treated effluent is utilized to displace freshwater usage",
    )

    # Annual VFA revenue
    m.fs.costing.annual_vfa_revenue = Expression(
        expr=(m.fs.zo_costing.vfa_product_cost * m.fs.costing.annual_vfa_production),
        doc="Annual VFA revenue",
    )

    # Annual Nitrogen revenue
    m.fs.costing.annual_ammonia_revenue = Expression(
        expr=(
            m.fs.zo_costing.ammonia_product_cost
            * m.fs.costing.annual_nitrogen_production
        ),
        doc="Annual ammonia revenue",
    )

    # Annual Phosphorus revenue
    m.fs.costing.annual_phosphorus_revenue = Expression(
        expr=(
            m.fs.zo_costing.phosphorus_product_cost
            * m.fs.costing.annual_phosphate_production
        ),
        doc="Annual phosphorus revenue",
    )

    # Annual Hydrogen revenue
    m.fs.costing.annual_hydrogen_revenue = Expression(
        expr=(
            m.fs.zo_costing.hydrogen_product_cost
            * m.fs.costing.annual_hydrogen_production
        ),
        doc="Annual hydrogen revenue",
    )

    # Annual sludge disposal cost
    m.fs.costing.annual_disposal_cost = Expression(
        expr=(m.fs.zo_costing.waste_disposal_cost * m.fs.costing.annual_sludge_waste),
        doc="Annual sludge disposal cost",
    )
    # Combine results from costing packages and calculate overall metrics
    @m.Expression()
    def total_capital_cost(b):
        return (
            pyunits.convert(
                m.fs.zo_costing.total_capital_cost, to_units=m.fs.costing.base_currency
            )
            + m.fs.mixer_costing.total_investment_cost
        )

    @m.Expression()
    def total_operating_cost(b):
        return (
            pyunits.convert(
                m.fs.zo_costing.total_fixed_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
            + pyunits.convert(
                m.fs.zo_costing.total_variable_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
            + pyunits.convert(
                m.fs.mixer_costing.total_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
            + m.fs.costing.annual_disposal_cost
        )

    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.total_operating_cost
        )
    )
    # Levelized cost of water
    m.fs.costing.LCOW = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_vfa_revenue
            - m.fs.costing.annual_hydrogen_revenue
            - m.fs.costing.annual_ammonia_revenue
            - m.fs.costing.annual_phosphorus_revenue
        )
        / m.fs.costing.annual_water_production,
        doc="Levelized Cost of Treated Water",
    )
    # Levelized cost of hydrogen
    m.fs.costing.LCOH2 = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_vfa_revenue
            - m.fs.costing.annual_water_revenue
            - m.fs.costing.annual_ammonia_revenue
            - m.fs.costing.annual_phosphorus_revenue
        )
        / m.fs.costing.annual_hydrogen_production,
        doc="Levelized Cost of Hydrogen Production",
    )
    # Levelized cost of nitrogen
    m.fs.costing.LCON = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_vfa_revenue
            - m.fs.costing.annual_hydrogen_revenue
            - m.fs.costing.annual_water_revenue
            - m.fs.costing.annual_phosphorus_revenue
        )
        / m.fs.costing.annual_nitrogen_production,
        doc="Levelized Cost of Nitrogen Production",
    )
    # Levelized cost of VFAs
    m.fs.costing.LCOVFA = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_water_revenue
            - m.fs.costing.annual_hydrogen_revenue
            - m.fs.costing.annual_ammonia_revenue
            - m.fs.costing.annual_phosphorus_revenue
        )
        / m.fs.costing.annual_vfa_production,
        doc="Levelized Cost of VFA Production",
    )
    # Levelized cost of phosphorus
    m.fs.costing.LCOP = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_vfa_revenue
            - m.fs.costing.annual_hydrogen_revenue
            - m.fs.costing.annual_ammonia_revenue
            - m.fs.costing.annual_water_revenue
        )
        / m.fs.costing.annual_phosphate_production,
        doc="Levelized Cost of Phosphorus Production",
    )

    # Levelized cost of COD removal
    m.fs.costing.LCOCOD = Expression(
        expr=(
            m.fs.costing.total_annualized_cost
            - m.fs.costing.annual_vfa_revenue
            - m.fs.costing.annual_hydrogen_revenue
            - m.fs.costing.annual_ammonia_revenue
            - m.fs.costing.annual_phosphorus_revenue
            - m.fs.costing.annual_water_revenue
        )
        / m.fs.costing.annual_cod_removed,
        doc="Levelized Cost of COD Removal",
    )

    # TODO add SEC with respect to influent
    # m.fs.costing.add_electricity_intensity(m.fs.product_H2O.properties[0].flow_vol)


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
