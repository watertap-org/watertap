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
    MunicipalDrinkingZO,
    WaterPumpingStationZO,
    PumpZO,
    CoagulationFlocculationZO,
    SedimentationZO,
    OzoneZO,
    FixedBedZO,
    GACZO,
    UVZO,
    IonExchangeZO,
    ChlorinationZO,
    StorageTankZO,
    BackwashSolidsHandlingZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)

    initialize_system(m)  # initialization needed for ozone unit

    results = solve(m)
    display_results(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    results = solve(m)
    display_costing(m)
    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(solute_list=["tds", "tss", "toc"])

    # unit models
    m.fs.feed = FeedZO(property_package=m.fs.prop)
    m.fs.intake_pump = WaterPumpingStationZO(
        property_package=m.fs.prop, database=m.db, process_subtype="raw"
    )
    m.fs.coag_and_floc = CoagulationFlocculationZO(
        property_package=m.fs.prop, database=m.db
    )
    m.fs.sedimentation = SedimentationZO(property_package=m.fs.prop, database=m.db)
    m.fs.ozonation = OzoneZO(property_package=m.fs.prop, database=m.db)
    m.fs.gravity_basin = FixedBedZO(
        property_package=m.fs.prop, database=m.db, process_subtype="gravity_basin"
    )
    m.fs.gac = GACZO(
        property_package=m.fs.prop, database=m.db, process_subtype="pressure_vessel"
    )
    m.fs.backwash_pump = WaterPumpingStationZO(
        property_package=m.fs.prop, database=m.db, process_subtype="treated"
    )
    m.fs.uv = UVZO(property_package=m.fs.prop, database=m.db)
    m.fs.anion_exchange = IonExchangeZO(
        property_package=m.fs.prop, database=m.db, process_subtype="anion_exchange"
    )
    m.fs.chlorination = ChlorinationZO(property_package=m.fs.prop, database=m.db)
    m.fs.storage = StorageTankZO(property_package=m.fs.prop, database=m.db)
    m.fs.recharge_pump = WaterPumpingStationZO(
        property_package=m.fs.prop, database=m.db, process_subtype="treated"
    )
    m.fs.product = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.intake_pump.inlet)
    m.fs.s02 = Arc(source=m.fs.intake_pump.outlet, destination=m.fs.coag_and_floc.inlet)
    m.fs.s03 = Arc(
        source=m.fs.coag_and_floc.outlet, destination=m.fs.sedimentation.inlet
    )
    m.fs.s04 = Arc(source=m.fs.sedimentation.treated, destination=m.fs.ozonation.inlet)
    m.fs.s05 = Arc(source=m.fs.ozonation.treated, destination=m.fs.gravity_basin.inlet)
    m.fs.s06 = Arc(source=m.fs.gravity_basin.treated, destination=m.fs.gac.inlet)
    m.fs.s07 = Arc(source=m.fs.gac.treated, destination=m.fs.uv.inlet)
    m.fs.s08 = Arc(source=m.fs.gac.byproduct, destination=m.fs.backwash_pump.inlet)
    m.fs.s09 = Arc(source=m.fs.uv.treated, destination=m.fs.anion_exchange.inlet)
    m.fs.s10 = Arc(
        source=m.fs.anion_exchange.treated, destination=m.fs.chlorination.inlet
    )
    m.fs.s11 = Arc(source=m.fs.chlorination.treated, destination=m.fs.storage.inlet)
    m.fs.s12 = Arc(source=m.fs.storage.outlet, destination=m.fs.recharge_pump.inlet)
    m.fs.s13 = Arc(source=m.fs.recharge_pump.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 0.9224 * pyunits.m**3 / pyunits.s
    conc_mass_tds = 0.63 * pyunits.kg / pyunits.m**3
    conc_mass_tss = 0.006525 * pyunits.kg / pyunits.m**3
    conc_mass_toc = 0.004 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "toc"].fix(conc_mass_toc)
    solve(m.fs.feed)

    # intake pump
    m.fs.intake_pump.load_parameters_from_database()
    m.fs.intake_pump.electricity.fix(93.2)

    # coagulation and flocculation
    m.fs.coag_and_floc.load_parameters_from_database(use_default_removal=True)

    # sedimentation
    m.fs.sedimentation.load_parameters_from_database(use_default_removal=True)

    # # ozonation
    m.fs.ozonation.load_parameters_from_database(use_default_removal=True)

    # fixed bed gravity basin
    m.fs.gravity_basin.load_parameters_from_database(use_default_removal=True)

    # granular activated carbon
    m.fs.gac.load_parameters_from_database(use_default_removal=True)

    # backwash pump
    m.fs.backwash_pump.load_parameters_from_database()
    m.fs.backwash_pump.electricity.fix(37.3)

    # uv aop
    m.fs.uv.load_parameters_from_database(use_default_removal=True)
    m.fs.uv.uv_reduced_equivalent_dose.fix(200)
    m.fs.uv.uv_transmittance_in.fix(0.90)

    # anion exchange
    m.fs.anion_exchange.load_parameters_from_database(use_default_removal=True)
    m.fs.anion_exchange.removal_frac_mass_comp[0, "tds"].fix(0.9)

    # chlorination
    m.fs.chlorination.load_parameters_from_database(use_default_removal=True)

    # storage
    m.fs.storage.load_parameters_from_database(use_default_removal=True)
    m.fs.storage.storage_time.fix(6)

    # recharge pump
    m.fs.recharge_pump.load_parameters_from_database()
    m.fs.recharge_pump.electricity.fix(186.4)


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
    unit_list = [
        "feed",
        "intake_pump",
        "coag_and_floc",
        "sedimentation",
        "ozonation",
        "gravity_basin",
        "gac",
        "backwash_pump",
        "uv",
        "anion_exchange",
        "chlorination",
        "storage",
        "recharge_pump",
        "product",
    ]

    for u in unit_list:
        m.fs.component(u).report()


def add_costing(m):
    m.fs.costing = ZeroOrderCosting()
    # typing aid
    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.intake_pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.coag_and_floc.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.sedimentation.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.ozonation.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.gravity_basin.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.gac.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.backwash_pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.uv.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.anion_exchange.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.chlorination.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.storage.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.recharge_pump.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product.properties[0].flow_vol)
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
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018)),
        )

    print("\nUtility Costs\n")
    for f in m.fs.costing.flow_types:
        print(
            f,
            " :   ",
            value(
                pyunits.convert(
                    m.fs.costing.aggregate_flow_costs[f],
                    to_units=pyunits.USD_2018 / pyunits.year,
                )
            ),
        )

    print("")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2018)
    )
    print(f"Total Capital Costs: {total_capital_cost:.2f} M$")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.2f} M$/year")
    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m^3")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m^3")


if __name__ == "__main__":
    m, results = main()
