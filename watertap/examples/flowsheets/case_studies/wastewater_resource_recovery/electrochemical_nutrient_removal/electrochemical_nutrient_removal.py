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
import idaes.core.util.scaling as iscale
from idaes.generic_models.unit_models import Product
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.unit_models.zero_order import (
    FeedZO,
    ElectroNPZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting


def main():
    m = build()
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)
    initialize_system(m)

    solve(m)

    display_results(m)

    add_costing(m)

    solve(m)
    # display_costing(m)
    return m


def build():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={"solute_list": ["nitrogen", "phosphorus", "struvite"]}
    )

    # components
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})
    m.fs.product_H2O = Product(default={"property_package": m.fs.prop})
    m.fs.product_struvite = Product(default={"property_package": m.fs.prop})
    m.fs.electroNP = ElectroNPZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        }
    )

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.electroNP.inlet)
    m.fs.s02 = Arc(source=m.fs.electroNP.treated, destination=m.fs.product_H2O.inlet)
    m.fs.s03 = Arc(
        source=m.fs.electroNP.byproduct, destination=m.fs.product_struvite.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
    # feed
    flow_vol = 37.9 / 3600 * pyunits.m**3 / pyunits.s
    conc_nitrogen = 0.715 * pyunits.kg / pyunits.m**3
    conc_phosphorus = 0.715 * pyunits.kg / pyunits.m**3
    conc_struvite = 0 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "nitrogen"].fix(conc_nitrogen)
    m.fs.feed.conc_mass_comp[0, "phosphorus"].fix(conc_phosphorus)
    m.fs.feed.conc_mass_comp[0, "struvite"].fix(conc_struvite)
    solve(m.fs.feed)

    # electroNP
    m.fs.electroNP.load_parameters_from_database(use_default_removal=True)


def add_costing(m):
    # yaml file path
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "case_1617.yaml",
    )

    # add costing block
    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})
    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.electroNP.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    # TODO- verify flow basis for electricity intensity and LCOW
    m.fs.costing.add_electricity_intensity(m.fs.product_struvite.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product_struvite.properties[0].flow_vol)

    # initialize
    m.fs.costing.initialize()


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
    unit_list = ["feed", "product_H2O", "product_struvite", "electroNP"]
    for u in unit_list:
        m.fs.component(u).report()


if __name__ == "__main__":
    m = main()
