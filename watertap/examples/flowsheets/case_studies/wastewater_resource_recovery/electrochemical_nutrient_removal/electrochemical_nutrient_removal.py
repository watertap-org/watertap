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

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    ElectroNPZO,
)


def main():
    m = build()
    # set_operating_conditions(m)
    # assert_degrees_of_freedom(m,0)
    assert_units_consistent(m)
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
    m.fs.permeate = Product(default={"property_package": m.fs.prop})
    m.fs.retentate = Product(default={"property_package": m.fs.prop})
    m.fs.electroNP = ElectroNPZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        }
    )

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.electroNP.inlet)
    m.fs.s02 = Arc(source=m.fs.electroNP.treated, destination=m.fs.permeate.inlet)
    m.fs.s03 = Arc(source=m.fs.electroNP.byproduct, destination=m.fs.retentate.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
    # feed
    flow_vol = 37.9 / 3600 * pyunits.m**3 / pyunits.s
    conc_n = 715 * pyunits.kg / pyunits.m**3
    conc_p = 715 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "nitrogen"].fix(conc_n)
    m.fs.feed.conc_mass_comp[0, "phosphorus"].fix(conc_p)
    solve(m.fs.feed)

    # electroNP
    m.fs.electroNP.load_parameters_from_database(use_default_removal=True)


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


if __name__ == "__main__":
    m = main()
