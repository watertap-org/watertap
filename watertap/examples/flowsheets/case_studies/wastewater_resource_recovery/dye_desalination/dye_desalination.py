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
    NanofiltrationZO,
)
from watertap.core.zero_order_costing import ZeroOrderCostingData as ZeroOrderCosting


def main():
    m = build()
    return m


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(default={"solute_list": ["dye", "tds"]})

    # unit model
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})

    m.fs.nanofiltration = NanofiltrationZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": "rHGO_dye_rejection",
        }
    )

    m.fs.permeate = Product(default={"property_package": m.fs.prop})
    m.fs.retentate = Product(default={"property_package": m.fs.prop})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.nanofiltration.inlet)
    m.fs.s02 = Arc(source=m.fs.nanofiltration.treated, destination=m.fs.permeate.inlet)
    m.fs.s03 = Arc(
        source=m.fs.nanofiltration.byproduct, destination=m.fs.retentate.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


if __name__ == "__main__":
    model = main()
