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
from pyomo.environ import (ConcreteModel,
                           value,
                           TransformationFactory,
                           units as pyunits,
                           assert_optimal_termination,
                           Block)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import (propagate_state,
                                            fix_state_vars,
                                            revert_state_vars)
from idaes.core.util.exceptions import ConfigurationError
from idaes.generic_models.unit_models.translator import Translator
from idaes.generic_models.unit_models import Mixer, Separator, Product
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.util.infeasible import print_infeasible_bounds, print_infeasible_constraints, print_close_to_bounds
from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (FeedZO,
                                             MunicipalDrinkingZO,
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
                                             BackwashSolidsHandlingZO)
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)

    initialize_system(m)

    results = solve(m, check_termination=False)  # initialization not needed for flowsheet of only zero order models
    display_results(m)
    # m.fs.ozonation.display()
    # print_infeasible_constraints(m)
    # print_infeasible_bounds(m)
    # print_close_to_bounds(m)
    # assert_optimal_termination(results)

    #
    # add_costing(m)
    # initialize_costing(m)
    # assert_degrees_of_freedom(m, 0)
    #
    # solve(m, tee=True)
    # display_costing(m)
    return m


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(default={"solute_list": ["tds", "tss", "toc"]})

    # unit models
    m.fs.feed = FeedZO(default={'property_package': m.fs.prop})
    m.fs.ozonation = OzoneZO(default={
        "property_package": m.fs.prop,
        "database": m.db})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.ozonation.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    # set unit model values
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def set_operating_conditions(m):

    # ---specifications---
    # feed
    flow_vol = 0.9224 * pyunits.m**3/pyunits.s
    conc_mass_tds = 0.63 * pyunits.kg/pyunits.m**3
    conc_mass_tss = 0.006525 * pyunits.kg/pyunits.m**3
    conc_mass_toc = 0.004 * pyunits.kg / pyunits.m ** 3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "toc"].fix(conc_mass_toc)
    solve(m.fs.feed)

    # # ozonation
    m.db.get_unit_operation_parameters("ozonation")
    m.fs.ozonation.load_parameters_from_database(use_default_removal=True)
    # m.fs.ozonation.display()
    # for j in ['H2O', 'tds', 'toc', 'tss']:
    #     m.fs.ozonation.properties_in[0].flow_mass_comp[j].fix(m.fs.feed.properties[0].flow_mass_comp[j])


def initialize_system(m):
    solve(m.fs.feed)
    m.fs.ozonation.properties_in[0].display()
    propagate_state(m.fs.s01)
    m.fs.ozonation.properties_in[0].display()
    assert False
    # flags = fix_state_vars(m.fs.ozonation.properties_in)
    # solve(m.fs.ozonation)
    # revert_state_vars(m.fs.ozonation.properties_in, flags)


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display_results(m):
    unit_list = ['feed', 'intake_pump', 'coag_and_flow',
                 'sedimentation', 'ozonation', 'gravity_basin',
                 'gac', 'backwash_pump', 'uv', 'anion_exchange',
                 'chlorination', 'storage', 'recharge_pump']

    for u in unit_list:
        if hasattr(m.fs, u):
            unit = getattr(m.fs, u)
            unit.report()


if __name__ == "__main__":
    main()
