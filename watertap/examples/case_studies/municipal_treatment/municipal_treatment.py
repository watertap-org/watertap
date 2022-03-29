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
                                             BackwashSolidsHandlingZO)
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)

    initialize_system(m)  # initialization needed for ozone unit

    results = solve(m)
    display_results(m)

    # TODO: add costing, currently costing functions only pass
    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)

    # results = solve(m)
    display_costing(m)
    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(default={"solute_list": ["tds", "tss", "toc"]})

    # unit models
    m.fs.feed = FeedZO(default={'property_package': m.fs.prop})
    m.fs.intake_pump = WaterPumpingStationZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "raw"})
    m.fs.coag_and_floc = CoagulationFlocculationZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.sedimentation = SedimentationZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.ozonation = OzoneZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.gravity_basin = FixedBedZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "gravity_basin"})
    m.fs.gac = GACZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "pressure_vessel"})
    m.fs.backwash_pump = WaterPumpingStationZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "treated"})
    m.fs.uv = UVZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.anion_exchange = IonExchangeZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "anion_exchange"})
    m.fs.chlorination = ChlorinationZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.storage = StorageTankZO(default={
        "property_package": m.fs.prop,
        "database": m.db})
    m.fs.recharge_pump = WaterPumpingStationZO(default={
        "property_package": m.fs.prop,
        "database": m.db,
        "process_subtype": "treated"})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.intake_pump.inlet)
    m.fs.s02 = Arc(source=m.fs.intake_pump.outlet, destination=m.fs.coag_and_floc.inlet)
    m.fs.s03 = Arc(source=m.fs.coag_and_floc.outlet, destination=m.fs.sedimentation.inlet)
    m.fs.s04 = Arc(source=m.fs.sedimentation.treated, destination=m.fs.ozonation.inlet)
    m.fs.s05 = Arc(source=m.fs.ozonation.treated, destination=m.fs.gravity_basin.inlet)
    m.fs.s06 = Arc(source=m.fs.gravity_basin.treated, destination=m.fs.gac.inlet)
    m.fs.s07 = Arc(source=m.fs.gac.treated, destination=m.fs.uv.inlet)
    m.fs.s08 = Arc(source=m.fs.gac.byproduct, destination=m.fs.backwash_pump.inlet)
    m.fs.s09 = Arc(source=m.fs.uv.treated, destination=m.fs.anion_exchange.inlet)
    m.fs.s10 = Arc(source=m.fs.anion_exchange.treated, destination=m.fs.chlorination.inlet)
    m.fs.s11 = Arc(source=m.fs.chlorination.treated, destination=m.fs.storage.inlet)
    m.fs.s12 = Arc(source=m.fs.storage.outlet, destination=m.fs.recharge_pump.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
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

    # chlorination
    m.fs.chlorination.load_parameters_from_database(use_default_removal=True)

    # storage
    m.fs.storage.load_parameters_from_database(use_default_removal=True)
    m.fs.storage.storage_time.fix(6)

    # recharge pump
    m.fs.recharge_pump.load_parameters_from_database()
    m.fs.recharge_pump.electricity.fix(186.4)


def initialize_system(m):
    unit_list = ['feed', 'intake_pump', 'coag_and_floc',
                 'sedimentation', 'ozonation', 'gravity_basin',
                 'gac', 'backwash_pump', 'uv', 'anion_exchange',
                 'chlorination', 'storage', 'recharge_pump']
    unit_inlet_arc_dic = {'feed': None,
                          'intake_pump': 's01',
                          'coag_and_floc': 's02',
                          'sedimentation': 's03',
                          'ozonation': 's04',
                          'gravity_basin': 's05',
                          'gac': 's06',
                          'uv': 's07',
                          'backwash_pump': 's08',
                          'anion_exchange': 's09',
                          'chlorination': 's10',
                          'storage': 's11',
                          'recharge_pump': 's12'}
    unit_inlet_properties_dic = {'feed': 'properties',
                          'intake_pump': 'properties',
                          'coag_and_floc': 'properties',
                          'sedimentation': 'properties_in',
                          'ozonation': 'properties_in',
                          'gravity_basin': 'properties_in',
                          'gac': 'properties_in',
                          'uv': 'properties_in',
                          'backwash_pump': 'properties',
                          'anion_exchange': 'properties_in',
                          'chlorination': 'properties_in',
                          'storage': 'properties',
                          'recharge_pump': 'properties'}

    for u_str in unit_list:
        if unit_inlet_arc_dic[u_str] is not None:
            unit = getattr(m.fs, u_str)
            arc = getattr(m.fs, unit_inlet_arc_dic[u_str])
            inlet_blk = getattr(unit, unit_inlet_properties_dic[u_str])
            propagate_state(arc)
            if unit_inlet_properties_dic[u_str] == 'properties_in':
                flags = fix_state_vars(inlet_blk)
                solve(unit)
                revert_state_vars(inlet_blk, flags)


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display_results(m):
    unit_list = ['feed', 'intake_pump', 'coag_and_floc',
                 'sedimentation', 'ozonation', 'gravity_basin',
                 'gac', 'backwash_pump', 'uv', 'anion_exchange',
                 'chlorination', 'storage', 'recharge_pump']

    for u in unit_list:
        if hasattr(m.fs, u):
            unit = getattr(m.fs, u)
            unit.report()


def add_costing(m):
    # TODO: add costing
    pass


def initialize_costing(m):
    # TODO: initialize costing
    pass


def display_costing(m):
     # TODO: display costing
    pass


if __name__ == "__main__":
    main()
