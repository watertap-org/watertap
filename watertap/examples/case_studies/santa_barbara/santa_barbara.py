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
                           Constraint,
                           Expression,
                           Objective,
                           Param,
                           TransformationFactory,
                           units as pyunits,
                           assert_optimal_termination,
                           Block)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models.translator import Translator
from idaes.generic_models.unit_models import Mixer, Separator, Product
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.seawater_prop_pack as prop_SW
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                     ConcentrationPolarizationType,
                                                     MassTransferCoefficient,
                                                     PressureChangeType)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pump_isothermal import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (FeedZO,
                                             SWOnshoreIntakeZO,
                                             ChemicalAdditionZO,
                                             ChlorinationZO,
                                             StaticMixerZO,
                                             StorageTankZO,
                                             MediaFiltrationZO,
                                             BackwashSolidsHandlingZO,
                                             CartridgeFiltrationZO,
                                             UVAOPZO,
                                             CO2AdditionZO,
                                             MunicipalDrinkingZO,
                                             LandfillZO)



def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    solve(m)
    # initialize_system(m, solver=solver)
    m.fs.display()
    assert False
    #
    # # simulate
    # solve(m, solver=solver)
    # print('\n***---Simulation results---***')
    # display_system(m)
    # display_design(m)
    # display_state(m)

def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.prop_prtrt = prop_ZO.WaterParameterBlock(default={"solute_list": ["tds", "tss"]})
    density = 1023.5 * pyunits.kg / pyunits.m ** 3
    m.fs.prop_prtrt.dens_mass_default = density
    m.fs.prop_psttrt = prop_ZO.WaterParameterBlock(default={"solute_list": ["tds"]})
    m.fs.prop_SW = prop_SW.SeawaterParameterBlock()

    # block structure
    prtrt = m.fs.pretreatment = Block()
    desal = m.fs.desalination = Block()
    psttrt = m.fs.posttreatment = Block()

    # unit models
    m.fs.feed = FeedZO(default={'property_package': m.fs.prop_prtrt})
    # pretreatment
    prtrt.intake = SWOnshoreIntakeZO(default={'property_package': m.fs.prop_prtrt})
    prtrt.ferric_chloride_addition = ChemicalAdditionZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db,
            "process_subtype": "ferric_chloride"})
    prtrt.chlorination = ChlorinationZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db})
    prtrt.static_mixer = StaticMixerZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db})
    prtrt.storage_tank_1 = StorageTankZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db})
    prtrt.media_filtration = MediaFiltrationZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db})
    prtrt.backwash_handling = BackwashSolidsHandlingZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db})
    prtrt.anti_scalant_addition = ChemicalAdditionZO(default={
            "property_package": m.fs.prop_prtrt,
            "database": m.db,
            "process_subtype": "anti-scalant"})
    prtrt.cartridge_filtration = CartridgeFiltrationZO(default={
        "property_package": m.fs.prop_prtrt,
        "database": m.db})

    # desalination
    # m.fs.S1 = Separator(default={
    #     "property_package": m.fs.properties,
    #     "outlet_list": ['P1', 'PXR']})
    # m.fs.P1 = Pump(default={'property_package': m.fs.properties})
    # m.fs.PXR = PressureExchanger(default={'property_package': m.fs.properties})
    # m.fs.P2 = Pump(default={'property_package': m.fs.properties})
    # m.fs.M1 = Mixer(default={
    #     "property_package": m.fs.properties,
    #     "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
    #     "inlet_list": ['P1', 'P2']})
    # m.fs.RO = ReverseOsmosis0D(default={
    #     "property_package": m.fs.properties,
    #     "has_pressure_change": True,
    #     "pressure_change_type": PressureChangeType.calculated,
    #     "mass_transfer_coefficient": MassTransferCoefficient.calculated,
    #     "concentration_polarization_type": ConcentrationPolarizationType.calculated,
    # })

    # posttreatment
    psttrt.storage_tank_2 = StorageTankZO(default={
            "property_package": m.fs.prop_psttrt,
            "database": m.db})
    psttrt.uv_aop = UVAOPZO(default={
        "property_package": m.fs.prop_psttrt,
        "database": m.db,
        "process_subtype": "hydrogen_peroxide"})
    psttrt.co2_addition = CO2AdditionZO(default={
        "property_package": m.fs.prop_psttrt,
        "database": m.db})
    psttrt.lime_addition = ChemicalAdditionZO(default={
            "property_package": m.fs.prop_psttrt,
            "database": m.db,
            "process_subtype": "lime"})
    psttrt.storage_tank_3 = StorageTankZO(default={
        "property_package": m.fs.prop_psttrt,
        "database": m.db})

    # product and disposal
    m.fs.municipal = MunicipalDrinkingZO(default={
        "property_package": m.fs.prop_psttrt,
        "database": m.db})
    m.fs.landfill = LandfillZO(default={
        "property_package": m.fs.prop_prtrt,
        "database": m.db})

    # translator blocks
    m.fs.tb_prtrt_psttrt = Translator(
        default={"inlet_property_package": m.fs.prop_prtrt,
                 "outlet_property_package": m.fs.prop_psttrt})

    @m.fs.tb_prtrt_psttrt.Constraint(['H2O', 'tds'])
    def eq_flow_mass_comp(blk, j):
        return blk.properties_in[0].flow_mass_comp[j] == blk.properties_out[0].flow_mass_comp[j]

    # connections
    prtrt.s01 = Arc(source=m.fs.feed.outlet, destination=prtrt.intake.inlet)
    prtrt.s02 = Arc(source=prtrt.intake.outlet, destination=prtrt.ferric_chloride_addition.inlet)
    prtrt.s03 = Arc(source=prtrt.ferric_chloride_addition.outlet, destination=prtrt.chlorination.inlet)
    prtrt.s04 = Arc(source=prtrt.chlorination.treated, destination=prtrt.static_mixer.inlet)
    prtrt.s05 = Arc(source=prtrt.static_mixer.outlet, destination=prtrt.storage_tank_1.inlet)
    prtrt.s06 = Arc(source=prtrt.storage_tank_1.outlet, destination=prtrt.media_filtration.inlet)
    prtrt.s07 = Arc(source=prtrt.media_filtration.byproduct, destination=prtrt.backwash_handling.inlet)
    prtrt.s08 = Arc(source=prtrt.media_filtration.treated, destination=prtrt.anti_scalant_addition.inlet)
    prtrt.s09 = Arc(source=prtrt.anti_scalant_addition.outlet, destination=prtrt.cartridge_filtration.inlet)
    prtrt.s10 = Arc(source=prtrt.cartridge_filtration.treated, destination=m.fs.tb_prtrt_psttrt.inlet)
    prtrt.s11 = Arc(source=prtrt.backwash_handling.byproduct, destination=m.fs.landfill.inlet)

    psttrt.s01 = Arc(source=m.fs.tb_prtrt_psttrt.outlet, destination=psttrt.storage_tank_2.inlet)
    psttrt.s02 = Arc(source=psttrt.storage_tank_2.outlet, destination=psttrt.uv_aop.inlet)
    psttrt.s03 = Arc(source=psttrt.uv_aop.treated, destination=psttrt.co2_addition.inlet)
    psttrt.s04 = Arc(source=psttrt.co2_addition.outlet, destination=psttrt.lime_addition.inlet)
    psttrt.s05 = Arc(source=psttrt.lime_addition.outlet, destination=psttrt.storage_tank_3.inlet)
    psttrt.s06 = Arc(source=psttrt.storage_tank_3.outlet, destination=m.fs.municipal.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling

    # calculate and propagate scaling factors
    # iscale.calculate_scaling_factors(m)

    return m

def set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=None):
    if solver is None:
        solver = get_solver()

    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # ---specifications---
    # feed
    flow_vol = 0.3092 * pyunits.m**3/pyunits.s
    conc_mass_tds = 35 * pyunits.kg/pyunits.m**3
    conc_mass_tss = 0.03 * pyunits.kg/pyunits.m**3
    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)

    # ---pretreatment---
    # intake

    # ferric chloride
    m.db.get_unit_operation_parameters("chemical_addition")
    prtrt.ferric_chloride_addition.load_parameters_from_database()
    prtrt.ferric_chloride_addition.chemical_dosage.fix(20)

    # chlorination
    m.db.get_unit_operation_parameters("chlorination")
    prtrt.chlorination.load_parameters_from_database(use_default_removal=True)

    # static mixer
    m.db.get_unit_operation_parameters("static_mixer")
    prtrt.static_mixer.load_parameters_from_database(use_default_removal=True)

    # storage tank
    m.db.get_unit_operation_parameters("storage_tank")
    prtrt.storage_tank_1.load_parameters_from_database(use_default_removal=True)
    prtrt.storage_tank_1.storage_time.fix(2)

    # media filtration
    m.db.get_unit_operation_parameters("media_filtration")
    prtrt.media_filtration.load_parameters_from_database(use_default_removal=True)

    # backwash handling
    m.db.get_unit_operation_parameters("backwash_solids_handling")
    prtrt.backwash_handling.load_parameters_from_database(use_default_removal=True)

    # anti-scalant
    prtrt.anti_scalant_addition.load_parameters_from_database()
    prtrt.anti_scalant_addition.chemical_dosage.fix(5)

    # cartridge filtration
    m.db.get_unit_operation_parameters("cartridge_filtration")
    prtrt.cartridge_filtration.load_parameters_from_database(use_default_removal=True)

    # ---desalination---

    # ---posttreatment---
    # storage tank 2
    psttrt.storage_tank_2.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_2.storage_time.fix(1)

    # uv aop
    m.db.get_unit_operation_parameters("uv_aop")
    psttrt.uv_aop.load_parameters_from_database(use_default_removal=True)
    psttrt.uv_aop.uv_reduced_equivalent_dose.fix(350)  # TODO: check this was the right thing to fix
    psttrt.uv_aop.uv_transmittance_in.fix(0.95)  # TODO: check this was the right thing to fix

    # uv aop
    m.db.get_unit_operation_parameters("co2_addition")
    psttrt.co2_addition.load_parameters_from_database(use_default_removal=True)

    # lime
    psttrt.lime_addition.load_parameters_from_database()
    psttrt.lime_addition.chemical_dosage.fix(2.3)

    # storage tank 3
    psttrt.storage_tank_3.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_3.storage_time.fix(1)

    # ---product and disposal---
    m.db.get_unit_operation_parameters("municipal_drinking")
    m.fs.municipal.load_parameters_from_database()

    m.db.get_unit_operation_parameters("landfill")
    m.fs.landfill.load_parameters_from_database()


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # ---initialize feed block---
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    m.fs.intake.initialize(optarg=optarg)


if __name__ == "__main__":
    main()
