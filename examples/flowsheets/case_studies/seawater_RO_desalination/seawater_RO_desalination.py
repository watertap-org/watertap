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
    Block,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import (
    propagate_state,
    fix_state_vars,
    revert_state_vars,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models import Mixer, Separator, Product
from idaes.models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import UnitModelCostingBlock

import watertap.property_models.seawater_prop_pack as prop_SW
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
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
    LandfillZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


def build_flowsheet(erd_type=None):
    m = build(erd_type=erd_type)
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    return m


def solve_flowsheet(flowsheet=None):
    m = flowsheet.parent_block()  # UI block is 'm.fs' but funcs below use 'm'
    initialize_system(m)
    assert_degrees_of_freedom(m, 0)
    solve(m)
    display_results(m)
    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)
    solve(m)


def main(erd_type="pressure_exchanger"):
    m = build_flowsheet(erd_type=erd_type)
    # m = build(erd_type=erd_type)
    #
    # set_operating_conditions(m)
    # assert_degrees_of_freedom(m, 0)

    initialize_system(m)
    assert_degrees_of_freedom(m, 0)

    solve(m)
    display_results(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)

    solve(m, tee=True)
    display_costing(m)

    return m


def build(erd_type=None):
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()
    m.erd_type = erd_type

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_prtrt = prop_ZO.WaterParameterBlock(solute_list=["tds", "tss"])
    density = 1023.5 * pyunits.kg / pyunits.m**3
    m.fs.prop_prtrt.dens_mass_default = density
    m.fs.prop_psttrt = prop_ZO.WaterParameterBlock(solute_list=["tds"])
    m.fs.prop_desal = prop_SW.SeawaterParameterBlock()

    # block structure
    prtrt = m.fs.pretreatment = Block()
    desal = m.fs.desalination = Block()
    psttrt = m.fs.posttreatment = Block()

    # unit models
    m.fs.feed = FeedZO(property_package=m.fs.prop_prtrt)
    # pretreatment
    prtrt.intake = SWOnshoreIntakeZO(property_package=m.fs.prop_prtrt)
    prtrt.ferric_chloride_addition = ChemicalAdditionZO(
        property_package=m.fs.prop_prtrt,
        database=m.db,
        process_subtype="ferric_chloride",
    )
    prtrt.chlorination = ChlorinationZO(property_package=m.fs.prop_prtrt, database=m.db)
    prtrt.static_mixer = StaticMixerZO(property_package=m.fs.prop_prtrt, database=m.db)
    prtrt.storage_tank_1 = StorageTankZO(
        property_package=m.fs.prop_prtrt, database=m.db
    )
    prtrt.media_filtration = MediaFiltrationZO(
        property_package=m.fs.prop_prtrt, database=m.db
    )
    prtrt.backwash_handling = BackwashSolidsHandlingZO(
        property_package=m.fs.prop_prtrt, database=m.db
    )
    prtrt.anti_scalant_addition = ChemicalAdditionZO(
        property_package=m.fs.prop_prtrt, database=m.db, process_subtype="anti-scalant"
    )
    prtrt.cartridge_filtration = CartridgeFiltrationZO(
        property_package=m.fs.prop_prtrt, database=m.db
    )

    # desalination
    desal.P1 = Pump(property_package=m.fs.prop_desal)
    desal.RO = ReverseOsmosis0D(
        property_package=m.fs.prop_desal,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    desal.RO.width.setub(5000)
    desal.RO.area.setub(20000)
    if erd_type == "pressure_exchanger":
        desal.S1 = Separator(
            property_package=m.fs.prop_desal, outlet_list=["P1", "PXR"]
        )
        desal.M1 = Mixer(
            property_package=m.fs.prop_desal,
            momentum_mixing_type=MomentumMixingType.equality,
            inlet_list=["P1", "P2"],
        )
        desal.PXR = PressureExchanger(property_package=m.fs.prop_desal)
        desal.P2 = Pump(property_package=m.fs.prop_desal)
    elif erd_type == "pump_as_turbine":
        desal.ERD = EnergyRecoveryDevice(property_package=m.fs.prop_desal)
    else:
        raise ConfigurationError(
            "erd_type was {}, but can only "
            "be pressure_exchanger or pump_as_turbine"
            "".format(erd_type)
        )

    # posttreatment
    psttrt.storage_tank_2 = StorageTankZO(
        property_package=m.fs.prop_psttrt, database=m.db
    )
    psttrt.uv_aop = UVAOPZO(
        property_package=m.fs.prop_psttrt,
        database=m.db,
        process_subtype="hydrogen_peroxide",
    )
    psttrt.co2_addition = CO2AdditionZO(
        property_package=m.fs.prop_psttrt, database=m.db
    )
    psttrt.lime_addition = ChemicalAdditionZO(
        property_package=m.fs.prop_psttrt, database=m.db, process_subtype="lime"
    )
    psttrt.storage_tank_3 = StorageTankZO(
        property_package=m.fs.prop_psttrt, database=m.db
    )

    # product and disposal
    m.fs.municipal = MunicipalDrinkingZO(
        property_package=m.fs.prop_psttrt, database=m.db
    )
    m.fs.landfill = LandfillZO(property_package=m.fs.prop_prtrt, database=m.db)
    m.fs.disposal = Product(property_package=m.fs.prop_desal)

    # translator blocks
    m.fs.tb_prtrt_desal = Translator(
        inlet_property_package=m.fs.prop_prtrt, outlet_property_package=m.fs.prop_desal
    )

    @m.fs.tb_prtrt_desal.Constraint(["H2O", "tds"])
    def eq_flow_mass_comp(blk, j):
        if j == "tds":
            jj = "TDS"
        else:
            jj = j
        return (
            blk.properties_in[0].flow_mass_comp[j]
            == blk.properties_out[0].flow_mass_phase_comp["Liq", jj]
        )

    m.fs.tb_desal_psttrt = Translator(
        inlet_property_package=m.fs.prop_desal, outlet_property_package=m.fs.prop_psttrt
    )

    @m.fs.tb_desal_psttrt.Constraint(["H2O", "TDS"])
    def eq_flow_mass_comp(blk, j):
        if j == "TDS":
            jj = "tds"
        else:
            jj = j
        return (
            blk.properties_in[0].flow_mass_phase_comp["Liq", j]
            == blk.properties_out[0].flow_mass_comp[jj]
        )

    # connections
    m.fs.s_feed = Arc(source=m.fs.feed.outlet, destination=prtrt.intake.inlet)
    prtrt.s01 = Arc(
        source=prtrt.intake.outlet, destination=prtrt.ferric_chloride_addition.inlet
    )
    prtrt.s02 = Arc(
        source=prtrt.ferric_chloride_addition.outlet,
        destination=prtrt.chlorination.inlet,
    )
    prtrt.s03 = Arc(
        source=prtrt.chlorination.treated, destination=prtrt.static_mixer.inlet
    )
    prtrt.s04 = Arc(
        source=prtrt.static_mixer.outlet, destination=prtrt.storage_tank_1.inlet
    )
    prtrt.s05 = Arc(
        source=prtrt.storage_tank_1.outlet, destination=prtrt.media_filtration.inlet
    )
    prtrt.s06 = Arc(
        source=prtrt.media_filtration.byproduct,
        destination=prtrt.backwash_handling.inlet,
    )
    prtrt.s07 = Arc(
        source=prtrt.media_filtration.treated,
        destination=prtrt.anti_scalant_addition.inlet,
    )
    prtrt.s08 = Arc(
        source=prtrt.anti_scalant_addition.outlet,
        destination=prtrt.cartridge_filtration.inlet,
    )
    m.fs.s_prtrt_tb = Arc(
        source=prtrt.cartridge_filtration.treated, destination=m.fs.tb_prtrt_desal.inlet
    )
    m.fs.s_landfill = Arc(
        source=prtrt.backwash_handling.byproduct, destination=m.fs.landfill.inlet
    )

    if erd_type == "pressure_exchanger":
        m.fs.s_tb_desal = Arc(
            source=m.fs.tb_prtrt_desal.outlet, destination=desal.S1.inlet
        )
        desal.s01 = Arc(source=desal.S1.P1, destination=desal.P1.inlet)
        desal.s02 = Arc(source=desal.P1.outlet, destination=desal.M1.P1)
        desal.s03 = Arc(source=desal.M1.outlet, destination=desal.RO.inlet)
        desal.s04 = Arc(
            source=desal.RO.retentate, destination=desal.PXR.high_pressure_inlet
        )
        desal.s05 = Arc(source=desal.S1.PXR, destination=desal.PXR.low_pressure_inlet)
        desal.s06 = Arc(
            source=desal.PXR.low_pressure_outlet, destination=desal.P2.inlet
        )
        desal.s07 = Arc(source=desal.P2.outlet, destination=desal.M1.P2)
        m.fs.s_disposal = Arc(
            source=desal.PXR.high_pressure_outlet, destination=m.fs.disposal.inlet
        )
    elif erd_type == "pump_as_turbine":
        m.fs.s_tb_desal = Arc(
            source=m.fs.tb_prtrt_desal.outlet, destination=desal.P1.inlet
        )
        desal.s01 = Arc(source=desal.P1.outlet, destination=desal.RO.inlet)
        desal.s02 = Arc(source=desal.RO.retentate, destination=desal.ERD.inlet)
        m.fs.s_disposal = Arc(source=desal.ERD.outlet, destination=m.fs.disposal.inlet)
    m.fs.s_desal_tb = Arc(
        source=desal.RO.permeate, destination=m.fs.tb_desal_psttrt.inlet
    )

    m.fs.s_tb_psttrt = Arc(
        source=m.fs.tb_desal_psttrt.outlet, destination=psttrt.storage_tank_2.inlet
    )
    psttrt.s01 = Arc(
        source=psttrt.storage_tank_2.outlet, destination=psttrt.uv_aop.inlet
    )
    psttrt.s02 = Arc(
        source=psttrt.uv_aop.treated, destination=psttrt.co2_addition.inlet
    )
    psttrt.s03 = Arc(
        source=psttrt.co2_addition.outlet, destination=psttrt.lime_addition.inlet
    )
    psttrt.s04 = Arc(
        source=psttrt.lime_addition.outlet, destination=psttrt.storage_tank_3.inlet
    )
    m.fs.s_municipal = Arc(
        source=psttrt.storage_tank_3.outlet, destination=m.fs.municipal.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.prop_desal.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.prop_desal.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "TDS")
    )
    # set unit model values
    iscale.set_scaling_factor(desal.P1.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.RO.area, 1e-4)
    if erd_type == "pressure_exchanger":
        iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.low_pressure_side.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.high_pressure_side.work, 1e-5)
    elif erd_type == "pump_as_turbine":
        iscale.set_scaling_factor(desal.ERD.control_volume.work, 1e-5)
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # ---specifications---
    # feed
    flow_vol = 0.3092 * pyunits.m**3 / pyunits.s
    conc_mass_tds = 35 * pyunits.kg / pyunits.m**3
    conc_mass_tss = 0.03 * pyunits.kg / pyunits.m**3
    temperature = 298 * pyunits.K
    pressure = 1e5 * pyunits.Pa

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    solve(m.fs.feed)

    m.fs.tb_prtrt_desal.properties_out[0].temperature.fix(temperature)
    m.fs.tb_prtrt_desal.properties_out[0].pressure.fix(pressure)

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
    # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    desal.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = 70e5 * pyunits.Pa
    desal.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    # RO unit
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    desal.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    desal.RO.width.fix(1000)  # stage width [m]
    desal.RO.area.fix(
        flow_vol * 4.5e4 * pyunits.s / pyunits.m
    )  # stage area [m2] TODO: replace with actual area

    if m.erd_type == "pressure_exchanger":
        # splitter (no degrees of freedom)

        # pressure exchanger, 1 degree of freedom (efficiency)
        desal.PXR.efficiency_pressure_exchanger.fix(0.95)

        # pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        desal.P2.efficiency_pump.fix(0.80)

        # mixer, no degrees of freedom
    elif m.erd_type == "pump_as_turbine":
        # ERD, 2 degrees of freedom (efficiency, outlet pressure)
        desal.ERD.efficiency_pump.fix(0.95)
        desal.ERD.control_volume.properties_out[0].pressure.fix(
            101325
        )  # atmospheric pressure [Pa]

    # ---posttreatment---
    # storage tank 2
    psttrt.storage_tank_2.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_2.storage_time.fix(1)

    # uv aop
    m.db.get_unit_operation_parameters("uv_aop")
    psttrt.uv_aop.load_parameters_from_database(use_default_removal=True)
    psttrt.uv_aop.uv_reduced_equivalent_dose.fix(
        350
    )  # TODO: check this was the right thing to fix
    psttrt.uv_aop.uv_transmittance_in.fix(
        0.95
    )  # TODO: check this was the right thing to fix

    # co2 addition
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


def initialize_system(m):
    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # initialize feed
    solve(m.fs.feed)

    # initialize pretreatment
    propagate_state(m.fs.s_feed)
    flags = fix_state_vars(prtrt.intake.properties)
    solve(prtrt)
    revert_state_vars(prtrt.intake.properties, flags)

    # initialize desalination
    propagate_state(m.fs.s_prtrt_tb)
    m.fs.tb_prtrt_desal.properties_out[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.tb_prtrt_desal.properties_in[0].flow_mass_comp["H2O"]
    )
    m.fs.tb_prtrt_desal.properties_out[0].flow_mass_phase_comp["Liq", "TDS"] = value(
        m.fs.tb_prtrt_desal.properties_in[0].flow_mass_comp["tds"]
    )

    desal.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_comp["H2O"]
    )
    desal.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"] = value(
        m.fs.feed.properties[0].flow_mass_comp["tds"]
    )
    desal.RO.feed_side.properties_in[0].temperature = value(
        m.fs.tb_prtrt_desal.properties_out[0].temperature
    )
    desal.RO.feed_side.properties_in[0].pressure = value(
        desal.P1.control_volume.properties_out[0].pressure
    )
    desal.RO.initialize()

    propagate_state(m.fs.s_tb_desal)
    if m.erd_type == "pressure_exchanger":
        flags = fix_state_vars(desal.S1.mixed_state)
        solve(desal)
        revert_state_vars(desal.S1.mixed_state, flags)
    elif m.erd_type == "pump_as_turbine":
        flags = fix_state_vars(desal.P1.control_volume.properties_in)
        solve(desal)
        revert_state_vars(desal.P1.control_volume.properties_in, flags)

    # initialize posttreatment
    propagate_state(m.fs.s_desal_tb)
    m.fs.tb_desal_psttrt.properties_out[0].flow_mass_comp["H2O"] = value(
        m.fs.tb_desal_psttrt.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.tb_desal_psttrt.properties_out[0].flow_mass_comp["tds"] = value(
        m.fs.tb_desal_psttrt.properties_in[0].flow_mass_phase_comp["Liq", "TDS"]
    )

    propagate_state(m.fs.s_tb_psttrt)
    flags = fix_state_vars(psttrt.storage_tank_2.properties)
    solve(psttrt)
    revert_state_vars(psttrt.storage_tank_2.properties, flags)


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display_results(m):
    m.fs.feed.report()
    m.fs.pretreatment.intake.report()
    m.fs.pretreatment.ferric_chloride_addition.report()
    m.fs.pretreatment.chlorination.report()
    m.fs.pretreatment.static_mixer.report()
    m.fs.pretreatment.storage_tank_1.report()
    m.fs.pretreatment.media_filtration.report()
    m.fs.pretreatment.backwash_handling.report()
    m.fs.pretreatment.anti_scalant_addition.report()
    m.fs.pretreatment.cartridge_filtration.report()
    if m.erd_type == "pressure_exchanger":
        m.fs.desalination.S1.report()
        m.fs.desalination.P1.report()
        m.fs.desalination.P2.report()
        m.fs.desalination.M1.report()
        m.fs.desalination.RO.report()
        m.fs.desalination.PXR.report()
    elif m.erd_type == "pump_as_turbine":
        m.fs.desalination.P1.report()
        m.fs.desalination.RO.report()
        m.fs.desalination.ERD.report()
    m.fs.posttreatment.storage_tank_2.report()
    m.fs.posttreatment.uv_aop.report()
    m.fs.posttreatment.co2_addition.report()
    m.fs.posttreatment.lime_addition.report()
    m.fs.posttreatment.storage_tank_3.report()
    m.fs.municipal.report()
    m.fs.landfill.report()
    m.fs.disposal.report()


def add_costing(m):
    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # Add costing package for zero-order units
    m.fs.zo_costing = ZeroOrderCosting()
    m.fs.ro_costing = WaterTAPCosting()

    # Add costing to zero order units
    # Pre-treatment units
    # This really looks like it should be a feed block in its own right
    # prtrt.intake.costing = UnitModelCostingBlock(default={
    #     "flowsheet_costing_block": m.fs.zo_costing})

    prtrt.ferric_chloride_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.chlorination.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.static_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.storage_tank_1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.media_filtration.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.backwash_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.anti_scalant_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    prtrt.cartridge_filtration.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )

    # RO Train
    # RO equipment is costed using more detailed costing package
    desal.P1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.ro_costing,
        costing_method_arguments={"cost_electricity_flow": False},
    )
    desal.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.ro_costing)
    if m.erd_type == "pressure_exchanger":
        # desal.S1.costing = UnitModelCostingBlock(default={
        #     "flowsheet_costing_block": m.fs.ro_costing})
        desal.M1.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.ro_costing
        )
        desal.PXR.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.ro_costing
        )
        desal.P2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.ro_costing,
            costing_method_arguments={"cost_electricity_flow": False},
        )
    elif m.erd_type == "pump_as_turbine":
        pass
        # desal.ERD.costing = UnitModelCostingBlock(default={
        #     "flowsheet_costing_block": m.fs.ro_costing})
    else:
        raise ConfigurationError(
            f"erd_type was {m.erd_type}, costing only implemented "
            "for pressure_exchanger or pump_as_turbine"
        )

    # For non-zero order unit operations, we need to register costed flows
    # separately.
    # However, to keep costs consistent, we will register these with the ZO
    # Costing package
    m.fs.zo_costing.cost_flow(desal.P1.work_mechanical[0], "electricity")
    if m.erd_type == "pressure_exchanger":
        m.fs.zo_costing.cost_flow(desal.P2.work_mechanical[0], "electricity")
    elif m.erd_type == "pump_as_turbine":
        pass
        # m.fs.zo_costing.cost_flow(
        #     desal.ERD.work_mechanical[0], "electricity")
    else:
        raise ConfigurationError(
            f"erd_type was {m.erd_type}, costing only implemented "
            "for pressure_exchanger or pump_as_turbine"
        )

    # Post-treatment units
    psttrt.storage_tank_2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    psttrt.uv_aop.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    psttrt.co2_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    psttrt.lime_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    psttrt.storage_tank_3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )

    # Product and disposal
    m.fs.municipal.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )
    m.fs.landfill.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.zo_costing
    )

    # Aggregate unit level costs and calculate overall process costs
    m.fs.zo_costing.cost_process()
    m.fs.ro_costing.cost_process()

    # Combine results from costing packages and calculate overall metrics
    @m.Expression()
    def total_capital_cost(b):
        return (
            pyunits.convert(
                m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2018
            )
            + m.fs.ro_costing.total_investment_cost
        )

    @m.Expression()
    def total_operating_cost(b):
        return (
            pyunits.convert(
                m.fs.zo_costing.total_fixed_operating_cost,
                to_units=pyunits.USD_2018 / pyunits.year,
            )
            + pyunits.convert(
                m.fs.zo_costing.total_variable_operating_cost,
                to_units=pyunits.USD_2018 / pyunits.year,
            )
            + m.fs.ro_costing.total_operating_cost
        )

    @m.Expression()
    def LCOW(b):
        return (
            b.total_capital_cost * b.fs.zo_costing.capital_recovery_factor
            + b.total_operating_cost
        ) / (
            pyunits.convert(
                b.fs.municipal.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.fs.zo_costing.utilization_factor
        )

    assert_units_consistent(m)


def initialize_costing(m):
    m.fs.zo_costing.initialize()
    m.fs.ro_costing.initialize()


def display_costing(m):
    # m.fs.zo_costing.display()
    # m.fs.ro_costing.display()

    m.total_capital_cost.display()
    m.total_operating_cost.display()
    m.LCOW.display()

    print("\nUnit Capital Costs\n")
    for u in m.fs.zo_costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018)),
        )
    for u in m.fs.ro_costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018)),
        )

    print("\nUtility Costs\n")
    for f in m.fs.zo_costing.flow_types:
        print(
            f,
            " :   ",
            value(
                pyunits.convert(
                    m.fs.zo_costing.aggregate_flow_costs[f],
                    to_units=pyunits.USD_2018 / pyunits.year,
                )
            ),
        )


def export_to_ui():
    from watertap.ui.fsapi import FlowsheetInterface

    def noop(*args, **kwargs):
        return

    return FlowsheetInterface(
        name="Seawater RO",
        description="Seawater RO desalination",
        do_export=noop,
        do_build=noop,
        do_solve=noop,
    )


if __name__ == "__main__":
    main(erd_type="pressure_exchanger")
    # main(erd_type='pump_as_turbine')
