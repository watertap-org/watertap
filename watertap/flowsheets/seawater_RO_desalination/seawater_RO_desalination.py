#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
    Block,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.initialization import (
    propagate_state,
    fix_state_vars,
    revert_state_vars,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import UnitModelCostingBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    DensityCalculation,
    MaterialFlowBasis,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
)
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom, check_solve

from watertap.core.wt_database import Database
from watertap.unit_models.zero_order import (
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
from watertap.costing.zero_order_costing import ZeroOrderCosting

# Set up logger
_log = idaeslog.getLogger(__name__)


def main(erd_type="pressure_exchanger", RO_1D=False):
    m = build(erd_type=erd_type, RO_1D=RO_1D)

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    initialize_system(m)

    assert_degrees_of_freedom(m, 0)

    solve(
        m, checkpoint=f" solve flowsheet after initializing {erd_type} system", tee=True
    )
    display_results(m)

    add_costing(m)
    m.fs.costing.initialize()
    assert_degrees_of_freedom(m, 0)

    solve(m, tee=True, checkpoint=f" solve {erd_type} flowsheet with costing")
    display_costing(m)

    return m


def build(erd_type=None, RO_1D=False):
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()
    m.erd_type = erd_type

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["tds", "tss"],
        diffusivity_data={("Liq", "tds"): 1.47e-9, ("Liq", "tss"): 1e-9},
        mw_data={"tds": 31.4e-3, "tss": 100e-3},
        material_flow_basis=MaterialFlowBasis.mass,
        ignore_neutral_charge=True,
        density_calculation=DensityCalculation.seawater,
    )

    # block structure
    prtrt = m.fs.pretreatment = Block()
    desal = m.fs.desalination = Block()
    psttrt = m.fs.posttreatment = Block()

    # unit models
    m.fs.feed = Feed(property_package=m.fs.properties)
    # pretreatment
    prtrt.intake = SWOnshoreIntakeZO(property_package=m.fs.properties, database=m.db)
    prtrt.ferric_chloride_addition = ChemicalAdditionZO(
        property_package=m.fs.properties,
        database=m.db,
        process_subtype="ferric_chloride",
    )
    prtrt.chlorination = ChlorinationZO(property_package=m.fs.properties, database=m.db)
    prtrt.static_mixer = StaticMixerZO(property_package=m.fs.properties, database=m.db)
    prtrt.storage_tank_1 = StorageTankZO(
        property_package=m.fs.properties, database=m.db
    )
    prtrt.media_filtration = MediaFiltrationZO(
        property_package=m.fs.properties, database=m.db
    )
    prtrt.backwash_handling = BackwashSolidsHandlingZO(
        property_package=m.fs.properties, database=m.db
    )
    prtrt.anti_scalant_addition = ChemicalAdditionZO(
        property_package=m.fs.properties, database=m.db, process_subtype="anti-scalant"
    )
    prtrt.cartridge_filtration = CartridgeFiltrationZO(
        property_package=m.fs.properties, database=m.db
    )

    # desalination
    desal.P1 = Pump(property_package=m.fs.properties)
    if RO_1D:
        desal.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
        )
    else:
        desal.RO = ReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
        )
    desal.RO.width.setub(5000)
    desal.RO.area.setub(20000)
    if erd_type == "pressure_exchanger":
        desal.S1 = Separator(
            property_package=m.fs.properties, outlet_list=["P1", "PXR"]
        )
        desal.M1 = Mixer(
            property_package=m.fs.properties,
            momentum_mixing_type=MomentumMixingType.equality,
            inlet_list=["P1", "P2"],
        )
        desal.PXR = PressureExchanger(property_package=m.fs.properties)
        desal.P2 = Pump(property_package=m.fs.properties)
    elif erd_type == "pump_as_turbine":
        desal.ERD = EnergyRecoveryDevice(property_package=m.fs.properties)
    else:
        raise ConfigurationError(
            "erd_type was {}, but can only "
            "be pressure_exchanger or pump_as_turbine"
            "".format(erd_type)
        )

    # posttreatment
    psttrt.storage_tank_2 = StorageTankZO(
        property_package=m.fs.properties, database=m.db
    )
    psttrt.uv_aop = UVAOPZO(
        property_package=m.fs.properties,
        database=m.db,
        process_subtype="hydrogen_peroxide",
    )
    psttrt.co2_addition = CO2AdditionZO(property_package=m.fs.properties, database=m.db)
    psttrt.lime_addition = ChemicalAdditionZO(
        property_package=m.fs.properties, database=m.db, process_subtype="lime"
    )
    psttrt.storage_tank_3 = StorageTankZO(
        property_package=m.fs.properties, database=m.db
    )

    # product and disposal
    m.fs.municipal = MunicipalDrinkingZO(
        property_package=m.fs.properties, database=m.db
    )
    m.fs.landfill = LandfillZO(property_package=m.fs.properties, database=m.db)
    m.fs.disposal = Product(property_package=m.fs.properties)

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

    m.fs.s_landfill = Arc(
        source=prtrt.backwash_handling.byproduct, destination=m.fs.landfill.inlet
    )

    if erd_type == "pressure_exchanger":
        prtrt.s_09 = Arc(
            source=prtrt.cartridge_filtration.treated, destination=desal.S1.inlet
        )
        desal.s01 = Arc(source=desal.S1.P1, destination=desal.P1.inlet)
        desal.s02 = Arc(source=desal.P1.outlet, destination=desal.M1.P1)
        desal.s03 = Arc(source=desal.M1.outlet, destination=desal.RO.inlet)
        desal.s04 = Arc(source=desal.RO.retentate, destination=desal.PXR.brine_inlet)
        desal.s05 = Arc(source=desal.S1.PXR, destination=desal.PXR.feed_inlet)
        desal.s06 = Arc(source=desal.PXR.feed_outlet, destination=desal.P2.inlet)
        desal.s07 = Arc(source=desal.P2.outlet, destination=desal.M1.P2)
        m.fs.s_disposal = Arc(
            source=desal.PXR.brine_outlet, destination=m.fs.disposal.inlet
        )
    elif erd_type == "pump_as_turbine":
        prtrt.s_09 = Arc(
            source=prtrt.cartridge_filtration.treated, destination=desal.P1.inlet
        )
        desal.s01 = Arc(source=desal.P1.outlet, destination=desal.RO.inlet)
        desal.s02 = Arc(source=desal.RO.retentate, destination=desal.ERD.inlet)
        m.fs.s_disposal = Arc(source=desal.ERD.outlet, destination=m.fs.disposal.inlet)
    desal.s_permeate_to_storage = Arc(
        source=desal.RO.permeate, destination=psttrt.storage_tank_2.inlet
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
    # set unit model values
    iscale.set_scaling_factor(desal.P1.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.RO.area, 1e-4)
    if erd_type == "pressure_exchanger":
        iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.feed_side.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.brine_side.work, 1e-5)
    elif erd_type == "pump_as_turbine":
        iscale.set_scaling_factor(desal.ERD.control_volume.work, 1e-5)

    if erd_type == "pressure_exchanger":
        desal.S1.mixed_state[0].flow_vol_phase
        desal.RO.feed_side.properties[0, 1].flow_vol_phase
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

    m.fs.feed.temperature[0].fix(temperature)
    m.fs.feed.pressure[0].fix(pressure)
    iscale.set_scaling_factor(
        m.fs.feed.properties[0].flow_vol_phase["Liq"], value(10 / flow_vol)
    )
    iscale.set_scaling_factor(
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "tds"],
        value(10 / conc_mass_tds),
    )
    iscale.set_scaling_factor(
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "tss"],
        value(10 / conc_mass_tss),
    )

    m.fs.feed.properties.calculate_state(
        var_args={
            ("conc_mass_phase_comp", ("Liq", "tds")): conc_mass_tds,
            ("conc_mass_phase_comp", ("Liq", "tss")): conc_mass_tss,
            ("flow_vol_phase", "Liq"): flow_vol,
        },
        hold_state=True,
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "tds"]),
        index=("Liq", "tds"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "tss"]),
        index=("Liq", "tss"),
    )

    iscale.calculate_scaling_factors(m)

    # ---pretreatment---
    # intake
    m.db.get_unit_operation_parameters("sw_onshore_intake")
    prtrt.intake.load_parameters_from_database()
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
    for u in (prtrt.ferric_chloride_addition, prtrt.anti_scalant_addition):
        iscale.set_scaling_factor(u.chemical_flow_vol, 1e6)
        iscale.constraint_scaling_transform(u.chemical_flow_vol_constraint, 1e6)

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
    desal.RO.B_comp[0, "tss"].fix(1e-10)  # membrane salt permeability coefficient [m/s]

    desal.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.feed_side.spacer_porosity.fix(0.9)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    desal.RO.width.fix(1000)  # stage width [m]
    desal.RO.area.fix(
        flow_vol * 4.5e4 * pyunits.s / pyunits.m
    )  # stage area [m2] TODO: replace with actual area
    m.fs.desalination.RO.recovery_mass_phase_comp.setlb(None)
    m.fs.desalination.RO.flux_mass_phase_comp.setlb(None)
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
    solve(m.fs.feed, checkpoint="solve flowsheet after initializing feed")

    # initialize pretreatment
    propagate_state(m.fs.s_feed)
    flags = fix_state_vars(prtrt.intake.properties)
    solve(prtrt, checkpoint="solve flowsheet after initializing pre-treatment")
    revert_state_vars(prtrt.intake.properties, flags)

    # initialize desalination
    propagate_state(prtrt.s_09)

    if m.erd_type == "pressure_exchanger":
        comps = m.fs.properties.component_list
        desal.S1.split_fraction[0, "PXR"].fix(0.55)
        desal.S1.initialize()
        desal.S1.split_fraction[0, "PXR"].unfix()
        propagate_state(desal.s01)
        desal.P1.initialize()
        propagate_state(desal.s02)
        p1_out = desal.P1.control_volume.properties_out[0]
        for j in comps:
            desal.M1.P2_state[0].flow_mass_phase_comp["Liq", j].set_value(
                value(desal.S1.PXR_state[0].flow_mass_phase_comp["Liq", j])
            )
        desal.M1.P2_state[0].temperature.set_value(value(p1_out.temperature))
        desal.M1.P2_state[0].pressure.set_value(value(p1_out.pressure))
        desal.M1.initialize()
        propagate_state(desal.s03)
        desal.RO.initialize()
        split_vol = value(
            desal.RO.feed_side.properties[0, 1].flow_vol_phase["Liq"]
            / desal.S1.mixed_state[0].flow_vol_phase["Liq"]
        )
        desal.S1.split_fraction[0, "PXR"].fix(split_vol)
        desal.S1.initialize()
        desal.S1.split_fraction[0, "PXR"].unfix()
        propagate_state(desal.s01)
        desal.P1.initialize()
        propagate_state(desal.s02)
        propagate_state(desal.s05)
        propagate_state(desal.s04)
        desal.PXR.initialize()
        propagate_state(desal.s06)

        desal.P2.control_volume.properties_out[0].pressure.fix(value(p1_out.pressure))
        desal.P2.initialize()
        desal.P2.control_volume.properties_out[0].pressure.unfix()
        propagate_state(desal.s07)
        desal.M1.initialize()
        propagate_state(desal.s03)

        flags = fix_state_vars(desal.S1.mixed_state)
        solve(
            desal,
            checkpoint=f"solve flowsheet after initializing desalination with {m.erd_type}",
        )
        revert_state_vars(desal.S1.mixed_state, flags)
    else:
        desal.P1.initialize()
        propagate_state(desal.s01)
        desal.RO.initialize(outlvl=idaeslog.DEBUG)
        propagate_state(desal.s02)
        desal.ERD.initialize()
        propagate_state(m.fs.s_disposal)

    # initialize posttreatment
    propagate_state(desal.s_permeate_to_storage)

    flags = fix_state_vars(psttrt.storage_tank_2.properties)
    solve(psttrt, checkpoint="solve flowsheet after initializing post-treatment")
    revert_state_vars(psttrt.storage_tank_2.properties, flags)

    propagate_state(m.fs.s_municipal)
    m.fs.municipal.initialize()
    propagate_state(m.fs.s_disposal)
    m.fs.disposal.initialize()
    propagate_state(m.fs.s_landfill)
    m.fs.landfill.initialize()
    propagate_state(desal.s_permeate_to_storage)


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
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

    # Add costing package
    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2023
    # Add costing to zero order units
    # Pre-treatment units
    # Intake unit really looks like it should be a feed block in its own right
    prtrt.intake.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    prtrt.ferric_chloride_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.chlorination.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.static_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.storage_tank_1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.media_filtration.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.backwash_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.anti_scalant_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    prtrt.cartridge_filtration.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

    # RO Train
    # RO equipment is costed using more detailed costing package
    desal.P1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_electricity_flow": True},
    )
    desal.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    if m.erd_type == "pressure_exchanger":
        # NOTE: Costing for the S1 splitter is neglected. Keeping the commented line below for awareness.
        # desal.S1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        desal.M1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        desal.PXR.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        desal.P2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_electricity_flow": True},
        )
    elif m.erd_type == "pump_as_turbine":
        desal.ERD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    else:
        raise ConfigurationError(
            f"erd_type was {m.erd_type}, costing only implemented "
            "for pressure_exchanger or pump_as_turbine"
        )

    # Post-treatment units
    psttrt.storage_tank_2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    psttrt.uv_aop.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    psttrt.co2_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    psttrt.lime_addition.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )
    psttrt.storage_tank_3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

    # Product and disposal
    m.fs.municipal.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.landfill.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Aggregate unit level costs and calculate overall process costs
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.municipal.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.municipal.properties[0].flow_vol)
    assert_units_consistent(m)


def display_costing(m):
    m.fs.costing.total_capital_cost.display()
    m.fs.costing.total_operating_cost.display()
    m.fs.costing.LCOW.display()
    m.fs.costing.specific_energy_consumption.display()

    print("\nUnit Capital Costs\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=m.fs.costing.base_currency)),
        )

    print("\nUtility Costs\n")
    for f in m.fs.costing.used_flows:
        print(
            f,
            " :   ",
            value(
                pyunits.convert(
                    m.fs.costing.aggregate_flow_costs[f],
                    to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
                )
            ),
        )


if __name__ == "__main__":
    m = main(erd_type="pressure_exchanger", RO_1D=True)
