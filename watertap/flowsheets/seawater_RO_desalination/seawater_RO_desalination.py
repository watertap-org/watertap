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
import os
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
    Block,
    Var,
    Constraint,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
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
from watertap.core.solvers import get_solver

# Set up logger
_log = idaeslog.getLogger(__name__)


def main(erd_type="pressure_exchanger", RO_1D=False, **kwargs):

    flow_vol = kwargs.get("flow_vol", 0.3092)
    conc_mass_tds = kwargs.get("conc_mass_tds", 35)
    conc_mass_tss = kwargs.get("conc_mass_tss", 0.03)

    m = build(
        erd_type=erd_type,
        RO_1D=RO_1D,
        flow_vol=flow_vol,
        conc_mass_tds=conc_mass_tds,
        conc_mass_tss=conc_mass_tss,
    )

    scale_model(m)
    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    initialize_system(m)

    assert_degrees_of_freedom(m, 0)

    solve(
        m,
        checkpoint=f" initializing {erd_type} system",
        tee=False,
    )
    # display_results(m)

    add_costing(m)
    m.fs.costing.initialize()
    assert_degrees_of_freedom(m, 0)

    solve(m, tee=True, checkpoint=f" solve {erd_type} flowsheet with costing")
    display_costing(m)

    return m


def build(
    erd_type=None, RO_1D=False, flow_vol=0.3092, conc_mass_tds=35, conc_mass_tss=0.03
):

    # Define flow conditions and attach to model

    rho = 1023.027 * pyunits.kg / pyunits.m**3
    m = ConcreteModel()
    m.flow_vol = flow_vol * pyunits.m**3 / pyunits.s

    m.flow_mass_water = pyunits.convert(
        m.flow_vol * rho, to_units=pyunits.kg / pyunits.s
    )
    m.conc_mass_tds = conc_mass_tds * pyunits.g / pyunits.liter
    m.flow_mass_tds = pyunits.convert(
        m.flow_vol * m.conc_mass_tds, to_units=pyunits.kg / pyunits.s
    )
    m.conc_mass_tss = conc_mass_tss * pyunits.g / pyunits.liter
    m.flow_mass_tss = pyunits.convert(
        m.flow_vol * m.conc_mass_tss, to_units=pyunits.kg / pyunits.s
    )
    m.temperature = 298.15 * pyunits.K
    m.pressure = 1 * pyunits.bar

    m.db = Database()
    m.erd_type = erd_type

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.rho = Var(initialize=rho, units=pyunits.kg / pyunits.m**3)
    m.fs.rho.fix()

    m.fs.flow_vol = Var(initialize=flow_vol, units=pyunits.m**3 / pyunits.s)
    m.fs.flow_vol.fix()

    m.fs.conc_mass_tds = Var(initialize=conc_mass_tds, units=pyunits.g / pyunits.L)
    m.fs.conc_mass_tds.fix()

    m.fs.conc_mass_tss = Var(initialize=conc_mass_tss, units=pyunits.g / pyunits.L)
    m.fs.conc_mass_tss.fix()

    m.fs.properties = MCASParameterBlock(
        solute_list=["TDS", "TSS"],
        diffusivity_data={("Liq", "TDS"): 1.47e-9, ("Liq", "TSS"): 1e-9},
        mw_data={"TDS": 31.4e-3, "TSS": 100e-3},
        material_flow_basis=MaterialFlowBasis.mass,
        ignore_neutral_charge=True,
        density_calculation=DensityCalculation.seawater,
    )

    # Block structure
    prtrt = m.fs.pretreatment = Block()
    desal = m.fs.desalination = Block()
    psttrt = m.fs.posttreatment = Block()

    # Unit models
    m.fs.feed = Feed(property_package=m.fs.properties)
    # Touch properties to calculate_state and scale properly
    m.fs.feed.properties[0].flow_vol_phase
    m.fs.feed.properties[0].conc_mass_phase_comp

    # Pretreatment
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

    # Desalination
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

    # Posttreatment
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

    # Product and disposal
    m.fs.municipal = MunicipalDrinkingZO(
        property_package=m.fs.properties, database=m.db
    )
    m.fs.landfill = LandfillZO(property_package=m.fs.properties, database=m.db)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Connections
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

    return m


def scale_model(m):

    desal = m.fs.desalination

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.flow_mass_water),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.flow_mass_tds),
        index=("Liq", "TDS"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.flow_mass_tss),
        index=("Liq", "TSS"),
    )

    iscale.set_scaling_factor(desal.P1.control_volume.work, 1e-5)
    iscale.set_scaling_factor(desal.RO.area, 1e-4)

    if m.erd_type == "pressure_exchanger":
        iscale.set_scaling_factor(desal.P2.control_volume.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.feed_side.work, 1e-5)
        iscale.set_scaling_factor(desal.PXR.brine_side.work, 1e-5)

    elif m.erd_type == "pump_as_turbine":
        iscale.set_scaling_factor(desal.ERD.control_volume.work, 1e-5)

    if m.erd_type == "pressure_exchanger":
        desal.S1.mixed_state[0].flow_vol_phase
        desal.RO.feed_side.properties[0, 1].flow_vol_phase

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)


def set_operating_conditions(m):
    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    m.fs.feed.properties.calculate_state(
        var_args={
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.conc_mass_tds,
            ("conc_mass_phase_comp", ("Liq", "TSS")): m.conc_mass_tss,
            ("flow_vol_phase", "Liq"): m.flow_vol,
            ("temperature", None): m.temperature,
            ("pressure", None): m.pressure,
        },
        hold_state=False,
    )

    m.fs.feed.properties[0].pressure.fix()
    m.fs.feed.properties[0].temperature.fix()
    # re-fix density after calculating state
    rho = value(m.fs.feed.properties[0].dens_mass_phase["Liq"])
    m.fs.rho.fix(rho)

    m.fs.feed.flow_mass_water_constr = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        == pyunits.convert(m.fs.flow_vol * m.fs.rho, to_units=pyunits.kg / pyunits.s)
    )
    m.fs.feed.flow_mass_tds_constr = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        == pyunits.convert(
            m.fs.flow_vol * m.fs.conc_mass_tds, to_units=pyunits.kg / pyunits.s
        )
    )
    m.fs.feed.flow_mass_tss_constr = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TSS"]
        == pyunits.convert(
            m.fs.flow_vol * m.fs.conc_mass_tss, to_units=pyunits.kg / pyunits.s
        )
    )

    # ---Pretreatment---
    # Intake
    m.db.get_unit_operation_parameters("sw_onshore_intake")
    prtrt.intake.load_parameters_from_database()
    # Ferric chloride
    m.db.get_unit_operation_parameters("chemical_addition")
    prtrt.ferric_chloride_addition.load_parameters_from_database()
    prtrt.ferric_chloride_addition.chemical_dosage.fix(20)

    # Chlorination
    m.db.get_unit_operation_parameters("chlorination")
    prtrt.chlorination.load_parameters_from_database(use_default_removal=True)

    # Static mixer
    m.db.get_unit_operation_parameters("static_mixer")
    prtrt.static_mixer.load_parameters_from_database(use_default_removal=True)

    # Storage tank
    m.db.get_unit_operation_parameters("storage_tank")
    prtrt.storage_tank_1.load_parameters_from_database(use_default_removal=True)
    prtrt.storage_tank_1.storage_time.fix(2)

    # Media filtration
    m.db.get_unit_operation_parameters("media_filtration")
    prtrt.media_filtration.load_parameters_from_database(use_default_removal=True)

    # Backwash handling
    m.db.get_unit_operation_parameters("backwash_solids_handling")
    prtrt.backwash_handling.load_parameters_from_database(use_default_removal=True)

    # Anti-scalant
    prtrt.anti_scalant_addition.load_parameters_from_database()
    prtrt.anti_scalant_addition.chemical_dosage.fix(5)

    # Cartridge filtration
    m.db.get_unit_operation_parameters("cartridge_filtration")
    prtrt.cartridge_filtration.load_parameters_from_database(use_default_removal=True)

    # ---Desalination---
    # Pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    desal.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = 70e5 * pyunits.Pa
    desal.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    # RO Unit
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    desal.RO.B_comp[0, "TSS"].fix(1e-10)  # membrane salt permeability coefficient [m/s]

    desal.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.feed_side.spacer_porosity.fix(0.9)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]

    # NOTE: 3234.153 parameter comes from original SW desal flowsheet
    # which had a flow rate of 0.3092 and resulted in stage width of 1000 m
    # 1000 / 0.3092 = 3234.153
    width_guess = value(m.flow_vol) * 3234.153
    if width_guess > desal.RO.width.ub:
        desal.RO.width.setub(value(width_guess) * 2)
    desal.RO.width.fix(width_guess)  # stage width [m]

    # NOTE: 45000 parameter comes from original SW desal flowsheet
    # which had a flow rate of 0.3092 and resulted in total area of 13914 m2
    # 13914 / 0.3092 = 45000
    area_guess = value(m.flow_vol) * 45000  # rough estimate for stage area [m2]
    if area_guess > desal.RO.area.ub:
        desal.RO.area.setub(value(area_guess) * 2)
    # stage area [m2] TODO: replace with actual area
    desal.RO.area.fix(area_guess)

    desal.RO.recovery_mass_phase_comp.setlb(0)
    desal.RO.flux_mass_phase_comp.setlb(0)

    if m.erd_type == "pressure_exchanger":
        # Splitter (no degrees of freedom)

        # Pressure exchanger, 1 degree of freedom (efficiency)
        desal.PXR.efficiency_pressure_exchanger.fix(0.95)

        # Pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        desal.P2.efficiency_pump.fix(0.80)

        # Mixer, no degrees of freedom
    elif m.erd_type == "pump_as_turbine":
        # ERD, 2 degrees of freedom (efficiency, outlet pressure)
        desal.ERD.efficiency_pump.fix(0.95)
        # Atmospheric pressure [Pa]
        desal.ERD.control_volume.properties_out[0].pressure.fix(101325)

    # ---Posttreatment---
    # Storage tank 2
    psttrt.storage_tank_2.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_2.storage_time.fix(1)

    # UV/AOP
    m.db.get_unit_operation_parameters("uv_aop")
    psttrt.uv_aop.load_parameters_from_database(use_default_removal=True)
    psttrt.uv_aop.uv_reduced_equivalent_dose.fix(350)
    psttrt.uv_aop.uv_transmittance_in.fix(0.95)

    # CO2 addition
    m.db.get_unit_operation_parameters("co2_addition")
    psttrt.co2_addition.load_parameters_from_database(use_default_removal=True)

    # Lime
    psttrt.lime_addition.load_parameters_from_database()
    psttrt.lime_addition.chemical_dosage.fix(2.3)

    # Storage tank 3
    psttrt.storage_tank_3.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_3.storage_time.fix(1)

    # ---Product and disposal---
    m.db.get_unit_operation_parameters("municipal_drinking")
    m.fs.municipal.load_parameters_from_database()

    m.db.get_unit_operation_parameters("landfill")
    m.fs.landfill.load_parameters_from_database()


def initialize_pretreatment(m):
    prtrt = m.fs.pretreatment

    prtrt.intake.initialize()
    propagate_state(prtrt.s01)
    prtrt.ferric_chloride_addition.initialize()
    propagate_state(prtrt.s02)
    prtrt.chlorination.initialize()
    propagate_state(prtrt.s03)
    prtrt.static_mixer.initialize()
    propagate_state(prtrt.s04)
    prtrt.storage_tank_1.initialize()
    propagate_state(prtrt.s05)
    prtrt.media_filtration.initialize()
    propagate_state(prtrt.s06)
    prtrt.backwash_handling.initialize()
    propagate_state(prtrt.s07)
    prtrt.anti_scalant_addition.initialize()
    propagate_state(prtrt.s08)
    prtrt.cartridge_filtration.initialize()


def initialize_posttreatment(m):
    psttrt = m.fs.posttreatment

    psttrt.storage_tank_2.initialize()
    propagate_state(psttrt.s01)
    psttrt.uv_aop.initialize()
    propagate_state(psttrt.s02)
    psttrt.co2_addition.initialize()
    propagate_state(psttrt.s03)
    psttrt.lime_addition.initialize()
    propagate_state(psttrt.s04)
    psttrt.storage_tank_3.initialize()


def initialize_system(m):

    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # Initialize feed
    solve(m.fs.feed, checkpoint=" initializing feed")

    # Initialize pretreatment
    propagate_state(m.fs.s_feed)

    initialize_pretreatment(m)

    # Initialize desalination
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
            checkpoint=f" initializing desalination with {m.erd_type}",
        )
        revert_state_vars(desal.S1.mixed_state, flags)
    else:
        desal.P1.initialize()
        propagate_state(desal.s01)
        desal.RO.initialize(outlvl=idaeslog.DEBUG)
        propagate_state(desal.s02)
        desal.ERD.initialize()
        propagate_state(m.fs.s_disposal)

    propagate_state(desal.s_permeate_to_storage)

    # Initialize posttreatment
    initialize_posttreatment(m)

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

    # Zero order costing
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "seawater_RO_desalination_global_costing.yaml",
    )

    # Add costing package
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
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

    header = "{:<35} | {:<20} | {:<25}"
    row_fmt = "{:<35} | {:<20,.2f} | {:<25}"
    divider = "-" * 87

    # System costing
    print(divider)
    print(header.format("Overall Metric", "Value", "Units"))
    print(divider)
    print(
        row_fmt.format(
            "Total Capital Cost",
            value(m.fs.costing.total_capital_cost),
            str(pyunits.get_units(m.fs.costing.total_operating_cost)),
        )
    )
    print(
        row_fmt.format(
            "Total Operating Cost",
            value(m.fs.costing.total_operating_cost),
            str(pyunits.get_units(m.fs.costing.total_operating_cost)),
        )
    )
    print(
        row_fmt.format(
            "LCOW", value(m.fs.costing.LCOW), str(pyunits.get_units(m.fs.costing.LCOW))
        )
    )
    print(
        row_fmt.format(
            "Specific Energy Consumption",
            value(m.fs.costing.specific_energy_consumption),
            str(pyunits.get_units(m.fs.costing.specific_energy_consumption)),
        )
    )
    print(divider + "\n")

    # Capital cost by unit process
    print(divider)
    print(header.format("Unit Process", f"CAPEX", "Units"))
    print(divider)
    for u in m.fs.costing._registered_unit_costing:
        cost_val = value(u.capital_cost)
        p_name = u.name.split(".")[-2]
        u_str = str(pyunits.get_units(u.capital_cost))
        print(row_fmt.format(p_name, cost_val, u_str))
    print(divider + "\n")

    # Operating costs by material/energy flow
    print(divider)
    print(header.format("Material/Energy Flow", "Value", "Units"))
    print(divider)
    for used_flow in m.fs.costing.used_flows:
        flow_val = value(m.fs.costing.aggregate_flow_costs[used_flow])
        print(
            row_fmt.format(
                used_flow,
                flow_val,
                str(pyunits.get_units(m.fs.costing.aggregate_flow_costs[used_flow])),
            )
        )
    print(divider)


if __name__ == "__main__":
    m = main(erd_type="pressure_exchanger", RO_1D=True)
