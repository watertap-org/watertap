from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)

from idaes.core.util.initialization import propagate_state
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
)

from idaes.models.unit_models import Product, Feed
from pyomo.environ import TransformationFactory
from pyomo.network import Arc
from watertap.unit_models.pressure_changer import Pump

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock

from watertap.unit_models.pseudo_steady_state import DeadVolume0D
from watertap.unit_models.pseudo_steady_state.reverse_osmosis_1D_with_holdup import (
    ReverseOsmosis1DwithHoldUp,
)
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
)

from watertap.unit_models.pseudo_steady_state.flushing import Flushing
from watertap.unit_models.pseudo_steady_state.conduit import Conduit
from pyomo.util.calc_var_value import calculate_variable_from_constraint


def build_ro_systems(use_ro_with_hold_up=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.operation_mode = "filtration"
    m.fs.ro_model_with_hold_up = use_ro_with_hold_up
    m.fs.properties = NaClParameterBlock()

    # Low concentration feed to the system
    m.fs.raw_feed = Feed(property_package=m.fs.properties)

    # Dead volume is included in both operation_modes
    m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

    # Touch relevant parameters
    m.fs.raw_feed.properties[0].flow_vol_phase
    m.fs.raw_feed.properties[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase

    # Feed pump
    m.fs.P1 = Pump(property_package=m.fs.properties)
    # Recirculation pump
    m.fs.P2 = Pump(property_package=m.fs.properties)
    # Mixer to combine feed and recycle streams
    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets=2,
        momentum_mixing_type=MomentumMixingType.equality,
    )
    if use_ro_with_hold_up:
        m.fs.RO = ReverseOsmosis1DwithHoldUp(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            module_type="spiral_wound",
            has_full_reporting=True,
        )
    else:
        m.fs.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            module_type="spiral_wound",
            has_full_reporting=True,
        )

    m.fs.product = Product(property_package=m.fs.properties)

    # Add connections
    m.fs.raw_feed_to_P1 = Arc(source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet)

    m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
    m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)

    m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

    m.fs.RO_permeate_to_product = Arc(
        source=m.fs.RO.permeate, destination=m.fs.product.inlet
    )

    m.fs.RO_retentate_to_dead_volume = Arc(
        source=m.fs.RO.retentate, destination=m.fs.dead_volume.inlet
    )
    m.fs.dead_volume_to_P2 = Arc(
        source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Touch relevant parameters
    m.fs.M1.inlet_1_state[0].flow_vol_phase
    m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
    m.fs.M1.inlet_2_state[0].flow_vol_phase
    m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
    m.fs.M1.mixed_state[0].flow_vol_phase
    m.fs.M1.mixed_state[0].conc_mass_phase_comp
    return m


def intialize_ro_systems(m):
    propagate_state(m.fs.raw_feed_to_P1)

    # m.fs.P2.outlet.pressure[0].fix(m.fs.P1.outlet.pressure[0].value)
    m.fs.P1.initialize()

    propagate_state(m.fs.P1_to_M1)
    copy_inlet_state_for_mixer(m)

    m.fs.M1.initialize()

    propagate_state(m.fs.M1_to_RO)
    m.fs.RO.initialize()

    propagate_state(m.fs.RO_permeate_to_product)
    propagate_state(m.fs.RO_retentate_to_dead_volume)

    m.fs.dead_volume.initialize()

    propagate_state(m.fs.dead_volume_to_P2)
    m.fs.P2.initialize()
    m.fs.P2.outlet.pressure[0].unfix()

    propagate_state(m.fs.P2_to_M1)

    m.fs.product.properties[0].flow_vol_phase["Liq"]
    m.fs.product.initialize()


def copy_inlet_state_for_mixer(m):
    for idx, obj in m.fs.M1.inlet_2.flow_mass_phase_comp.items():
        obj.value = m.fs.M1.inlet_1.flow_mass_phase_comp[idx].value


def build_flushing_with_RO(use_ro_with_hold_up=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.operation_mode = "flushing_with_filtration"
    m.fs.ro_model_with_hold_up = use_ro_with_hold_up
    m.fs.properties = NaClParameterBlock()

    # Low concentration feed to the system
    m.fs.raw_feed = Feed(property_package=m.fs.properties)
    # Low concentration feed to the system
    m.fs.conduit_feed = Feed(property_package=m.fs.properties)

    # Dead volume is included in both operation_modes
    m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

    # Touch relevant parameters
    m.fs.raw_feed.properties[0].flow_vol_phase
    m.fs.raw_feed.properties[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase

    # Feed pump
    m.fs.P1 = Pump(property_package=m.fs.properties)
    # Recirculation pump
    m.fs.P2 = Pump(property_package=m.fs.properties)
    # Mixer to combine feed and recycle streams
    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets=2,
        momentum_mixing_type=MomentumMixingType.equality,
    )
    if use_ro_with_hold_up:
        m.fs.RO = ReverseOsmosis1DwithHoldUp(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            module_type="spiral_wound",
            has_full_reporting=True,
        )
    else:
        m.fs.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            module_type="spiral_wound",
            has_full_reporting=True,
        )

    m.fs.product = Product(property_package=m.fs.properties)

    # Add connections
    m.fs.raw_feed_to_P1 = Arc(source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet)
    m.fs.conduit_feed_to_P2 = Arc(
        source=m.fs.conduit_feed.outlet, destination=m.fs.P2.inlet
    )

    m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
    m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)

    m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)
    m.fs.M1.mixed_state[0].conc_mass_phase_comp[...]
    m.fs.RO_permeate_to_product = Arc(
        source=m.fs.RO.permeate, destination=m.fs.product.inlet
    )

    m.fs.RO_retentate_to_dead_volume = Arc(
        source=m.fs.RO.retentate, destination=m.fs.dead_volume.inlet
    )
    # m.fs.ro_product_constraint = Constraint(
    #     expr=m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
    #     == m.fs.product.properties[0].flow_vol_phase["Liq"]
    # )
    # immitates a recycle
    # m.fs.ro_equal_mass_constraint = Constraint(
    #     expr=sum(
    #         m.fs.conduit_feed.properties[0].flow_mass_phase_comp["Liq", comp]
    #         for comp in m.fs.properties.component_list
    #     )
    #     == sum(
    #         m.fs.dead_volume.outlet.flow_mass_phase_comp[0, "Liq", comp]
    #         for comp in m.fs.properties.component_list
    #     )
    # )
    m.fs.ro_equal_mass_constraint = Constraint(
        expr=m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"]
        == m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase["Liq"]
    )
    m.fs.conduit_pressure_constraint = Constraint(
        expr=m.fs.conduit_feed.properties[0].pressure == m.fs.RO.retentate.pressure[0]
    )

    iscale.constraint_scaling_transform(m.fs.conduit_pressure_constraint, 1e-5)
    # m.fs.ro_product_constraint.deactivate()
    iscale.constraint_scaling_transform(m.fs.ro_equal_mass_constraint, 1 / 0.001)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Touch relevant parameters
    m.fs.M1.inlet_1_state[0].flow_vol_phase
    m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
    m.fs.M1.inlet_2_state[0].flow_vol_phase
    m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
    m.fs.M1.mixed_state[0].flow_vol_phase
    m.fs.M1.mixed_state[0].conc_mass_phase_comp
    return m


def initialize_flushing_with_RO(m):
    propagate_state(m.fs.raw_feed_to_P1)

    propagate_state(m.fs.conduit_feed_to_P2)
    # m.fs.P2.outlet.pressure[0].fix(m.fs.P1.outlet.pressure[0].value)
    m.fs.P1.initialize()
    m.fs.P2.initialize()

    m.fs.P1.outlet.pressure[0].unfix()
    m.fs.P2.outlet.pressure[0].unfix()
    propagate_state(m.fs.P1_to_M1)

    propagate_state(m.fs.P2_to_M1)
    m.fs.M1.initialize()

    propagate_state(m.fs.M1_to_RO)
    m.fs.RO.initialize()
    m.fs.conduit_feed.properties[0].pressure = m.fs.RO.retentate.pressure[0].value
    m.fs.conduit_feed.properties[0].pressure.unfix()
    propagate_state(m.fs.RO_permeate_to_product)
    propagate_state(m.fs.RO_retentate_to_dead_volume)

    m.fs.dead_volume.initialize()
    m.fs.product.properties[0].flow_vol_phase["Liq"]

    m.fs.product.initialize()


def build_flushing_unit(
    operation_mode=None,
):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.operation_mode = "flushing"
    m.fs.properties = NaClParameterBlock()

    # Low concentration feed to the system
    m.fs.raw_feed = Feed(property_package=m.fs.properties)

    # Dead volume is included in both operation_modes
    m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

    # Touch relevant parameters
    m.fs.raw_feed.properties[0].flow_vol_phase
    m.fs.raw_feed.properties[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp
    m.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P2 = Pump(property_package=m.fs.properties)

    # Add connections
    m.fs.raw_feed_to_P1 = Arc(source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet)
    m.fs.P1_to_dead_volume = Arc(
        source=m.fs.P1.outlet, destination=m.fs.dead_volume.inlet
    )

    m.fs.dead_volume_to_P2 = Arc(
        source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
    )

    m.fs.flushing = Flushing()

    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def initialize_flushing_unit(m):
    propagate_state(m.fs.raw_feed_to_P1)
    m.fs.P1.initialize()

    propagate_state(m.fs.P1_to_dead_volume)
    m.fs.dead_volume.initialize()

    propagate_state(m.fs.dead_volume_to_P2)
    m.fs.P2.initialize()
    m.fs.P2.outlet.pressure[0].fix(101325 * pyunits.Pa)

    calculate_variable_from_constraint(
        m.fs.flushing.pre_flushing_concentration,
        m.fs.pre_flushing_conc_constraint,
    )
    calculate_variable_from_constraint(
        m.fs.flushing.mean_residence_time,
        m.fs.dead_volume_residence_time_constraint,
    )
    m.fs.flushing.mean_residence_time.fix()

    m.fs.flushing.flushing_time.fix()
    m.fs.flushing.flushing_time = value(m.fs.flushing.mean_residence_time)
    m.fs.flushing.pre_flushing_concentration.fix()
    m.fs.flushing.post_flushing_concentration.fix(
        value(m.fs.flushing.pre_flushing_concentration)
    )
    m.fs.flushing.post_flushing_concentration.unfix()
    m.fs.flushing.flushing_efficiency.unfix()
    m.fs.flushing.initialize()


def build_flushing_unit_only(mp, start_period=1, end_period=None):

    m_start = mp.get_active_process_blocks()[start_period]

    mp.flushing = Flushing()
    mp.flushing.start_period = start_period
    mp.flushing.end_period = end_period
    if end_period is None:
        blks = mp.get_active_process_blocks()[start_period:]
        m_end = mp.get_active_process_blocks()[-1]
    else:
        blks = mp.get_active_process_blocks()[start_period:end_period]
        m_end = mp.get_active_process_blocks()[end_period]

    # @mp.Constraint(list(range(start_period, len(blks))))
    # def eq_flushing_time_constraint(m, i):
    #     return blks[0].fs.RO.feed_side.accumulation_time[0] == (
    #         blks[i].fs.RO.feed_side.accumulation_time[0]
    #     )

    # for blk in mp.get_active_process_blocks():
    #     print(blk.name, blk.fs.operation_mode)
    for blk in blks:
        # print(blk.name, blk.fs.operation_mode)
        blk.fs.RO.feed_side.accumulation_time.unfix()
        blk.fs.RO.feed_side.accumulation_time[0].setlb(1e-8)

    # Constraints
    @mp.Constraint()
    def dead_volume_residence_time_constraint(m):
        return mp.flushing.mean_residence_time == (
            pyunits.convert(
                (
                    m_start.fs.RO.feed_side.volume
                    + m_start.fs.dead_volume.volume[0, "Liq"]
                )
                / (
                    m_start.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
                    + m_start.fs.conduit_feed.properties[0].flow_vol_phase["Liq"]
                ),
                to_units=pyunits.s,
            )
        )

    # Calculate pre-flushing/ dead volume delta state concentration. Concentration before
    # flushing should be the delta state concentration
    @mp.Constraint()
    def pre_flushing_conc_constraint(m):
        return (
            mp.flushing.pre_flushing_concentration
            == m_start.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                0, 1, "Liq", "NaCl"
            ]
            * m_start.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, 1, "Liq"]
        )

    # Concentration after flushing should be the dead volume properties out concentration
    @mp.Constraint()
    def post_flushing_conc_constraint(m):
        return (
            m_end.fs.RO.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
            == mp.flushing.post_flushing_concentration
        )

    @mp.Constraint()
    def flushing_concentration_constraint(m):
        return (
            m_start.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"]
            == mp.flushing.flushing_feed_concentration
        )

    iscale.constraint_scaling_transform(mp.pre_flushing_conc_constraint, 1)
    # iscale.constraint_scaling_transform(mp.post_flushing_conc_constraint, 1)
    iscale.set_scaling_factor(mp.flushing.pre_flushing_concentration, 1)
    iscale.set_scaling_factor(mp.flushing.post_flushing_concentration, 1)

    return mp


def initialize_flushing(mp):
    # if mp.flushing.end_period is None:
    #     blks = mp.get_active_process_blocks()[mp.flushing.start_period :]
    #     m_end = mp.get_active_process_blocks()[-1]
    # else:
    #     blks = mp.get_active_process_blocks()[
    #         mp.flushing.start_period : mp.flushing.end_period
    #     ]
    #     m_end = mp.get_active_process_blocks()[mp.flushing.end_period]
    calculate_variable_from_constraint(
        mp.flushing.pre_flushing_concentration,
        mp.pre_flushing_conc_constraint,
    )
    calculate_variable_from_constraint(
        mp.flushing.post_flushing_concentration,
        mp.post_flushing_conc_constraint,
    )

    calculate_variable_from_constraint(
        mp.flushing.mean_residence_time,
        mp.dead_volume_residence_time_constraint,
    )

    calculate_variable_from_constraint(
        mp.flushing.flushing_feed_concentration,
        mp.flushing_concentration_constraint,
    )
    mp.flushing.mean_residence_time.fix()

    mp.flushing.flushing_time.unfix()
    mp.flushing.flushing_time = value(mp.flushing.mean_residence_time)
    mp.flushing.pre_flushing_concentration.fix()
    mp.flushing.post_flushing_concentration.fix()
    mp.flushing.flushing_efficiency = 0.9
    mp.flushing.flushing_efficiency.unfix()
    mp.flushing.pre_flushing_concentration.unfix()
    mp.flushing.flushing_time.unfix()
    mp.flushing.flushing_time.setub(1000)
    mp.flushing.flushing_time.setlb(1e-8)
    mp.flushing.mean_residence_time.unfix()

    mp.flushing.pre_flushing_concentration.unfix()
    mp.flushing.post_flushing_concentration.unfix()
    mp.flushing.flushing_feed_concentration.unfix()


def build_conduit(mp):
    mp.conduit = Conduit()

    return mp
