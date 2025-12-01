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
from pyomo.environ import TransformationFactory
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

from watertap.unit_models.pressure_changer import Pump
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
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
)
from watertap.unit_models.pseudo_steady_state import DeadVolume0D
from watertap.unit_models.pseudo_steady_state.flushing import Flushing

from watertap.costing import (
    WaterTAPCosting,
)

from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
)

import watertap.flowsheets.ccro.utils.utils as cc_utils
from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration

solver = get_solver()


def create_ccro_multiperiod(
    n_time_points=10, include_costing=True, cc_configuration=None
):
    """
    Create multiperiod model for CCRO system
    """
    if isinstance(cc_configuration, CCROConfiguration) is False:
        user_config = cc_configuration
        cc_configuration = CCROConfiguration()
        cc_configuration.update(user_config)

    watertap_solver = get_solver()

    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_ccro_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        solver=watertap_solver,
        outlvl=logging.WARNING,
    )

    mp.n_time_points = n_time_points
    mp.include_costing = include_costing

    operation_mode_list = ["filtration"] * (n_time_points - 1) + ["flushing"]
    flowsheet_options = {
        t: {"time_blk": t, "operation_mode": operation_mode_list[t]}
        for t in range(n_time_points)
    }

    # Build instances of the process model for each time period
    mp.build_multi_period_model(model_data_kwargs=flowsheet_options)

    if include_costing:
        mp.costing = WaterTAPCosting()
        for t, m in enumerate(mp.get_active_process_blocks(), 1):
            if t == 1:
                cc_utils.register_costed_unit(
                    mp, m.fs.P1, register_electricity_flow_only=True
                )
                cc_utils.register_costed_unit(mp, m.fs.RO)
            elif t == n_time_points - 1:
                cc_utils.register_costed_unit(mp, m.fs.P2)
                cc_utils.register_costed_unit(mp, m.fs.P1)
            # Last time period is flushing
            elif t == n_time_points:
                cc_utils.register_costed_unit(
                    mp, m.fs.P2, register_electricity_flow_only=True
                )
                cc_utils.register_costed_unit(
                    mp, m.fs.P1, register_electricity_flow_only=True
                )
    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if t == 1:
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            unfix_dof(
                m, unfix_dead_volume_state=False, cc_configuration=cc_configuration
            )
            old_m = m
        # Last time period is flushing
        if t == n_time_points:
            cc_utils.copy_time_period_links(
                old_m,
                m,
                [
                    {
                        "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
                        "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]',
                    },
                    {
                        "new_model_var": 'fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]',
                        "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]',
                    },
                ],
            )
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            unfix_dof(
                m, unfix_dead_volume_state=False, cc_configuration=cc_configuration
            )
        else:
            cc_utils.copy_state(old_m, m)
            cc_utils.copy_time_period_links(
                old_m,
                m,
                [
                    {
                        "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
                        "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]',
                    },
                    {
                        "new_model_var": 'fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]',
                        "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]',
                    },
                ],
            )

        print(f"Period {t} DOF:", degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0

        results = solve(model=m, tee=True)
        assert_optimal_termination(results)
        unfix_dof(m, unfix_dead_volume_state=True, cc_configuration=cc_configuration)
        old_m = m

    print("--- Completed sequential initialization ---")
    add_multiperiod_variables(mp)
    add_multiperiod_constraints(mp)

    if include_costing:
        print("Adding costing")
        mp.costing.cost_process()
        mp.costing.add_LCOW(mp.avg_product_flow_rate)
        mp.costing.add_specific_energy_consumption(mp.avg_product_flow_rate, name="SEC")

        mp.costing.initialize()

    scale_multiperiod_model(mp)
    calculate_scaling_factors(mp)

    print("Multi-period DOF:", degrees_of_freedom(mp))
    return mp


def build_ccro_system(time_blk=None, operation_mode=None):
    """
    Build the CCRO steady state flowsheet
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.operation_mode = operation_mode

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

    if m.fs.operation_mode == "filtration":

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

        m.fs.RO = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=4,
            module_type="spiral_wound",
            has_full_reporting=True,
        )

        m.fs.product = Product(property_package=m.fs.properties)

        # Add connections
        m.fs.raw_feed_to_P1 = Arc(
            source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet
        )

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

    elif m.fs.operation_mode == "flushing":

        m.fs.P1 = Pump(property_package=m.fs.properties)
        m.fs.P2 = Pump(property_package=m.fs.properties)

        # Add connections
        m.fs.raw_feed_to_P1 = Arc(
            source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet
        )
        m.fs.P1_to_dead_volume = Arc(
            source=m.fs.P1.outlet, destination=m.fs.dead_volume.inlet
        )

        m.fs.dead_volume_to_P2 = Arc(
            source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
        )

        m.fs.flushing = Flushing()

        TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def initialize_system(m, **kwargs):
    """
    Initialize the model by fixing the values of certain variables.
    """

    if m.fs.operation_mode == "filtration":

        propagate_state(m.fs.raw_feed_to_P1)

        m.fs.P1.outlet.pressure[0].fix(
            m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"].value * 2 + 2e5
        )
        m.fs.P2.outlet.pressure[0].fix(
            m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"].value * 2 + 2e5
        )
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

    if m.fs.operation_mode == "flushing":

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

        m.fs.flushing.flushing_time.unfix()
        m.fs.flushing.pre_flushing_concentration.fix()
        m.fs.flushing.flushing_efficiency.fix()
        m.fs.flushing.pre_flushing_concentration.display()
        m.fs.flushing.display()
        m.fs.flushing.initialize()

    return m


def copy_inlet_state_for_mixer(m):
    for idx, obj in m.fs.M1.inlet_2.flow_mass_phase_comp.items():
        obj.value = m.fs.M1.inlet_1.flow_mass_phase_comp[idx].value


def fix_dof_and_initialize(m, cc_configuration=None, **kwargs):
    """
    Fix DOF for MP model and initialize steady-state models.
    """

    set_operating_conditions(m=m, cc_configuration=cc_configuration, **kwargs)
    initialize_system(m=m, **kwargs)


def get_variable_pairs(t1, t2):
    """
    Get variable pairs for connecting two time periods.
    1. dead_volume mass fraction to delta_state
    2. dead_volume density to delta_state
    """

    return [
        (
            t2.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"],
            t1.fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp[
                "Liq", "NaCl"
            ],
        ),
        (
            t2.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"],
            t1.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"],
        ),
    ]


def set_operating_conditions(m, cc_configuration=None, **kwargs):

    if m.fs.operation_mode == "filtration":

        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
        m.fs.raw_feed.properties[0].temperature.fix(293.15)

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
            cc_configuration["raw_feed_flowrate"]
        )
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
            cc_configuration["raw_feed_conc"]
        )

        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        # Dead Volume operating conditions
        # Fix volume
        m.fs.dead_volume.volume.fix(cc_configuration["dead_volume"])
        m.fs.dead_volume.delta_state.volume.fix(cc_configuration["dead_volume"])

        m.fs.dead_volume.accumulation_time.fix(cc_configuration["accumulation_time"])

        # Fixing the flow rate of the dead volume delta state to the feed to calculate the mass fraction and density
        # I found fixing mass fraction and density is easiest way to get initial state
        # we will also use these as connection points between current and future state.

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m.fs.raw_feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value
        )
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
        )

        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(cc_configuration["p1_eff"])

        cc_utils.set_pump_operating_pressure(
            unit=m.fs.P1,
            configuration_options=cc_configuration,
            osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
        )

        # Pump 2 operating conditions
        m.fs.P2.efficiency_pump.fix(cc_configuration["p2_eff"])

        # Set RO configuration parameters
        m.fs.RO.A_comp.fix(cc_configuration["A_comp"])
        m.fs.RO.B_comp.fix(cc_configuration["B_comp"])
        m.fs.RO.area.fix(cc_configuration["membrane_area"])
        m.fs.RO.length.unfix()
        m.fs.RO.width.unfix()
        m.fs.RO.feed_side.velocity[0, 0].fix(0.15)
        m.fs.RO.feed_side.channel_height.fix(cc_configuration["channel_height"])
        m.fs.RO.feed_side.spacer_porosity.fix(cc_configuration["spacer_porosity"])

        m.fs.RO.permeate.pressure[0].fix(101325 * pyunits.Pa)

        m.fs.RO.feed_side.friction_factor_darcy.setub(200)
        # m.fs.RO.flux_mass_phase_comp.setub(1)
        # m.fs.RO.feed_side.cp_modulus.setub(50)
        # m.fs.RO.feed_side.cp_modulus.setlb(0.1)
        m.fs.RO.deltaP.setlb(None)

        scale_filtration_system(m)

    elif m.fs.operation_mode == "flushing":

        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(cc_configuration["p1_eff"])

        cc_utils.set_pump_operating_pressure(
            unit=m.fs.P1,
            configuration_options=cc_configuration,
            osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
        )

        # Pump 2 operating conditions - Add only for costing function. No work is done by this pump
        m.fs.P2.efficiency_pump.fix(cc_configuration["p2_eff"])

        # Concentration of the flushing water is the raw feed concentration
        if cc_configuration["flushing_conc"] == "raw_feed_conc":
            m.fs.flushing.flushing_feed_concentration.fix(
                cc_configuration["raw_feed_conc"]
            )
        else:
            m.fs.flushing.flushing_feed_concentration.fix(
                cc_configuration["flushing_conc"]
            )

        # Calculating values for the dead volume and delta state
        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
        m.fs.raw_feed.properties[0].temperature.fix(cc_configuration["temperature"])

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

        # Reassign the raw feed flowrate and concentration
        if cc_configuration["flushing_flowrate"] == "raw_feed_flowrate":
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                cc_configuration["raw_feed_flowrate"]
            )
        else:
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                cc_configuration["flushing_flowrate"]
            )
        if cc_configuration["flushing_conc"] == "raw_feed_conc":
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                cc_configuration["raw_feed_conc"]
            )
        else:
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                cc_configuration["flushing_conc"]
            )
        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        # Dead Volume operating conditions
        # Fix volume
        m.fs.dead_volume.volume.fix(cc_configuration["dead_volume"])
        m.fs.dead_volume.delta_state.volume.fix(cc_configuration["dead_volume"])
        m.fs.dead_volume.accumulation_time.fix(cc_configuration["accumulation_time"])
        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix()
        m.fs.flushing.flushing_efficiency.fix(cc_configuration["flushing_efficiency"])

        # Constraints
        @m.fs.Constraint()
        def dead_volume_residence_time_constraint(m):
            return m.flushing.mean_residence_time == (
                pyunits.convert(
                    m.dead_volume.volume[0, "Liq"]
                    / m.raw_feed.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.s,
                )
            )

        # Calculate pre-flushing/ dead volume delta state concentration. Concentration before
        # flushing should be the delta state concentration
        @m.fs.Constraint()
        def pre_flushing_conc_constraint(m):
            return (
                m.flushing.pre_flushing_concentration
                == m.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
                * m.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            )

        # Concentration after flushing should be the dead volume properties out concentration
        @m.fs.Constraint()
        def post_flushing_conc_constraint(m):
            return (
                m.dead_volume.dead_volume.mass_phase_comp[0, "Liq", "NaCl"]
                == m.flushing.post_flushing_concentration
                * m.dead_volume.volume[0, "Liq"]
            )

        scale_flushing_system(m)

    return m


def scale_flushing_system(m=None):
    """
    Scale flushing model configuration
    """

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)

    set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 10)
    set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 10)

    constraint_scaling_transform(m.fs.post_flushing_conc_constraint, 1e-1)
    constraint_scaling_transform(m.fs.pre_flushing_conc_constraint, 1e-1)
    set_scaling_factor(m.fs.flushing.pre_flushing_concentration, 1e-1)
    set_scaling_factor(m.fs.flushing.post_flushing_concentration, 1e-1)

    calculate_scaling_factors(m)


def scale_filtration_system(m):
    """
    Scale filtration model configuration
    """

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    set_scaling_factor(m.fs.RO.area, 1e-2)

    set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 10)
    set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 10)

    calculate_scaling_factors(m)


def scale_multiperiod_model(mp):

    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)
    iscale.set_scaling_factor(mp.total_cycle_time, 1e-3)
    iscale.set_scaling_factor(mp.total_filtration_time, 1e-3)
    iscale.set_scaling_factor(mp.final_concentration, 1e-1)

    iscale.constraint_scaling_transform(mp.total_filtration_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.total_cycle_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.final_concentration_constraint, 1e-1)
    iscale.constraint_scaling_transform(mp.global_dead_volume_constraint, 1e2)

    for c in mp.ro_membrane_area_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.ro_membrane_length_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.equal_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_delta_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_recycle_rate.values():
        iscale.constraint_scaling_transform(c, 1e2)


def solve(model=None, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(model)}")
    results = solver.solve(model, tee=tee)
    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)

        raise RuntimeError(msg)
    else:
        print(msg)
        # debug(model, solver=solver, automate_rescale=False, resolve=False)
        # check_jac(model)
        assert False
    return results


def unfix_dof(m, unfix_dead_volume_state=True, cc_configuration=None, **kwargs):
    """
    Unfix linking variables in MP model
    """

    if unfix_dead_volume_state:
        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()

    if m.fs.operation_mode == "filtration":

        m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix(
            cc_configuration["recycle_flowrate"]
        )

    elif m.fs.operation_mode == "flushing":
        m.fs.dead_volume.accumulation_time.unfix()
        m.fs.flushing.flushing_efficiency.fix(cc_configuration["flushing_efficiency"])
        m.fs.flushing.flushing_time.unfix()
        m.fs.flushing.mean_residence_time.unfix()
        m.fs.dead_volume.dead_volume.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.flushing.pre_flushing_concentration.unfix()
        m.fs.flushing.post_flushing_concentration.unfix()


def add_multiperiod_variables(mp):
    """
    Add variables to the multiperiod model.
    """
    mp.dead_volume_to_area_ratio = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3 / pyunits.m**2,
        doc="Global dead volume over all time periods",
    )
    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)

    mp.dead_volume_to_area_multiplier = Var(
        initialize=2,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Global dead volume over all time periods",
    )
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)

    mp.dead_volume_to_area_multiplier.fix(1)
    mp.dead_volume_to_area_ratio.fix(
        value(1 * pyunits.m * (3.14 * 0.1016**2) * pyunits.m**2 / (7.2 * pyunits.m**2))
    )

    # Permeate and feed
    mp.overall_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Overall water recovery over all time periods",
    )
    mp.apparent_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Apparent water recovery over all time periods",  # from BLM report on OCWD pilot
    )

    mp.final_concentration = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.L,
        doc="Final concentration of the product stream",
    )

    mp.total_permeate = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total permeate over all time periods",
    )

    mp.total_feed = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total feed over all time periods",
    )

    # mp.total_brine = Var(
    #     initialize=0.5,
    #     domain=NonNegativeReals,
    #     units=pyunits.m**3,
    #     doc="Total brine over all time periods",
    # )

    mp.avg_product_flow_rate = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3 / pyunits.s,
        doc="Average permeate production over all time periods",
    )

    mp.total_filtration_time = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.s,
        doc="Total filtration time excluding flushing",
    )

    mp.total_cycle_time = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.s,
        doc="Total cycle time including flushing",
    )

    mp.cycle_time_ratio = Var(
        initialize=1,
        domain=NonNegativeReals,
        bounds=(0, 1.0001),
        units=pyunits.dimensionless,
        doc="Ratio of total cycle time to filtration time",
    )


def add_multiperiod_constraints(mp):
    """
    Add constraints to the multiperiod model.
    """

    # Get all filtration time blocks
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]

    # RO membrane area should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane area equality through all filtration periods",
    )
    def ro_membrane_area_constraint(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return blks[t].fs.RO.area == b0.fs.RO.area

    mp.ro_membrane_area_constraint.deactivate()

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane length equality through all filtration periods",
    )
    def ro_membrane_length_constraint(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return blks[t].fs.RO.length == b0.fs.RO.length

    mp.ro_membrane_length_constraint.deactivate()

    # Dead volume should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Dead volume equality through all filtration periods",
    )
    def equal_dead_volume_constraint(b, t):
        if t == b.TIME.first():
            return Constraint.Skip
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == b0.fs.dead_volume.volume[0, "Liq"]
        )

    mp.equal_dead_volume_constraint.deactivate()

    @mp.Constraint(
        mp.TIME, doc="Delta dead volume equality through all filtration periods"
    )
    def equal_delta_dead_volume_constraint(b, t):
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == blks[t].fs.dead_volume.delta_state.volume[0, "Liq"]
        )

    mp.equal_delta_dead_volume_constraint.deactivate()

    # Recycle rate should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Recycle rate equality through all filtration periods",
    )
    def equal_recycle_rate(b, t):
        if t in [b.TIME.first(), b.TIME.last()]:
            return Constraint.Skip
        return (
            blks[t].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
            == blks[0].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
        )

    mp.equal_recycle_rate.deactivate()

    @mp.Constraint(doc="Global dead volume constraint")
    def global_dead_volume_constraint(b):
        return (
            b0.fs.dead_volume.volume[0, "Liq"]
            == b0.fs.RO.area
            * b.dead_volume_to_area_ratio
            * b.dead_volume_to_area_multiplier
        )

    calculate_variable_from_constraint(
        b0.fs.dead_volume.volume[0, "Liq"], mp.global_dead_volume_constraint
    )

    mp.global_dead_volume_constraint.deactivate()

    # Density at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Density equality between end of flushing and start of filtration"
    )
    def cycle_end_density_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .dens_mass_phase["Liq"]
        )

    # Mass fraction at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Mass fraction equality between end of flushing and start of filtration"
    )
    def cycle_end_mass_frac_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
            == blks[-1]
            .fs.dead_volume.dead_volume.properties_out[0]
            .mass_frac_phase_comp["Liq", "NaCl"]
        )

    @mp.Constraint(doc="Total filtration time constraint")
    def total_filtration_time_constraint(b):
        return b.total_filtration_time == sum(
            blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(b.n_time_points - 1)
        )

    @mp.Constraint(doc="Total cycle time constraint")
    def total_cycle_time_constraint(b):
        return (
            b.total_cycle_time
            == b.total_filtration_time + blks[-1].fs.flushing.flushing_time
        )

    @mp.Constraint(doc="Cycle time ratio constraint")
    def cycle_time_ratio_constraint(b):
        return b.cycle_time_ratio == (b.total_filtration_time / b.total_cycle_time)

    @mp.Constraint(doc="Final concentration constraint")
    def final_concentration_constraint(b):
        return b.final_concentration == blks[-1].fs.flushing.pre_flushing_concentration

    # Total permeate
    @mp.Constraint(doc="Total permeate produced over all time periods")
    def total_permeate_constraint(b):
        return b.total_permeate == sum(
            blks[t].fs.product.properties[0].flow_vol_phase["Liq"]
            * blks[t].fs.dead_volume.accumulation_time[0]
            for t in range(b.n_time_points - 1)
        )

    # # Total brine -> Convert to expression
    # @mp.Constraint(doc="Total brine produced over all time periods")
    # def total_brine_constraint(mp):
    #     blks = list(mp.get_active_process_blocks())
    #     return mp.total_brine == sum(
    #         blks[t].fs.brine.properties[0].flow_vol_phase["Liq"]
    #         * blks[t].fs.dead_volume.accumulation_time[0]
    #         for t in range(n_time_points - 1)
    #     )

    @mp.Constraint(doc="Average product flow rate over all time periods")
    def eq_avg_product_flow_rate(mp):
        blks = mp.get_active_process_blocks()
        return mp.avg_product_flow_rate == pyunits.convert(
            mp.total_permeate
            / (
                sum(
                    blks[t].fs.dead_volume.accumulation_time[0]
                    for t in range(len(blks) - 1)
                )
                + blks[-1].fs.flushing.flushing_time
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    # Total feed
    @mp.Constraint(doc="Total feed volume over all time periods")
    def total_feed_constraint(b):
        return (
            b.total_feed
            == sum(
                blks[t].fs.raw_feed.properties[0].flow_vol_phase["Liq"]
                * blks[t].fs.dead_volume.accumulation_time[0]
                for t in range(b.n_time_points - 1)
            )
            + b0.fs.dead_volume.volume[0, "Liq"]
        )

    # Overall water recovery
    @mp.Constraint(doc="Overall water recovery for system")
    def overall_water_recovery_constraint(b):
        return b.total_permeate == b.overall_recovery * b.total_feed


def fix_overall_water_recovery(mp, overall_water_recovery):

    mp.overall_recovery.fix(overall_water_recovery)

    # Fixed for accumulation time for initialization
    for t, m in enumerate(mp.get_active_process_blocks()):
        m.fs.dead_volume.accumulation_time.unfix()
        # m.fs.dead_volume.accumulation_time.setlb(1)
        # m.fs.dead_volume.accumulation_time.setub(400)

        set_scaling_factor(m.fs.dead_volume.accumulation_time, 1e-2)

    # Equal accumulation time across all filtration periods
    @mp.Constraint(list(range(1, mp.n_time_points - 1)))
    def accumulation_time_cons(mp, t):
        blks = list(mp.get_active_process_blocks())
        return blks[t].fs.dead_volume.accumulation_time[0] == (
            blks[0].fs.dead_volume.accumulation_time[0]
        )


def setup_optimization(
    mp, overall_water_recovery=0.5, max_cycle_time_hr=1, recycle_flow_bounds=(0.1, 20)
):
    """
    Setup the multiperiod model for optimization.
    """
    fix_overall_water_recovery(mp, overall_water_recovery)
    mp.global_dead_volume_constraint.activate()
    mp.equal_dead_volume_constraint.activate()
    mp.equal_delta_dead_volume_constraint.activate()
    mp.ro_membrane_area_constraint.activate()
    mp.ro_membrane_length_constraint.activate()
    mp.total_cycle_time.setub(max_cycle_time_hr * pyunits.hours)
    # mp.total_cycle_time.setlb(0.98 * max_cycle_time_hr * pyunits.hours)
    mp.equal_recycle_rate.activate()
    print("DOF for optimization:", degrees_of_freedom(mp))
    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if m.fs.find_component("RO") is not None:
            m.fs.RO.length.unfix()

            m.fs.RO.width.unfix()
            m.fs.RO.length.setlb(0.1)

            m.fs.RO.width.setlb(0.1)

            # m.fs.RO.mixed_permeate[0].conc_mass_phase_comp.display()
            m.fs.RO.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].setub(0.5)
            m.fs.RO.area.setlb(50)
            m.fs.RO.area.unfix()
            m.fs.RO.flux_mass_phase_comp[0.0, 1.0, "Liq", "H2O"].setlb(
                0.001 * pyunits.kg / pyunits.m**2 / pyunits.hr
            )
            m.fs.RO.feed_side.velocity[0, 0].setub(0.3)
            m.fs.RO.feed_side.velocity[0, 0].setlb(0.05)
            m.fs.RO.feed_side.velocity[0, 0].unfix()
        m.fs.dead_volume.volume.unfix()
        m.fs.dead_volume.delta_state.volume[0, "Liq"].unfix()

        if t != mp.n_time_points:
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
                recycle_flow_bounds[0] * pyunits.L / pyunits.s
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
                recycle_flow_bounds[1] * pyunits.L / pyunits.s
            )

    if mp.include_costing:
        mp.cost_objective = Objective(expr=mp.costing.LCOW)

    print("DOF for optimization:", degrees_of_freedom(mp))


def print_results_table(mp, w=15):
    """
    Print multiperiod CCRO results in a tabular format in the terminal.
    w = column width
    """
    n = 13  # number of columns
    title = "CCRO MULTIPERIOD RESULTS"
    side = int(((n * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    # print("\n" + "=" * 120)
    # print("CCRO MULTIPERIOD RESULTS")
    print(f"{'-' * (n * w)}")

    # Header
    print(
        f"{'Period':<{w}s}{'Acc Time':<{w}s}{'Raw Feed':<{w}s}{'Permeate':<{w}s}{'SP Recovery':<{w}s}{'P1':<{w}s}{'P2':<{w}s}{'RO In':<{w}s}{'RO In':<{w}s}{'Dead Vol In':<{w}s}{'Dead Vol In':<{w}s}{'Delta State':<{w}s}{'Dead Vol':<{w}s}"
    )
    print(
        f"{'':<{w}s}{'(s)':<{w}s}{'(L/min)':<{w}s}{'(L/min)':<{w}s}{'(%)':<{w}s}{'(kW)':<{w}s}{'(kW)':<{w}s}{'(Pa)':<{w}s}{'(Psi)':<{w}s}{'(L/min)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}"
    )
    print(f"{'-' * (n * w)}")

    # Data rows
    for t, blks in enumerate(mp.get_active_process_blocks(), 1):
        # if t == len(blks):
        if blks.fs.operation_mode == "flushing":
            accumulation_time = blks.fs.flushing.flushing_time.value
        else:
            accumulation_time = blks.fs.dead_volume.accumulation_time[0].value

        raw_feed = pyunits.convert(
            blks.fs.raw_feed.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.L / pyunits.min,
        )()

        if blks.fs.find_component("product") is not None:
            permeate = pyunits.convert(
                blks.fs.product.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.L / pyunits.min,
            )()
        else:
            permeate = 0
        if blks.fs.operation_mode == "filtration":

            p2_out = value(
                pyunits.convert(blks.fs.P2.work_mechanical[0], to_units=pyunits.kW)
            )
            p1_out = value(
                pyunits.convert(blks.fs.P1.work_mechanical[0], to_units=pyunits.kW)
            )
            ro_pressure = blks.fs.M1.mixed_state[0].pressure.value
            ro_pressure_psi = value(
                pyunits.convert(
                    blks.fs.M1.mixed_state[0].pressure, to_units=pyunits.psi
                )
            )
            sp_recovery = blks.fs.RO.recovery_vol_phase[0.0, "Liq"].value

        if blks.fs.operation_mode == "flushing":
            mixer_out = value(
                pyunits.convert(
                    blks.fs.raw_feed.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.L / pyunits.min,
                )
            )
            sp_recovery = 0
            ro_pressure = 0

        dead_vol_in = pyunits.convert(
            blks.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyunits.L / pyunits.min,
        )()
        dead_vol_in_conc = (
            blks.fs.dead_volume.dead_volume.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )

        delta_state_conc = blks.fs.dead_volume.delta_state.conc_mass_phase_comp[
            "Liq", "NaCl"
        ]()
        dead_vol_out_conc = (
            blks.fs.dead_volume.dead_volume.properties_out[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )

        print(
            f"{t:<{w}d}{accumulation_time:<{w}.2f}{raw_feed:<{w}.6f}{permeate:<{w}.6f}{sp_recovery:<{w}.6f}{p1_out:<{w}.6f}{p2_out:<{w}.6f}{ro_pressure:<{w}.2f}{ro_pressure_psi:<{w}.2f}{dead_vol_in:<{w}.6f}{dead_vol_in_conc:<{w}.6f}{delta_state_conc:<{w}.6f}{dead_vol_out_conc:<{w}.6f}"
        )

    print(f"{'=' * (n * w)}")
    w = w * 2

    print(
        f"{'Total cycle time (s):':<{w}s}",
        mp.total_cycle_time.value,
    )
    print(
        f"{'Mean residence time (s):':<{w}s}",
        mp.get_active_process_blocks()[-1].fs.flushing.mean_residence_time.value,
    )
    print(
        f"{'Flushing time (s):':<{w}s}",
        mp.get_active_process_blocks()[-1].fs.flushing.flushing_time.value,
    )
    print(
        f"{'Flushing efficiency:':<{w}s}",
        mp.get_active_process_blocks()[-1].fs.flushing.flushing_efficiency.value,
    )
    print(
        f"{'Pre-flushing conc (kg/m3):':<{w}s}",
        mp.get_active_process_blocks()[-1].fs.flushing.pre_flushing_concentration.value,
    )
    print(
        f"{'Post-flushing conc (kg/m3):':<{w}s}",
        mp.get_active_process_blocks()[
            -1
        ].fs.flushing.post_flushing_concentration.value,
    )

    print(f"\n{'Overall recovery:':<{w}s}", mp.overall_recovery.value)
    print(f"{'Total feed (m3):':<{w}s}", mp.total_feed.value)
    print(f"{'Total permeate (m3):':<{w}s}", mp.total_permeate.value)

    print(
        f"{'Density at start of cycle:':<{w}s}",
        mp.get_active_process_blocks()[0]
        .fs.dead_volume.dead_volume.properties_out[0]
        .dens_mass_phase["Liq"]
        .value,
    )
    print(
        f"{'Density at the end of flushing:':<{w}s}",
        mp.get_active_process_blocks()[-1]
        .fs.dead_volume.dead_volume.properties_out[0]
        .dens_mass_phase["Liq"]
        .value,
    )
    print(
        f"{'Dead Volume:':<{w}s}",
        mp.get_active_process_blocks()[0].fs.dead_volume.volume[0, "Liq"].value,
    )

    print(
        f"{'Membrane Area:':<{w}s}", mp.get_active_process_blocks()[0].fs.RO.area.value
    )
    print(
        f"{'Membrane Length:':<{w}s}",
        mp.get_active_process_blocks()[0].fs.RO.length.value,
    )

    print(
        f"{'Membrane Width:':<{w}s}",
        mp.get_active_process_blocks()[0].fs.RO.width.value,
    )
    print(
        f"{'RO Inlet Velocity:':<{w}s}",
        mp.get_active_process_blocks()[0].fs.RO.feed_side.velocity[0, 0].value,
    )


if __name__ == "__main__":
    cc_config = CCROConfiguration()
    cc_config["raw_feed_conc"] = 5 * pyunits.g / pyunits.L
    mp = create_ccro_multiperiod(
        n_time_points=6, include_costing=True, cc_configuration=cc_config
    )
    # results = solve(mp, tee=False)
    # print_results_table(mp, w=15)

    setup_optimization(
        mp,
        overall_water_recovery=0.8,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
    )

    results = solve(mp)
    print_results_table(mp, w=16)
