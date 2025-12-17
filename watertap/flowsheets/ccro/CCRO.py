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
)

import watertap.flowsheets.ccro.utils.utils as cc_utils
from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration


from watertap.flowsheets.ccro.ccro_flowsheet_functions import (
    multi_period_constraints as ccro_mp_constraints,
    scaling as ccro_scaling,
    operating_conditions as ccro_operating_conditions,
)

import watertap.flowsheets.ccro.utils.ipoptv2 as ipt2


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
        unfix_dof_func=ccro_operating_conditions.unfix_dof,
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
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.RO,
                    register_electricity_cost=False,
                    register_capital_cost=True,
                )
            elif t == n_time_points - 1:
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=True,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=True,
                )
            # Last time period is flushing
            elif t == n_time_points:
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                )
    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if t == 1:
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            ccro_operating_conditions.unfix_dof(
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
            ccro_operating_conditions.unfix_dof(
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
        ccro_operating_conditions.unfix_dof(
            m, unfix_dead_volume_state=True, cc_configuration=cc_configuration
        )
        old_m = m

    print("--- Completed sequential initialization ---")
    ccro_mp_constraints.add_multiperiod_variables(mp, cc_configuration=cc_configuration)
    ccro_mp_constraints.add_multiperiod_constraints(
        mp, cc_configuration=cc_configuration
    )

    if include_costing:
        print("Adding costing")
        mp.costing.cost_process()
        mp.costing.add_LCOW(mp.avg_product_flow_rate)
        mp.costing.add_specific_energy_consumption(mp.avg_product_flow_rate, name="SEC")

        mp.costing.initialize()

    ccro_scaling.scale_multiperiod_model(mp)
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
            finite_elements=10,
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
        m.fs.flushing.flushing_time = value(m.fs.flushing.mean_residence_time)
        m.fs.flushing.pre_flushing_concentration.fix()
        m.fs.flushing.post_flushing_concentration.fix(
            value(m.fs.flushing.pre_flushing_concentration)
        )
        m.fs.flushing.post_flushing_concentration.unfix()
        m.fs.flushing.flushing_efficiency.fix()
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
    ccro_operating_conditions.set_operating_conditions(
        m=m, cc_configuration=cc_configuration, **kwargs
    )
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


def solve(
    model=None,
    solver=None,
    tee=True,
    raise_on_failure=True,
    use_ipoptv2=False,
):
    # ---solving---
    if solver is None and use_ipoptv2 == False:
        solver = get_solver()
    elif solver is None and use_ipoptv2:

        solver = ipt2.get_solver()
        # for constraint in model.component_data_objects(Constraint, active=True):
        #     sc = iscale.get_constraint_transform_applied_scaling_factor(constraint)
        #     if sc is not None:
        #         print(f"Constraint: {constraint.name}, Scaling: {sc}")
        #         iscale.constraint_scaling_transform_undo(constraint)
        #         iscale.set_scaling_factor(constraint, sc)
        # else:
        #     print(f"Constraint: {constraint.name}, No scaling factor found")

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


def setup_optimization(
    mp, overall_water_recovery=0.5, max_cycle_time_hr=1, recycle_flow_bounds=(0.1, 20)
):
    """
    Setup the multiperiod model for optimization.
    """
    ccro_mp_constraints.fix_overall_water_recovery(mp, overall_water_recovery)
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
        f"{'Period':<{w}s}{'Acc Time':<{w}s}{'Raw Feed':<{w}s}{'Permeate':<{w}s}{'Raw Feed':<{w}s}{'Permeate':<{w}s}{'SP Recovery':<{w}s}{'P1':<{w}s}{'P2':<{w}s}{'RO In':<{w}s}{'RO In':<{w}s}{'Dead Vol In':<{w}s}{'Dead Vol In':<{w}s}{'Dead Vol In':<{w}s}{'Delta State':<{w}s}{'Dead Vol':<{w}s}{'Dead dens delta':<{w}s}{'Dead dens in':<{w}s}{'Dead dens out':<{w}s}"
    )
    print(
        f"{'':<{w}s}{'(s)':<{w}s}{'(L/min)':<{w}s}{'(L/min)':<{w}s}{'(kg/min)':<{w}s}{'(kg/min)':<{w}s}{'(%)':<{w}s}{'(kW)':<{w}s}{'(kW)':<{w}s}{'(Pa)':<{w}s}{'(Psi)':<{w}s}{'(L/min)':<{w}s}{'(kg/min)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m3)':<{w}s}{'(kg/m3)':<{w}s}{'(kg/m3)':<{w}s}"
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
        raw_feed_mass = sum(
            pyunits.convert(
                blks.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", comp],
                to_units=pyunits.kg / pyunits.min,
            )()
            for comp in blks.fs.properties.component_list
        )
        if blks.fs.find_component("product") is not None:
            permeate = pyunits.convert(
                blks.fs.product.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.L / pyunits.min,
            )()
            permeate_mass = sum(
                pyunits.convert(
                    blks.fs.product.properties[0].flow_mass_phase_comp["Liq", comp],
                    to_units=pyunits.kg / pyunits.min,
                )()
                for comp in blks.fs.properties.component_list
            )
        else:
            permeate_mass = 0
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
        dead_vol_mass = sum(
            pyunits.convert(
                blks.fs.dead_volume.dead_volume.properties_in[0].flow_mass_phase_comp[
                    "Liq", comp
                ],
                to_units=pyunits.kg / pyunits.min,
            )()
            for comp in blks.fs.properties.component_list
        )
        dead_vol_in_conc = (
            blks.fs.dead_volume.dead_volume.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )
        dead_vol_dense_in = (
            blks.fs.dead_volume.dead_volume.properties_in[0]
            .dens_mass_phase["Liq"]
            .value
        )
        dead_vol_dense_out = (
            blks.fs.dead_volume.dead_volume.properties_out[0]
            .dens_mass_phase["Liq"]
            .value
        )
        dead_vol_dense_delta = blks.fs.dead_volume.delta_state.dens_mass_phase[
            0, "Liq"
        ].value
        delta_state_conc = blks.fs.dead_volume.delta_state.conc_mass_phase_comp[
            "Liq", "NaCl"
        ]()
        dead_vol_out_conc = (
            blks.fs.dead_volume.dead_volume.properties_out[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )

        print(
            f"{t:<{w}d}{accumulation_time:<{w}.2f}{raw_feed:<{w}.6f}{permeate:<{w}.6f}{raw_feed_mass:<{w}.6f}{permeate_mass:<{w}.6f}{sp_recovery:<{w}.6f}{p1_out:<{w}.6f}{p2_out:<{w}.6f}{ro_pressure:<{w}.2f}{ro_pressure_psi:<{w}.2f}{dead_vol_in:<{w}.6f}{dead_vol_mass:<{w}.6f}{dead_vol_in_conc:<{w}.6f}{delta_state_conc:<{w}.6f}{dead_vol_out_conc:<{w}.6f}{dead_vol_dense_delta:<{w}.6f}{dead_vol_dense_in:<{w}.6f}{dead_vol_dense_out:<{w}.6f}"
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
    from idaes.core.util.model_diagnostics import DiagnosticsToolbox

    cc_config = CCROConfiguration()
    cc_config["raw_feed_conc"] = 35 * pyunits.g / pyunits.L

    cc_config["accumulation_time"] = 1 * pyunits.second
    mp = create_ccro_multiperiod(
        n_time_points=51,
        include_costing=True,
        cc_configuration=cc_config,
    )
    setup_optimization(
        mp,
        overall_water_recovery=0.5,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
    )

    db = DiagnosticsToolbox(mp)
    db.display_constraints_with_large_residuals()
    # assert False
    results = solve(mp, use_ipoptv2=False)
    mp.overall_recovery.fix(0.55)
    results = solve(mp, use_ipoptv2=False)
    print_results_table(mp, w=16)
    blks = list(mp.get_active_process_blocks())
    blks[0].fs.raw_feed.properties[0].display()
    blks[-1].fs.raw_feed.properties[0].display()
