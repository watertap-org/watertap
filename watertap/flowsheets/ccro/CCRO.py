from pyomo.environ import (
    check_optimal_termination,
    value,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver
from watertap.costing import (
    WaterTAPCosting,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
import watertap.flowsheets.ccro.utils.utils as cc_utils
import watertap.flowsheets.ccro.ccro_flowsheet_functions.unit_operations as unit_operations
from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from watertap.flowsheets.ccro.ccro_flowsheet_functions import (
    multi_period_constraints as ccro_mp_constraints,
    scaling as ccro_scaling,
    operating_conditions as ccro_operating_conditions,
)
import watertap.flowsheets.ccro.utils.ipoptv2 as ipt2


def create_ccro_multiperiod(
    n_time_points=5,
    n_flushing_points=5,
    include_costing=True,
    cc_configuration=None,
    use_ro_with_hold_up=False,
):
    """
    Create multiperiod model for CCRO system
    """
    if isinstance(cc_configuration, CCROConfiguration) is False:
        user_config = cc_configuration
        cc_configuration = CCROConfiguration()
        cc_configuration.update(user_config)

    watertap_solver = get_solver()
    print(n_flushing_points + n_flushing_points)
    mp = MultiPeriodModel(
        n_time_points=n_time_points + n_flushing_points,
        process_model_func=build_ccro_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=ccro_operating_conditions.unfix_dof,
        solver=watertap_solver,
        outlvl=logging.WARNING,
    )

    mp.time_points = n_time_points
    mp.flushing_points = n_flushing_points
    mp.include_costing = include_costing
    operation_mode_list = []
    for n in range(n_time_points):
        operation_mode_list.append("filtration")
    if n_flushing_points == 1:
        operation_mode_list.append("flushing")
    else:
        for n in range(n_flushing_points):
            operation_mode_list.append("flushing_with_filtration")
    flowsheet_options = {
        t: {
            "time_blk": t,
            "operation_mode": operation_mode_list[t],
            "use_ro_with_hold_up": use_ro_with_hold_up,
        }
        for t in range(n_time_points + n_flushing_points)
    }
    for i, item in flowsheet_options.items():
        print(i, item)
    # Build instances of the process model for each time period
    mp.build_multi_period_model(model_data_kwargs=flowsheet_options)
    print("Built multi-period model with operation modes:")
    for blk in mp.get_active_process_blocks():
        print(blk.name, blk.fs.operation_mode)

    ccro_mp_constraints.add_multiperiod_variables(mp, cc_configuration=cc_configuration)
    if mp.flushing_points > 1:
        unit_operations.build_flushing_unit_only(mp, start_period=mp.time_points)
        unit_operations.build_conduit(mp)

    if include_costing:
        mp.costing = WaterTAPCosting()
        for t, m in enumerate(mp.get_active_process_blocks(), 1):
            # this will track fraction of power used over the cycle.
            if m.fs.operation_mode == "filtration":
                utilization_ratio = (
                    m.fs.dead_volume.accumulation_time[0] / mp.total_cycle_time
                )
            elif m.fs.operation_mode == "flushing_with_filtration":
                utilization_ratio = (
                    mp.flushing.flushing_time / mp.flushing_points
                ) / mp.total_cycle_time
            elif m.fs.operation_mode == "flushing":
                utilization_ratio = m.fs.flushing.flushing_time / mp.total_cycle_time
            if t == 1:
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                    utilization_factor=utilization_ratio,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                    utilization_factor=utilization_ratio,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.RO,
                    register_electricity_cost=False,
                    register_capital_cost=True,
                    utilization_factor=utilization_ratio,
                )
            elif (
                t == n_time_points - 1
            ):  # Last operating pressure - assume its highest !
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=True,
                    utilization_factor=utilization_ratio,
                    # costing_method_arguments={"pump_type": "low_pressure"},
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=True,
                    utilization_factor=utilization_ratio,
                )
            # Last time period is flushing
            elif t >= n_time_points:
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P2,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                    utilization_factor=utilization_ratio,
                )
                cc_utils.register_costed_unit(
                    mp,
                    m.fs.P1,
                    register_electricity_cost=True,
                    register_capital_cost=False,
                    utilization_factor=utilization_ratio,
                )

        if mp.flushing_points > 1:
            cc_utils.register_costed_unit(
                mp,
                mp.conduit,
                register_electricity_cost=False,
                register_capital_cost=True,
            )
    for t, m in enumerate(mp.get_active_process_blocks(), 1):

        # Last time period is flushing
        print(f"Period {t}:", n_time_points, n_flushing_points, m.fs.operation_mode)
        if t == 1:
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            ccro_operating_conditions.unfix_dof(
                m, unfix_dead_volume_state=False, cc_configuration=cc_configuration
            )
            old_m = m
        if (
            t == (n_time_points + n_flushing_points) and n_flushing_points == 1
        ):  # if we only have one flushing period its just dead volume being flushed!
            base_links = [
                {
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
                    "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]',
                },
                {
                    "new_model_var": 'fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]',
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]',
                },
            ]

            cc_utils.copy_time_period_links(old_m, m, base_links)
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            ccro_operating_conditions.unfix_dof(
                m, unfix_dead_volume_state=False, cc_configuration=cc_configuration
            )
        elif t > n_time_points and n_flushing_points > 1 and t == n_time_points + 1:
            base_links = [
                {
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
                    "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]',
                },
                {
                    "new_model_var": 'fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]',
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]',
                },
            ]
            if m.fs.ro_model_with_hold_up:
                idx = m.fs.RO.difference_elements
                for i in idx:
                    base_links.append(
                        {
                            "old_model_var": f'fs.RO.feed_side.properties[0, {i}].mass_frac_phase_comp["Liq", "NaCl"]',
                            "new_model_var": f'fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[0, {i}, "Liq", "NaCl" ]',
                        },
                    )
                    base_links.append(
                        {
                            "old_model_var": f'fs.RO.feed_side.properties[0, {i}].dens_mass_phase["Liq"]',
                            "new_model_var": f'fs.RO.feed_side.delta_state.node_dens_mass_phase[0, {i}, "Liq"]',
                        },
                    )
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            cc_utils.copy_time_period_links(old_m, m, base_links)
            # ccro_operating_conditions.unfix_dof(
            #     m, unfix_dead_volume_state=False, cc_configuration=cc_configuration
            # )
        else:
            base_links = [
                {
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
                    "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]',
                },
                {
                    "new_model_var": 'fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]',
                    "old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]',
                },
            ]
            if m.fs.ro_model_with_hold_up:
                idx = m.fs.RO.difference_elements
                for i in idx:
                    base_links.append(
                        {
                            "old_model_var": f'fs.RO.feed_side.properties[0, {i}].mass_frac_phase_comp["Liq", "NaCl"]',
                            "new_model_var": f'fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[0, {i}, "Liq", "NaCl" ]',
                        },
                    )
                    base_links.append(
                        {
                            "old_model_var": f'fs.RO.feed_side.properties[0, {i}].dens_mass_phase["Liq"]',
                            "new_model_var": f'fs.RO.feed_side.delta_state.node_dens_mass_phase[0, {i}, "Liq"]',
                        },
                    )
            print("copying time period links...")
            cc_utils.copy_state(old_m, m)
            cc_utils.copy_time_period_links(old_m, m, base_links)
        print(f"Period {t} DOF:", degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0

        results = solve(model=m, tee=True)
        print_close_to_bounds(m)
        assert_optimal_termination(results)
        ccro_operating_conditions.unfix_dof(
            m, unfix_dead_volume_state=True, cc_configuration=cc_configuration
        )
        old_m = m
    if mp.flushing_points > 1:
        unit_operations.initialize_flushing(mp)
    print("--- Completed sequential initialization ---")
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
    for m in mp.get_active_process_blocks():
        if m.fs.find_component("RO") is not None:
            print(
                m.name,
                "outle TDS:",
                m.fs.RO.feed_side.properties[0, 1]
                .conc_mass_phase_comp["Liq", "NaCl"]
                .value,
            )
    print("Multi-period DOF:", degrees_of_freedom(mp))
    return mp


def build_ccro_system(time_blk=None, operation_mode=None, use_ro_with_hold_up=True):
    """
    Build the CCRO steady state flowsheet
    """
    if operation_mode == "filtration":
        m = unit_operations.build_ro_systems(use_ro_with_hold_up=use_ro_with_hold_up)

    elif operation_mode == "flushing":
        m = unit_operations.build_flushing_unit()
        m.fs.ro_model_with_hold_up = (
            use_ro_with_hold_up  # need to add for volume updates
        )

    elif operation_mode == "flushing_with_filtration":
        m = unit_operations.build_flushing_with_RO(
            use_ro_with_hold_up=use_ro_with_hold_up
        )
    else:
        raise (ValueError("Invalid operation mode selected."))
    return m


def initialize_system(m, **kwargs):
    """
    Initialize the model by fixing the values of certain variables.
    """

    if m.fs.operation_mode == "filtration":
        unit_operations.intialize_ro_systems(m)
    if m.fs.operation_mode == "flushing":
        unit_operations.initialize_flushing_unit(m)
    if m.fs.operation_mode == "flushing_with_filtration":
        unit_operations.initialize_flushing_with_RO(m)

    return m


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
    pair_list = []
    if t1.fs.ro_model_with_hold_up:
        idx = t1.fs.RO.difference_elements
        if (
            t2.fs.operation_mode == "filtration"
            or t2.fs.operation_mode == "flushing_with_filtration"
        ):
            for i in idx:
                pair_list.append(
                    (
                        t2.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                            0, i, "Liq", "NaCl"
                        ],
                        t1.fs.RO.feed_side.properties[0, i].mass_frac_phase_comp[
                            "Liq", "NaCl"
                        ],
                    ),
                )
                pair_list.append(
                    (
                        t2.fs.RO.feed_side.delta_state.node_dens_mass_phase[
                            0, i, "Liq"
                        ],
                        t1.fs.RO.feed_side.properties[0, i].dens_mass_phase["Liq"],
                    ),
                )
    pair_list.append(
        (
            t2.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"],
            t1.fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp[
                "Liq", "NaCl"
            ],
        ),
    )
    pair_list.append(
        (
            t2.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"],
            t1.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"],
        ),
    )
    return pair_list


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
    print_close_to_bounds(model)
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)

        # raise RuntimeError(msg)
    return results


def setup_optimization(
    mp, overall_water_recovery=0.5, max_cycle_time_hr=1, recycle_flow_bounds=(0.1, 20)
):
    """
    Setup the multiperiod model for optimization.
    """
    ccro_mp_constraints.fix_overall_water_recovery(mp, overall_water_recovery)
    mp.global_dead_volume_constraint.activate()
    if mp.find_component("global_ro_volume_constraint") is not None:
        mp.global_ro_volume_constraint.activate()
        mp.equal_ro_volume_constraint.activate()
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
            m.fs.RO.area.setlb(50)
            m.fs.RO.area.unfix()
            m.fs.RO.flux_mass_phase_comp[0.0, 1.0, "Liq", "H2O"].setlb(
                0.001 * pyunits.kg / pyunits.m**2 / pyunits.hr
            )
            m.fs.RO.feed_side.velocity[0, 0].setub(0.3)
            m.fs.RO.feed_side.velocity[0, 0].setlb(0.05)
            m.fs.RO.feed_side.velocity[0, 0].unfix()
            if m.fs.ro_model_with_hold_up:
                for i in m.fs.RO.difference_elements:
                    m.fs.RO.feed_side.volume.unfix()
                    m.fs.RO.feed_side.accumulation_time.unfix()
                    # should be redundunt
                    m.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                        0, i, "Liq", "NaCl"
                    ].unfix()
                    m.fs.RO.feed_side.delta_state.node_dens_mass_phase[
                        0, i, "Liq"
                    ].unfix()
        m.fs.dead_volume.volume.unfix()
        m.fs.dead_volume.delta_state.volume[0, "Liq"].unfix()
        m.fs.dead_volume.accumulation_time[0].unfix()
        if m.fs.operation_mode == "filtration":
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
                recycle_flow_bounds[0] * pyunits.L / pyunits.s
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
                recycle_flow_bounds[1] * pyunits.L / pyunits.s
            )
        elif m.fs.operation_mode == "flushing_with_filtration":
            m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"].unfix()
            m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"].setlb(
                recycle_flow_bounds[0] * pyunits.L / pyunits.s
            )
            m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"].setub(
                recycle_flow_bounds[1] * pyunits.L / pyunits.s
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setlb(
                recycle_flow_bounds[0] * pyunits.L / pyunits.s
            )
            m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].setub(
                recycle_flow_bounds[1] * pyunits.L / pyunits.s
            )
        elif m.fs.operation_mode == "flushing":
            m.fs.flushing.flushing_efficiency.unfix()
            m.fs.flushing.flushing_efficiency.setub(0.95)
            m.fs.flushing.flushing_efficiency.setlb(0.1)
    if mp.find_component("flushing") is not None:
        mp.flushing.flushing_efficiency.unfix()
        mp.flushing.flushing_efficiency.setub(0.95)
        mp.flushing.flushing_efficiency.setlb(0.1)
    if mp.find_component("conduit") is not None:
        mp.conduit.volume.unfix()
    if mp.include_costing:
        mp.cost_objective = Objective(expr=mp.costing.LCOW)

    print("DOF for optimization:", degrees_of_freedom(mp))


def fix_optimization_dofs(
    mp,
    accumulation_time=None,
    overal_water_recovery=None,
    add_water_recovery_objective=False,
    add_initial_pressure_objective=False,
    membrane_area=None,
    membrane_length=None,
    recycle_rate=None,
    flushing_efficiency=None,
    target_pressure=245 * pyunits.psi,
):
    blk0 = mp.get_active_process_blocks()[0]
    blkfs = mp.get_active_process_blocks()[-1]
    blk0.fs.RO.length.fix()
    blk0.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix()
    if membrane_area is not None:
        blk0.fs.RO.area.fix(membrane_area)
    if membrane_length is not None:
        blk0.fs.RO.length.fix(membrane_length)
    if recycle_rate is not None:
        blk0.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix(
            recycle_rate
        )
    blk0.fs.RO.area.fix()
    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if t != mp.n_time_points:
            m.fs.RO.area.setlb(0)
            m.fs.RO.feed_side.velocity[0, 0].setub(None)
            m.fs.RO.feed_side.velocity[0, 0].setlb(None)
    if accumulation_time is not None:
        blk0.fs.dead_volume.accumulation_time[0].fix(accumulation_time)
    mp.total_cycle_time.setub(None)
    mp.overall_recovery.unfix()
    if add_water_recovery_objective:
        mp.min_water_recovery = Objective(
            expr=(mp.overall_recovery * 100 - overal_water_recovery * 100) ** 2
        )
        mp.min_water_recovery.pprint()
        mp.overall_recovery.unfix()
        mp.overall_recovery.setlb(overal_water_recovery + 0.1)
        mp.overall_recovery.setlb(overal_water_recovery - 0.1)
        if mp.find_component("cost_objective") is not None:
            mp.cost_objective.deactivate()
    elif overal_water_recovery is not None:
        mp.overall_recovery.fix(overal_water_recovery)
    if flushing_efficiency is not None:
        blkfs.fs.flushing.flushing_efficiency.fix(flushing_efficiency)
    else:
        blkfs.fs.flushing.flushing_efficiency.unfix()
        blkfs.fs.flushing.flushing_efficiency.setub(0.95)
        blkfs.fs.flushing.flushing_efficiency.setlb(0.25)
    if add_initial_pressure_objective:
        mp.min_initial_pressure = Objective(
            expr=(
                pyunits.convert(blk0.fs.P1.outlet.pressure[0], to_units=pyunits.bar)
                - pyunits.convert(target_pressure, to_units=pyunits.bar)
            )
            ** 2
        )
        if mp.find_component("cost_objective") is not None:
            mp.cost_objective.deactivate()
        mp.min_initial_pressure.pprint()

    print("DOF with MP fixed:", degrees_of_freedom(mp))


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
            to_units=pyunits.L / pyunits.s,
        )()
        if (
            blks.fs.operation_mode == "filtration"
            or blks.fs.operation_mode == "flushing_with_filtration"
        ):
            raw_feed_mass = sum(
                pyunits.convert(
                    blks.fs.RO.feed_side.properties[0, 0].flow_mass_phase_comp[
                        "Liq", comp
                    ],
                    to_units=pyunits.kg / pyunits.s,
                )()
                for comp in blks.fs.properties.component_list
            )
        else:
            raw_feed_mass = 0
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
        if (
            blks.fs.operation_mode == "filtration"
            or blks.fs.operation_mode == "flushing_with_filtration"
        ):

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
    if mp.flushing_points > 1:
        flushing_block = mp.flushing
    else:
        flushing_block = mp.get_active_process_blocks()[-1].fs.flushing
    print(
        f"{'Mean residence time (s):':<{w}s}",
        flushing_block.mean_residence_time.value,
    )
    print(
        f"{'Flushing time (s):':<{w}s}",
        flushing_block.flushing_time.value,
    )
    print(
        f"{'Flushing efficiency:':<{w}s}",
        flushing_block.flushing_efficiency.value,
    )
    print(
        f"{'Pre-flushing conc (kg/m3):':<{w}s}",
        flushing_block.pre_flushing_concentration.value,
    )
    print(
        f"{'Post-flushing conc (kg/m3):':<{w}s}",
        flushing_block.post_flushing_concentration.value,
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
    if mp.get_active_process_blocks()[0].fs.ro_model_with_hold_up:
        print(
            f"{'Total Hold-up Volume:':<{w}s}",
            mp.get_active_process_blocks()[0].fs.RO.feed_side.volume.value
            + mp.get_active_process_blocks()[0].fs.dead_volume.volume[0, "Liq"].value,
        )
        print(
            f"{'RO Hold-up Volume:':<{w}s}",
            mp.get_active_process_blocks()[0].fs.RO.feed_side.volume.value,
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
    if mp.find_component("costing"):
        print(
            f"{'Levelized Cost of Water ($/m3):':<{w}s}",
            value(mp.costing.LCOW),
        )


def validation_configs():
    config = CCROConfiguration()

    config["raw_feed_conc"] = 1.028 * pyunits.g / pyunits.L
    config["raw_feed_flowrate"] = pyunits.convert(
        6.4 * pyunits.gallon / pyunits.ft**2 / pyunits.day * (3 * 37.2) * pyunits.m**2,
        to_units=pyunits.L / pyunits.min,
    )
    config["temperature"] = (273.15 + 20) * pyunits.K
    config["membrane_area"] = 3 * 37.2 * pyunits.m**2
    config["membrane_length"] = 3 * pyunits.m
    config["osmotic_overpressure"] = 1.5 * pyunits.dimensionless
    config["A_comp"] = 1.37778e-11 * pyunits.meter / (pyunits.second * pyunits.Pa)
    config["B_comp"] = 3.2e-08 * pyunits.meter / pyunits.second
    config["channel_height"] = 0.00086 * pyunits.m
    config["spacer_porosity"] = 0.9 * pyunits.dimensionless
    config["dead_volume_to_area_ratio"] = (
        1.1 * 0.00086 * 0.9 * pyunits.m**3 / pyunits.m**2
    )
    config["accumulation_time"] = 10 * pyunits.second
    config["recycle_flowrate"] = 3.28 * pyunits.L / pyunits.s
    config["flushing_efficiency"] = 0.85 * pyunits.dimensionless
    config["pipe_to_module_ratio"] = 1 * pyunits.dimensionless
    config.display()
    return config


def validation_seedling_configs():
    config = CCROConfiguration()

    config["raw_feed_conc"] = 5.831466976 * pyunits.g / pyunits.L
    config["raw_feed_flowrate"] = pyunits.convert(
        15 * pyunits.L / pyunits.m**2 / pyunits.hr * (7.9) * pyunits.m**2,
        to_units=pyunits.L / pyunits.min,
    )
    config["temperature"] = (273.15 + 20) * pyunits.K
    config["membrane_area"] = 7.9 * pyunits.m**2
    config["membrane_length"] = 1 * pyunits.m
    config["osmotic_overpressure"] = 1.5 * pyunits.dimensionless
    config["A_comp"] = (
        5.963600814843386e-12 * pyunits.meter / (pyunits.second * pyunits.Pa)
    )
    config["B_comp"] = 3.0790017613480806e-08 * pyunits.meter / pyunits.second
    config["channel_height"] = 0.0008636000000000001 * pyunits.m
    config["spacer_porosity"] = 0.9 * pyunits.dimensionless
    config["dead_volume_to_area_ratio"] = (
        1.1 * 0.00086 * 0.9 * pyunits.m**3 / pyunits.m**2
    )
    config["pipe_to_module_ratio"] = 1 * pyunits.dimensionless
    config["accumulation_time"] = 10 * pyunits.second
    config["recycle_flowrate"] = 50 * pyunits.L / pyunits.min
    config["flushing_efficiency"] = 0.9 * pyunits.dimensionless
    config.display()
    return config


def check_jac(m, print_extreme_jacobian_values=True):
    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(m, min_scale=1e-8)
    try:
        cond_number = iscale.jacobian_cond(m, jac=jac_scaled) / 1e10
        print("--------------------------")
        print("COND NUMBER:", cond_number)
    except:
        print("Cond number failed")
        cond_number = None
    if print_extreme_jacobian_values:
        print("--------------------------")
        print("Extreme Jacobian entries:")
        extreme_entries = iscale.extreme_jacobian_entries(
            m, jac=jac_scaled, nlp=nlp, zero=1e-20, large=100
        )
        for val, var, con in extreme_entries:
            print(val, var.name, con.name)
        print("--------------------------")
        print("Extreme Jacobian columns:")
        extreme_cols = iscale.extreme_jacobian_columns(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, var in extreme_cols:
            print(val, var.name)
        print("------------------------")
        print("Extreme Jacobian rows:")
        extreme_rows = iscale.extreme_jacobian_rows(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, con in extreme_rows:
            print(val, con.name)
    for var, scale in iscale.badly_scaled_var_generator(m):
        print(
            "Badly scaled variable:",
            var.name,
            var.value,
            iscale.get_scaling_factor(var),
        )
    return cond_number


if __name__ == "__main__":
    from idaes.core.util.model_diagnostics import DiagnosticsToolbox

    cc_config = CCROConfiguration()
    cc_config["raw_feed_conc"] = 5 * pyunits.g / pyunits.L
    cc_config["membrane_area"] = 100 * pyunits.meter**2

    cc_config["accumulation_time"] = 1 * pyunits.second
    # cc_config = validation_seedling_configs()
    periods = 10
    flushing_periods = 5
    mp = create_ccro_multiperiod(
        n_time_points=periods,
        n_flushing_points=flushing_periods,
        include_costing=True,
        cc_configuration=cc_config,
        use_ro_with_hold_up=True,
    )
    setup_optimization(
        mp,
        overall_water_recovery=0.5,
        max_cycle_time_hr=10,
        recycle_flow_bounds=(1, 20),
    )
    from idaes.core.util.model_diagnostics import DiagnosticsToolbox

    blks = list(mp.get_active_process_blocks())
    mp.overall_recovery.fix()  # Unfixed with times fixed should get 1 DOF!
    first_block = blks[0]
    # first_block.fs.dead_volume.accumulation_time[0].fix(1 * pyunits.second)

    # first_block.fs.RO.area.fix()
    # first_block.fs.RO.length.fix()
    # first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix()
    # if flushing_periods > 1:
    #     flushing_block = blks[mp.time_points]
    #     # flushing_block.fs.dead_volume.accumulation_time[0].setub(5 * pyunits.second)
    #     flushing_block.fs.dead_volume.accumulation_time[0].fix(1 * pyunits.second)
    # mp.flushing.display()
    # dg = DiagnosticsToolbox(mp)
    # dg.report_structural_issues()
    # # dg.display_overconstrained_set()
    # flushing_block.fs.conduit_feed.properties[0].flow_vol_phase["Liq"].fix()

    # dg = DiagnosticsToolbox(mp)
    # dg.report_structural_issues()
    # dg.display_overconstrained_set()
    # mp.flushing.deactivate()
    # mp.flushing.flushing_efficiency.fix(0.5)
    # mp.flushing.flushing_time.fix(20 * pyunits.second)
    # fix_optimization_dofs(
    #     mp,
    #     overal_water_recovery=0.5,
    #     add_water_recovery_objective=True,
    #     membrane_area=200 * pyunits.meter**2,
    #     membrane_length=1 * pyunits.meter,
    #     recycle_rate=10 * pyunits.L / pyunits.s,
    #     flushing_efficiency=0.8,
    # )
    # check_jac(m)
    # # assert False
    # mp.find_component(
    #     "blocks[0].process.fs.P2.control_volume.properties_out[0].flow_vol_phase[Liq]"
    # ).fix(0.02)

    results = solve(mp, use_ipoptv2=False)
    print_results_table(mp, w=16)
    for m in mp.get_active_process_blocks():
        if m.fs.find_component("RO") is not None:
            print(
                m.name,
                "inlet RO TDS:",
                m.fs.RO.feed_side.properties[0, 0.1]
                .conc_mass_phase_comp["Liq", "NaCl"]
                .value,
                "outlet RO TDS:",
                m.fs.RO.feed_side.properties[0, 1]
                .conc_mass_phase_comp["Liq", "NaCl"]
                .value,
                m.fs.RO.feed_side.accumulation_time[0].value,
                m.fs.RO.feed_side.volume.value,
            )
    for r in [0.6, 0.7, 0.8, 0.9]:
        mp.overall_recovery.fix(r)  # Unfixed with times fixed should get 1 DOF!
        results = solve(mp, use_ipoptv2=False)
        print_results_table(mp, w=16)
        for m in mp.get_active_process_blocks():
            if m.fs.find_component("RO") is not None:
                print(
                    m.name,
                    "inlet RO TDS:",
                    m.fs.RO.feed_side.properties[0, 0.1]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .value,
                    "outlet RO TDS:",
                    m.fs.RO.feed_side.properties[0, 1]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .value,
                    m.fs.RO.feed_side.accumulation_time[0].value,
                    m.fs.RO.feed_side.volume.value,
                )
