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
from watertap.unit_models.pseudo_steady_state.reverse_osmosis_1D_with_holdup import (
    ReverseOsmosis1DwithHoldUp,
)
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

from watertap.unit_models.pseudo_steady_state.flushing import Flushing
import watertap.flowsheets.ccro.utils.ipoptv2 as ipt2
import pandas as pd


def load_validation_data_into_model(
    mp=None,
    file_path="validation_data/sine_500-900psi_60s_period.csv",
    data_point_skips=20,
):
    """
    Load validation data from a CSV file.
    Args:
        file_path: Path to the CSV file.
    Returns:
        DataFrame containing the validation data.
    """
    data = pd.read_csv(file_path)
    if mp is not None:
        for t, m in enumerate(mp.get_active_process_blocks()):
            dp = (t + 1) * data_point_skips

            run_time = (
                data["Runtime (min)"].to_numpy()[dp]
                - data["Runtime (min)"].to_numpy()[dp - data_point_skips]
            ) * pyunits.min
            print(
                f"Loading data point {data['Runtime (min)'].to_numpy()[dp]*60} into time period {t}"
            )
            m.fs.time_point = Var(
                initialize=run_time, domain=NonNegativeReals, units=pyunits.min
            )
            m.fs.time_point.fix(data["Runtime (min)"].to_numpy()[dp] * pyunits.min)
            flow_rate = (
                data["Feed Flowrate (L/min)"].to_numpy()[dp] * pyunits.L / pyunits.min
            )
            pressure = data["Feed Pressure (psi)"].to_numpy()[dp] * pyunits.psi
            m.fs.P1.outlet.pressure[0].fix(pressure)
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(flow_rate)

            if m.fs.ro_model_with_hold_up:
                m.fs.RO.feed_side.accumulation_time.fix(run_time)

    return data


def create_flush_setup(
    mp,
    starting_c=35,
    flushing_c=5,
    flush_flow=10 * pyunits.L / pyunits.s,
    feed_flow=1 * pyunits.L / pyunits.s,
    # time_step=1 * pyunits.s,
):
    """
    Create flushing setup for multiperiod CCRO model
    """

    dead_volume = mp.get_active_process_blocks()[0].fs.RO.feed_side.volume
    time_step = (
        value(
            pyunits.convert(
                dead_volume * 2 / (feed_flow + flush_flow), to_units=pyunits.seconds
            )
        )
        * pyunits.seconds
    )
    time_step = time_step / len(mp.get_active_process_blocks())
    print("Dead Volume:", value(dead_volume))
    print(f"Calculated time step for flushing: {time_step}")
    for t, m in enumerate(mp.get_active_process_blocks()):
        if t == 0:
            m.fs.P1.outlet.pressure[0].unfix()
            m.fs.RO.mixed_permeate[0].flow_vol_phase["Liq"].fix(feed_flow)
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                feed_flow + flush_flow
            )
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                starting_c
            )

            if m.fs.ro_model_with_hold_up:
                m.fs.RO.feed_side.accumulation_time.fix(time_step)
        else:
            m.fs.P1.outlet.pressure[0].unfix()
            m.fs.RO.mixed_permeate[0].flow_vol_phase["Liq"].fix(feed_flow)
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                feed_flow + flush_flow
            )
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                flushing_c
            )

            if m.fs.ro_model_with_hold_up:
                m.fs.RO.feed_side.accumulation_time.fix(time_step)


def create_ccro_multiperiod(
    n_time_points=10,
    cc_configuration=None,
    use_ro_with_hold_up=False,
    add_flushing=False,
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
        process_model_func=build_psro_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        solver=watertap_solver,
        outlvl=logging.WARNING,
    )

    mp.n_time_points = n_time_points

    flowsheet_options = {
        t: {
            "time_blk": t,
            "use_ro_with_hold_up": use_ro_with_hold_up,
        }
        for t in range(n_time_points)
    }

    # Build instances of the process model for each time period
    mp.build_multi_period_model(model_data_kwargs=flowsheet_options)

    for t, m in enumerate(mp.get_active_process_blocks(), 1):
        if t == 1:
            fix_dof_and_initialize(m, cc_configuration=cc_configuration)
            # unfix_dof(m)
            old_m = m
        # Last time period is flushing
        base_links = []

        cc_utils.copy_state(old_m, m)
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
            cc_utils.copy_time_period_links(old_m, m, base_links)
        print(f"Period {t} DOF:", degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0

        results = solve(model=m, tee=True)
        assert_optimal_termination(results)
        if t != 1:
            unfix_dof(m)
        old_m = m
        if t == 1:
            unfix_dof(m)
            create_stabilization_constraints(m)

    calculate_scaling_factors(mp)

    print("Multi-period DOF:", degrees_of_freedom(mp))
    return mp


def add_flushing_unit(mp, cc_configuration, start_period=1, end_period=None):

    m_start = mp.get_active_process_blocks()[start_period]

    mp.flushing = Flushing()
    if end_period is None:
        blks = mp.get_active_process_blocks()[start_period:]

        m_end = mp.get_active_process_blocks()[-1]
    else:
        blks = mp.get_active_process_blocks()[start_period:end_period]

        m_end = mp.get_active_process_blocks()[end_period]
    print(blks)
    print(m_start.name)
    print(m_end.name)
    print(list(range(start_period, len(blks))))

    @mp.Constraint(list(range(start_period, len(blks))))
    def eq_flushing_time_constraint(m, i):
        return blks[0].fs.RO.feed_side.accumulation_time[0] == (
            blks[i].fs.RO.feed_side.accumulation_time[0]
        )

    for blk in blks:
        print(blk.name)
        blk.fs.RO.feed_side.accumulation_time.unfix()
        blk.fs.RO.feed_side.accumulation_time[0].setlb(1e-8)

    # Constraints
    @mp.Constraint()
    def dead_volume_residence_time_constraint(m):
        return mp.flushing.mean_residence_time == (
            pyunits.convert(
                m_start.fs.RO.feed_side.volume
                / m_start.fs.raw_feed.properties[0].flow_vol_phase["Liq"],
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

    iscale.constraint_scaling_transform(mp.pre_flushing_conc_constraint, 1)
    # iscale.constraint_scaling_transform(mp.post_flushing_conc_constraint, 1)
    iscale.set_scaling_factor(mp.flushing.pre_flushing_concentration, 1)
    iscale.set_scaling_factor(mp.flushing.post_flushing_concentration, 1)
    calculate_variable_from_constraint(
        mp.flushing.pre_flushing_concentration,
        mp.pre_flushing_conc_constraint,
    )
    calculate_variable_from_constraint(
        mp.flushing.mean_residence_time,
        mp.dead_volume_residence_time_constraint,
    )
    mp.flushing.mean_residence_time.fix()

    mp.flushing.flushing_time.unfix()
    mp.flushing.flushing_time = value(mp.flushing.mean_residence_time)
    mp.flushing.pre_flushing_concentration.fix()
    mp.flushing.post_flushing_concentration.fix(
        value(
            m_end.fs.RO.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
    )

    mp.flushing.flushing_efficiency.unfix()
    mp.flushing.flushing_feed_concentration.fix(
        5 * pyunits.g / pyunits.L  # value(mp.flushing.pre_flushing_concentration) / 2
    )
    mp.flushing.display()
    mp.flushing.initialize()
    mp.flushing.pre_flushing_concentration.unfix()
    mp.flushing.flushing_time.unfix()
    mp.flushing.flushing_time.setub(1000)
    mp.flushing.flushing_time.setlb(1e-8)
    mp.flushing.mean_residence_time.unfix()
    mp.flushing.post_flushing_concentration.fix(6 * pyunits.g / pyunits.L)
    # mp.flushing.flushing_efficiency.fix(0.9)


def fix_dof_and_initialize(m, cc_configuration=None, **kwargs):
    """
    Fix DOF for MP model and initialize steady-state models.
    """
    # Feed block operating conditions
    m.fs.raw_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
    m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
    m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
    m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
    m.fs.raw_feed.properties[0].temperature.fix(cc_configuration["temperature"])
    m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
        cc_configuration["raw_feed_flowrate"]
    )
    m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
        cc_configuration["raw_feed_conc"]
    )

    solver = get_solver()
    solver.solve(m.fs.raw_feed)
    flow_mass = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    flow_nacl = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value
    print(f"Flow H2O: {flow_mass}, Flow NaCl: {flow_nacl}")
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_mass, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_nacl, index=("Liq", "NaCl")
    )

    # Dead Volume operating conditions
    # Fix volume
    if cc_configuration["dead_volume"] == "base_on_dead_volume_to_area_ratio":
        dead_volume = (
            cc_configuration["dead_volume_to_area_ratio"]
            * cc_configuration["membrane_area"]
        )
        print(f"Dead Volume set to {dead_volume}")
    else:
        dead_volume = cc_configuration["dead_volume"] = cc_configuration["dead_volume"]
    if m.fs.ro_model_with_hold_up:
        dead_volume = dead_volume * cc_configuration["pipe_to_module_ratio"]

    # Fixing the flow rate of the dead volume delta state to the feed to calculate the mass fraction and density
    # I found fixing mass fraction and density is easiest way to get initial state
    # we will also use these as connection points between current and future state.

    if m.fs.ro_model_with_hold_up:
        m.fs.RO.feed_side.volume.fix(
            cc_configuration["dead_volume_to_area_ratio"]
            * cc_configuration["membrane_area"]
        )
        m.fs.RO.feed_side.accumulation_time.fix(cc_configuration["accumulation_time"])
        idx = m.fs.RO.difference_elements
        for i in idx:
            m.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                0, i, "Liq", "NaCl"
            ].fix(m.fs.raw_feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value)
            m.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"].fix(
                m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
            )
        solver.solve(m.fs.RO.feed_side.delta_state)
    # Pump 1 operating conditions
    m.fs.P1.efficiency_pump.fix(cc_configuration["p1_eff"])

    cc_utils.set_pump_operating_pressure(
        unit=m.fs.P1,
        configuration_options=cc_configuration,
        osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
    )
    # Pump 2 operating conditions

    # Set RO configuration parameters
    m.fs.RO.A_comp.fix(cc_configuration["A_comp"])
    m.fs.RO.B_comp.fix(cc_configuration["B_comp"])
    m.fs.RO.area.fix(cc_configuration["membrane_area"])
    print(
        cc_configuration["membrane_length"],
    )
    if cc_configuration["membrane_length"] == "auto":
        m.fs.RO.length.unfix()
        m.fs.RO.feed_side.velocity[0, 0].fix(0.15)
    else:
        m.fs.RO.length.fix(cc_configuration["membrane_length"])
        m.fs.RO.feed_side.velocity[0, 0].unfix()
    m.fs.RO.width.unfix()
    m.fs.RO.feed_side.channel_height.fix(cc_configuration["channel_height"])
    m.fs.RO.feed_side.spacer_porosity.fix(cc_configuration["spacer_porosity"])

    m.fs.RO.permeate.pressure[0].fix(101325 * pyunits.Pa)

    m.fs.RO.feed_side.friction_factor_darcy.setub(200)
    # m.fs.RO.flux_mass_phase_comp.setub(1)
    # m.fs.RO.feed_side.cp_modulus.setub(50)
    # m.fs.RO.feed_side.cp_modulus.setlb(0.1)
    m.fs.RO.deltaP.setlb(None)

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)

    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    iscale.calculate_scaling_factors(m)
    initialize_system(m=m, **kwargs)


def unfix_dof(
    m,
):
    if m.fs.ro_model_with_hold_up:
        idx = m.fs.RO.difference_elements
        for i in idx:
            m.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                0, i, "Liq", "NaCl"
            ].unfix()
            m.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"].unfix()


def create_stabilization_constraints(m):
    if m.fs.ro_model_with_hold_up:

        @m.fs.RO.feed_side.Constraint(m.fs.RO.difference_elements)
        def frac_stabilization_constraint(b, i):
            return (
                b.delta_state.node_mass_frac_phase_comp[0, i, "Liq", "NaCl"]
                == b.properties[0, i].mass_frac_phase_comp["Liq", "NaCl"]
            )

        @m.fs.RO.feed_side.Constraint(m.fs.RO.difference_elements)
        def dense_stabilization_constraint(b, i):
            return (
                b.delta_state.node_dens_mass_phase[0, i, "Liq"]
                == b.properties[0, i].dens_mass_phase["Liq"]
            )


def build_psro_system(time_blk=None, use_ro_with_hold_up=True):
    """
    Build the CCRO steady state flowsheet
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_model_with_hold_up = use_ro_with_hold_up
    m.fs.properties = NaClParameterBlock()

    # Low concentration feed to the system
    m.fs.raw_feed = Feed(property_package=m.fs.properties)

    # Touch relevant parameters
    m.fs.raw_feed.properties[0].flow_vol_phase
    m.fs.raw_feed.properties[0].conc_mass_phase_comp
    # Feed pump
    m.fs.P1 = Pump(property_package=m.fs.properties)
    # Recirculation pump
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

    # Add connections
    m.fs.raw_feed_to_P1 = Arc(source=m.fs.raw_feed.outlet, destination=m.fs.P1.inlet)

    m.fs.P1_to_RO = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def initialize_system(m, **kwargs):
    """
    Initialize the model by fixing the values of certain variables.
    """

    propagate_state(m.fs.raw_feed_to_P1)

    # m.fs.P2.outlet.pressure[0].fix(m.fs.P1.outlet.pressure[0].value)
    m.fs.P1.initialize()

    propagate_state(m.fs.P1_to_RO)
    m.fs.RO.initialize()
    return m


def get_variable_pairs(t1, t2):
    """
    Get variable pairs for connecting two time periods.
    1. dead_volume mass fraction to delta_state
    2. dead_volume density to delta_state
    """
    pair_list = []
    if t1.fs.ro_model_with_hold_up:
        idx = t1.fs.RO.difference_elements
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
                    t2.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"],
                    t1.fs.RO.feed_side.properties[0, i].dens_mass_phase["Liq"],
                ),
            )
    return pair_list


def solve(
    model=None,
    solver=None,
    tee=True,
    raise_on_failure=True,
    use_ipoptv2=False,
    **kwargs,
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
    if model.find_component("flushing") is not None:
        model.flushing.display()
    if hasattr(model, "get_active_process_blocks"):
        for blk in model.get_active_process_blocks():
            print(
                blk.fs.RO.feed_side.properties[0, 1]
                .conc_mass_phase_comp["Liq", "NaCl"]
                .value,
                blk.fs.RO.feed_side.accumulation_time[0].value,
            )
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


def validation_configs():
    config = CCROConfiguration()

    config["raw_feed_conc"] = 34.8 * pyunits.g / pyunits.L
    config["raw_feed_flowrate"] = 21.28 * pyunits.L / pyunits.min
    config["temperature"] = (273.15 + 20) * pyunits.K
    config["membrane_area"] = 7.4 * pyunits.m**2
    config["membrane_length"] = 1 * pyunits.m
    config["osmotic_overpressure"] = 1.5 * pyunits.dimensionless
    config["A_comp"] = 4.4e-12 * pyunits.meter / (pyunits.second * pyunits.Pa)
    config["B_comp"] = 3.08e-08 * pyunits.meter / pyunits.second
    config["channel_height"] = 0.0007112 * pyunits.m
    config["spacer_porosity"] = 0.89 * pyunits.dimensionless
    config["dead_volume_to_area_ratio"] = (
        1.1 * 0.0007112 * 0.89 * pyunits.m**3 / pyunits.m**2
    )
    config["accumulation_time"] = 5 * pyunits.second
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

    cc_config = CCROConfiguration()  # validation_configs()
    cc_config["raw_feed_conc"] = 35 * pyunits.g / pyunits.L
    # cc_config["membrane_area"] = 100 * pyunits.m**2
    cc_config["raw_feed_flowrate"] = 6 * pyunits.L / pyunits.s  # validation_configs()
    # cc_config["membrane_length"] = 10 * pyunits.m
    cc_config.display()
    mp = create_ccro_multiperiod(
        n_time_points=50,
        cc_configuration=cc_config,
        use_ro_with_hold_up=True,
    )
    create_flush_setup(mp)
    # load_validation_data_into_model(
    #     mp,
    #     file_path="validation_data/sine_700-900psi_60s_period.csv",
    #     data_point_skips=40,
    # )
    solve(mp, use_ipoptv2=False)

    add_flushing_unit(mp, cc_configuration=cc_config)
    solve(mp, use_ipoptv2=False)
