from pyomo.environ import (
    check_optimal_termination,
    units as pyunits,
)
from pyomo.environ import TransformationFactory
from pyomo.network import Arc

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver
from watertap.flowsheets.ccro.multiperiod.model_state_tool import ModelState

from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock


__all__ = [
    "get_variable_pairs",
    "scale_flushing_system",
    "scale_filtration_system",
    "solve",
    "unfix_dof",
    "copy_state_prop_time_period_links",
    "copy_time_period_links",
    "fix_dof_and_initialize",
    "config_op_dict",
    "print_results_table",
    "fix_overall_water_recovery",
    "set_operating_conditions",
    "initialize_system",
]


atmospheric_pressure = 101325 * pyunits.Pa


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

    # set_scaling_factor(
    #     m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"], 1
    # )
    # set_scaling_factor(
    #     m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"], 1
    # )

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

    # set_scaling_factor(
    #     m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"], 1
    # )
    # set_scaling_factor(
    #     m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"], 1
    # )

    calculate_scaling_factors(m)


def solve(model=None, solver=None, tee=False, raise_on_failure=True):
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


def unfix_dof(m, unfix_dead_volume_state=True, op_dict=None, **kwargs):
    """
    Unfix linking variables in MP model
    """

    if unfix_dead_volume_state:
        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()

    if m.fs.configuration == "filtration":

        m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix(
            op_dict["recycle_flowrate"]
        )

        # if unfix_dead_volume_state:
        #     m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        #     m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()

    elif m.fs.configuration == "flushing":

        m.fs.flushing.flushing_time.unfix()
        m.fs.flushing.mean_residence_time.unfix()
        m.fs.flushing.flushing_efficiency.fix(op_dict["flushing_efficiency"])
        m.fs.dead_volume.dead_volume.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.flushing.pre_flushing_concentration.unfix()
        m.fs.flushing.post_flushing_concentration.unfix()

        # if unfix_dead_volume_state:
        #     m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        #     m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()


def copy_state_prop_time_period_links(m_old, m_new):
    copy_state(m_old, m_new)
    copy_time_period_links(m_old, m_new)


def copy_time_period_links(m_old, m_new):
    """
    Copy linking variables between time periods.
    """
    m_new.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
        m_old.fs.dead_volume.dead_volume.properties_out[0]
        .mass_frac_phase_comp["Liq", "NaCl"]
        .value
    )
    m_new.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
        m_old.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"].value
    )


def copy_state(old_model, new_model):
    model_state = ModelState()
    model_state.get_model_state(old_model)
    model_state.set_model_state(new_model)


def fix_dof_and_initialize(m, op_dict=None, **kwargs):
    """
    Fix DOF for MP model and initialize steady-state models.
    """

    set_operating_conditions(m=m, op_dict=op_dict, **kwargs)
    initialize_system(m=m, **kwargs)


def copy_inlet_state_for_mixer(m):
    for idx, obj in m.fs.M1.inlet_2.flow_mass_phase_comp.items():
        obj.value = m.fs.M1.inlet_1.flow_mass_phase_comp[idx].value * 1


def set_operating_conditions(m, op_dict=None, **kwargs):

    if m.fs.configuration == "filtration":

        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(atmospheric_pressure)
        m.fs.raw_feed.properties[0].temperature.fix(op_dict["temperature"])

        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(op_dict["p1_eff"])
        m.fs.P1.control_volume.properties_out[0].pressure.fix(
            op_dict["p1_pressure_start"]
        )

        # Pump 2 operating conditions
        m.fs.P2.efficiency_pump.fix(op_dict["p2_eff"])

        # Set RO configuration parameters
        m.fs.RO.A_comp.fix(op_dict["A_comp"])
        m.fs.RO.B_comp.fix(op_dict["B_comp"])
        m.fs.RO.area.fix(op_dict["membrane_area"])
        m.fs.RO.length.fix(op_dict["membrane_length"])
        m.fs.RO.feed_side.channel_height.fix(op_dict["channel_height"])
        m.fs.RO.feed_side.spacer_porosity.fix(op_dict["spacer_porosity"])

        m.fs.RO.permeate.pressure[0].fix(atmospheric_pressure)

        # m.fs.RO.feed_side.K.setlb(1e-6)
        m.fs.RO.feed_side.friction_factor_darcy.setub(200)
        # m.fs.RO.flux_mass_phase_comp.setub(1)
        # m.fs.RO.feed_side.cp_modulus.setub(50)
        # m.fs.RO.feed_side.cp_modulus.setlb(0.1)
        m.fs.RO.deltaP.setlb(None)

        # for e in m.fs.RO.permeate_side:
        #     if e[-1] != 0:
        #         m.fs.RO.permeate_side[e].pressure_osm_phase["Liq"].setlb(200)
        #         m.fs.RO.permeate_side[e].molality_phase_comp["Liq", "NaCl"].setlb(
        #             1e-8
        #         )

        # Single pass water recovery constraint
        # m.fs.single_pass_water_recovery_constraint = Constraint(
        #     expr=m.fs.RO.flow_vol_re
        #     == m.fs.single_pass_water_recovery
        # )

        # Dead Volume operating conditions
        # Fixed volume
        m.fs.dead_volume.volume.fix(op_dict["dead_volume"])
        m.fs.dead_volume.delta_state.volume.fix(op_dict["dead_volume"])

        m.fs.dead_volume.accumulation_time.fix(op_dict["accumulation_time"])

        # Fixing the flow rate of the dead volume delta state
        # Using the feed to calculate the mass fraction and density

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
            op_dict["recycle_flowrate"]
        )
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
            op_dict["recycle_conc_start"]
        )
        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].unfix()
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()

        # I found fixing mass fraction and density is easiest way to get initial state
        # we will also use these as connection points between current and future state.

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m.fs.raw_feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value
        )
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
        )

        # Reassign the raw feed flowrate and concentration
        m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            op_dict["raw_feed_flow_mass_water"]
        )
        m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            op_dict["raw_feed_flow_mass_salt"]
        )

        solver = get_solver()
        solver.solve(m.fs.raw_feed)
        scale_filtration_system(m)

    elif m.fs.configuration == "flushing":

        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(op_dict["p1_eff"])
        m.fs.P1.control_volume.properties_out[0].pressure.fix(
            op_dict["p1_pressure_start"]
        )

        # Pump 2 operating conditions - Add only for costing function. No work is done by this pump
        m.fs.P2.efficiency_pump.fix(op_dict["p2_eff"])

        # Concentration of the flushing water is the raw feed concentration
        m.fs.flushing.raw_feed_concentration.fix(op_dict["raw_feed_conc"])
        # m.fs.flushing.flushing_efficiency.fix(0.5)

        # m.fs.flushing.mean_residence_time.fix(pyunits.convert(op_dict["dead_volume"]/op_dict["flushing_flowrate"], to_units=pyunits.s))

        # m.fs.flushing.mean_residence_time_constr = Constraint(expr=pyunits.convert(m.fs.dead_volume.volume[0, "Liq"] / m.fs.raw_feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.s) == m.fs.flushing.mean_residence_time)
        
        # m.fs.flushing.number_tanks_in_series.set_value(3)
        # m.fs.flushing.accumulator_volume.set_value(dead_volume)
        # m.fs.flushing.flushing_flow_rate.set_value(flushing_flowrate)

        # Calculating values for the dead volume and delta state
        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(atmospheric_pressure)
        m.fs.raw_feed.properties[0].temperature.fix(op_dict["temperature"])

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

        # Fixed volume

        m.fs.dead_volume.volume.fix(op_dict["dead_volume"])
        m.fs.dead_volume.delta_state.volume.fix(op_dict["dead_volume"])

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix()

        # Constraints
        @m.fs.Constraint()
        def dead_volume_residence_time_constraint(m):
            return (
                m.flushing.mean_residence_time
                == (pyunits.convert(
                   m.dead_volume.volume[0, "Liq"]
                    / m.raw_feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.s)
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

        constraint_scaling_transform(m.fs.post_flushing_conc_constraint, 1e-1)
        constraint_scaling_transform(m.fs.pre_flushing_conc_constraint, 1e-1)
        set_scaling_factor(m.fs.flushing.pre_flushing_concentration, 1e-1)
        set_scaling_factor(m.fs.flushing.post_flushing_concentration, 1e-1)
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].unfix()
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()

        # Reassign the raw feed flowrate and concentration
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
            op_dict["flushing_flowrate"]
        )
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
            op_dict["raw_feed_conc"]
        )

        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        scale_flushing_system(m)

    return m


def initialize_system(m, **kwargs):
    """
    Initialize the model by fixing the values of certain variables.
    """

    if m.fs.configuration == "filtration":

        # m.fs.raw_feed.initialize()
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

    if m.fs.configuration == "flushing":

        propagate_state(m.fs.raw_feed_to_P1)
        m.fs.P1.initialize()

        propagate_state(m.fs.P1_to_dead_volume)
        m.fs.dead_volume.initialize()

        propagate_state(m.fs.dead_volume_to_P2)
        m.fs.P2.initialize()
        m.fs.P2.outlet.pressure[0].fix(atmospheric_pressure)

        calculate_variable_from_constraint(
            m.fs.flushing.pre_flushing_concentration,
            m.fs.pre_flushing_conc_constraint,
        )
        m.fs.flushing.pre_flushing_concentration.display()
        m.fs.flushing.initialize()
    return m


def config_op_dict(op_dict):
    """
    Configure operating dictionary to have mass flowrates instead of volumetric flowrates
    """
    op_dict2 = op_dict

    op_dict2["rho"] = op_dict["rho"] * pyunits.kg / pyunits.m**3
    op_dict2["raw_feed_conc"] = op_dict["raw_feed_conc"] * pyunits.gram / pyunits.liter
    op_dict2["temperature"] = op_dict["temperature"] * pyunits.degK
    op_dict2["raw_feed_flowrate"] = (
        op_dict["raw_feed_flowrate"] * pyunits.liter / pyunits.minute
    )
    op_dict2["p1_eff"] = op_dict["p1_eff"] * pyunits.dimensionless
    op_dict2["p2_eff"] = op_dict["p2_eff"] * pyunits.dimensionless
    op_dict2["p1_pressure_start"] = op_dict["p1_pressure_start"] * pyunits.psi
    op_dict2["p2_pressure_start"] = op_dict["p2_pressure_start"] * pyunits.psi
    op_dict2["A_comp"] = op_dict["A_comp"] * pyunits.m / (pyunits.s * pyunits.Pa)
    op_dict2["B_comp"] = op_dict["B_comp"] * pyunits.m / pyunits.s
    op_dict2["recycle_flowrate"] = (
        op_dict["recycle_flowrate"] * pyunits.liter / pyunits.minute
    )
    op_dict2["flushing_flowrate"] = (
        op_dict["raw_feed_flowrate"] + op_dict["recycle_flowrate"]
    )
    op_dict["dead_volume"] = op_dict["dead_volume"] * pyunits.m**3
    op_dict["accumulation_time"] = op_dict["accumulation_time"] * pyunits.second

    op_dict2["raw_feed_flow_mass_water"] = (
        op_dict2["raw_feed_flowrate"] * op_dict2["rho"]
    )

    op_dict2["raw_feed_flow_mass_salt"] = (
        op_dict2["raw_feed_flowrate"] * op_dict2["raw_feed_conc"]
    )

    return op_dict2


def print_results_table(mp, w=15):
    """
    Print multiperiod CCRO results in a tabular format in the terminal.
    w = column width
    """
    n = 12  # number of columns
    title = "CCRO MULTIPERIOD RESULTS"
    side = int(((n * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    # print("\n" + "=" * 120)
    # print("CCRO MULTIPERIOD RESULTS")
    print(f"{'-' * (n * w)}")

    # Header
    print(
        f"{'Period':<{w}s}{'Acc Time':<{w}s}{'Raw Feed':<{w}s}{'Permeate':<{w}s}{'SP Recovery':<{w}s}{'P2':<{w}s}{'RO In':<{w}s}{'RO In':<{w}s}{'Dead Vol In':<{w}s}{'Dead Vol In':<{w}s}{'Delta State':<{w}s}{'Dead Vol':<{w}s}"
    )
    print(
        f"{'':<{w}s}{'(s)':<{w}s}{'(L/min)':<{w}s}{'(L/min)':<{w}s}{'(%)':<{w}s}{'(L/min)':<{w}s}{'(Pa)':<{w}s}{'(Psi)':<{w}s}{'(L/min)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}{'(kg/m³)':<{w}s}"
    )
    print(f"{'-' * (n * w)}")

    # Data rows
    for t, blks in enumerate(mp.get_active_process_blocks(), 1):
        # if t == len(blks):
        if blks.fs.configuration == "flushing":
            accumulation_time = blks.fs.flushing.flushing_time.value
        else:
            accumulation_time = blks.fs.dead_volume.accumulation_time[0].value

        raw_feed = pyunits.convert(blks.fs.raw_feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.min)()
        
        if blks.fs.find_component("product") is not None:
            permeate = pyunits.convert(blks.fs.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.min)()
        else:
            permeate = 0
        if blks.fs.configuration == "filtration":
            # mixer_out = value(
            #     pyunits.convert(
            #         blks.fs.M1.mixed_state[0].flow_vol_phase["Liq"],
            #         to_units=pyunits.L / pyunits.min,
            #     )
            # )
            p2_out = value(
                pyunits.convert(
                    blks.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.L / pyunits.min,
                )
            )
            ro_pressure = blks.fs.M1.mixed_state[0].pressure.value
            ro_pressure_psi = value(
                pyunits.convert(
                    blks.fs.M1.mixed_state[0].pressure, to_units=pyunits.psi
                )
            )
            sp_recovery = blks.fs.RO.recovery_vol_phase[0.0, "Liq"].value

        if blks.fs.configuration == "flushing":
            mixer_out = value(
                pyunits.convert(
                    blks.fs.raw_feed.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.L / pyunits.min,
                )
            )
            sp_recovery = 0
            ro_pressure = 0

        dead_vol_in = (
            pyunits.convert(blks.fs.dead_volume.dead_volume.properties_in[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.min)()
        )
        dead_vol_in_conc = (
            blks.fs.dead_volume.dead_volume.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )
        # dead_vol_out = value(
        #     pyunits.convert(
        #         blks.fs.dead_volume.dead_volume.properties_out[0].flow_vol_phase["Liq"],
        #         to_units=pyunits.L / pyunits.min,
        #     )
        # )
        delta_state_conc = blks.fs.dead_volume.delta_state.conc_mass_phase_comp["Liq", "NaCl"]()
        dead_vol_out_conc = (
            blks.fs.dead_volume.dead_volume.properties_out[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
            .value
        )

        print(
            f"{t:<{w}d}{accumulation_time:<{w}.2f}{raw_feed:<{w}.6f}{permeate:<{w}.6f}{sp_recovery:<{w}.6f}{p2_out:<{w}.6f}{ro_pressure:<{w}.2f}{ro_pressure_psi:<{w}.2f}{dead_vol_in:<{w}.6f}{dead_vol_in_conc:<{w}.6f}{delta_state_conc:<{w}.6f}{dead_vol_out_conc:<{w}.6f}"
        )

    print(f"{'=' * (n * w)}")
    w = w * 2

    print(
        f"{'Total cycle time (s):':<{w}s}",
        mp.total_cycle_time.value,
    )
    print(
         f"{'Mean residence time (s):':<{w}s}",
        mp.get_active_process_blocks()[
            -1
        ].fs.flushing.mean_residence_time.value,
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
