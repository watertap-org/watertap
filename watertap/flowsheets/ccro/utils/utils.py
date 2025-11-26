from pyomo.environ import (
    check_optimal_termination,
    units as pyunits,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    Block,
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
import idaes.core.util.scaling as iscale
from watertap.flowsheets.ccro.multiperiod.model_state_tool import ModelState
from idaes.core import FlowsheetBlock, UnitModelCostingBlock


__all__ = [
    "scale_flushing_system",
    "scale_filtration_system",
    "solve",
    "config_op_dict",
    "print_results_table",
    "set_operating_conditions",
    "initialize_system",
    "get_variable_pairs",  
    "unfix_dof",
    "fix_dof_and_initialize",
    "add_multiperiod_variables",
    "add_multiperiod_constraints",
    "scale_multiperiod_model",
    "copy_state_prop_time_period_links",
    "copy_time_period_links",
    "copy_state",
    "copy_inlet_state_for_mixer",
    "register_costed_unit",
    "add_ccro_costing",
]

atmospheric_pressure = 101325 * pyunits.Pa

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


def composition_calculator(m, state_block = None):
    """
    This function returns the composition of the stream based on inputs:
    1. Volumetric flowrate
    2. Concentration
    3. Temperature
    4. Pressure
    """
    # Create minimal temporary model
    temp_block = Block()
    m.add_component("temp_state_block", temp_block)
    temp_block.props = m.fs.properties.build_state_block([0], defined_state=True)
    # Set state
    pressure = 101325 * pyunits.Pa
    temperature = 298.3 * pyunits.K
    flow_vol = 1.5 * pyunits.m**3/pyunits.s
    concentration = 0.5 * pyunits.kg/pyunits.m**3

    temp_block.props[0].pressure.fix(pressure)
    temp_block.props[0].temperature.fix(temperature)
    temp_block.props[0].flow_vol_phase["Liq"].fix(flow_vol)
    temp_block.props[0].conc_mass_phase_comp["Liq", "NaCl"].fix(concentration)
    
    # Solve
    solver = get_solver()
    solver.solve(temp_block)

    if state_block is not None:
        state_block = temp_block
        # Delete temp block
        m.del_component("temp_state_block")


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

        m.fs.RO.feed_side.friction_factor_darcy.setub(200)
        # m.fs.RO.flux_mass_phase_comp.setub(1)
        # m.fs.RO.feed_side.cp_modulus.setub(50)
        # m.fs.RO.feed_side.cp_modulus.setlb(0.1)
        m.fs.RO.deltaP.setlb(None)

        # Dead Volume operating conditions
        # Fix volume
        m.fs.dead_volume.volume.fix(op_dict["dead_volume"])
        m.fs.dead_volume.delta_state.volume.fix(op_dict["dead_volume"])

        m.fs.dead_volume.accumulation_time.fix(op_dict["accumulation_time"])

        # Fixing the flow rate of the dead volume delta state
        # Using the feed to calculate the mass fraction and density --> Updated to use a helper function

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
        
        # Reassign the raw feed flowrate and concentration
        m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            op_dict["raw_feed_flow_mass_water"]
        )
        m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            op_dict["raw_feed_flow_mass_salt"]
        )
        # m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
        #     op_dict["recycle_flowrate"]
        # )
        # m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
        #     op_dict["recycle_conc_start"]
        # )
        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        # m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].unfix()
        # m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()

        # I found fixing mass fraction and density is easiest way to get initial state
        # we will also use these as connection points between current and future state.

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m.fs.raw_feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value
        )
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
        )

        # Reassign the raw feed flowrate and concentration
        # m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        #     op_dict["raw_feed_flow_mass_water"]
        # )
        # m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        #     op_dict["raw_feed_flow_mass_salt"]
        # )

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
            # p2_out = value(
            #     pyunits.convert(
            #         blks.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"],
            #         to_units=pyunits.L / pyunits.min,
            #     )
            # )
            p2_out = value(
                pyunits.convert(
                    blks.fs.P2.work_mechanical[0], to_units=pyunits.kW
                )
            )
            p1_out = value(
                pyunits.convert(
                    blks.fs.P1.work_mechanical[0], to_units=pyunits.kW
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

    elif m.fs.configuration == "flushing":
        m.fs.flushing.flushing_efficiency.fix(op_dict["flushing_efficiency"])
        m.fs.flushing.flushing_time.unfix()
        m.fs.flushing.mean_residence_time.unfix()
        m.fs.dead_volume.dead_volume.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.flushing.pre_flushing_concentration.unfix()
        m.fs.flushing.post_flushing_concentration.unfix()


def fix_dof_and_initialize(m, op_dict=None, **kwargs):
    """
    Fix DOF for MP model and initialize steady-state models.
    """

    set_operating_conditions(m=m, op_dict=op_dict, **kwargs)
    initialize_system(m=m, **kwargs)


def add_multiperiod_variables(mp):
    """
    Add variables to the multiperiod model.
    """

    bs = list(mp.get_active_process_blocks())
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

    accumulation_time = mp.op_dict["accumulation_time"]

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
                accumulation_time * (len(blks) - 1) + blks[-1].fs.flushing.flushing_time
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
    pass


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


def copy_inlet_state_for_mixer(m):
    for idx, obj in m.fs.M1.inlet_2.flow_mass_phase_comp.items():
        obj.value = m.fs.M1.inlet_1.flow_mass_phase_comp[idx].value * 1


def register_costed_unit(
    mp, unit, costing_method_arguments={}, register_electricity_flow_only=False
):
    if register_electricity_flow_only:
        lb = unit.work_mechanical[0.0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        unit.work_mechanical.setlb(0)
        mp.costing.cost_flow(
            pyunits.convert(unit.work_mechanical[0.0], to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        unit.work_mechanical.setlb(lb)
    else:
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mp.costing,
            costing_method_arguments=costing_method_arguments,
        )


def add_ccro_costing(m=None, mp=None):
    """
    Add costing blocks to steady-state model.
    """

    m.fs.costing = WaterTAPCosting()
    costing_method_arguments = dict(mp=mp)

    m.fs.RO.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments=costing_method_arguments,
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(mp.avg_product_flow_rate)
    m.fs.costing.add_specific_energy_consumption(mp.avg_product_flow_rate, name="SEC")
