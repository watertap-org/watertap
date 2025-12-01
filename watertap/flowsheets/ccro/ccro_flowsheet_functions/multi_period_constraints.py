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

from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.scaling as iscale


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

        iscale.set_scaling_factor(m.fs.dead_volume.accumulation_time, 1e-2)

    # Equal accumulation time across all filtration periods
    @mp.Constraint(list(range(1, mp.n_time_points - 1)))
    def accumulation_time_cons(mp, t):
        blks = list(mp.get_active_process_blocks())
        return blks[t].fs.dead_volume.accumulation_time[0] == (
            blks[0].fs.dead_volume.accumulation_time[0]
        )
