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


def add_multiperiod_variables(mp, cc_configuration=None):
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
    mp.dead_volume_to_area_ratio.fix(cc_configuration["dead_volume_to_area_ratio"])
    mp.pipe_to_module_ratio = Var(
        initialize=0.2,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Global dead volume over all time periods",
    )
    iscale.set_scaling_factor(mp.pipe_to_module_ratio, 1)
    mp.pipe_to_module_ratio.fix(cc_configuration["pipe_to_module_ratio"])
    mp.dead_volume_to_area_multiplier.fix(1)
    # Permeate and feed
    mp.overall_recovery = Var(
        initialize=0.5,
        # bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Overall water recovery over all time periods",
    )

    # mp.apparent_recovery = Var(
    #     initialize=0.5,
    #     # bounds=(0, 1),
    #     domain=NonNegativeReals,
    #     units=pyunits.dimensionless,
    #     doc="Apparent water recovery over all time periods",  # from BLM report on OCWD pilot
    # )

    mp.recycle_loop_concentration = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.L,
        doc="Final concentration of the product stream",
    )

    mp.total_feed_vol = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total feed volume produced over all time periods",
    )

    mp.total_permeate_vol = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        bounds=(0, None),
        units=pyunits.m**3,
        doc="Total permeate volume produced over all time periods",
    )

    mp.total_permeate_salt = Var(
        initialize=1,
        domain=NonNegativeReals,
        bounds=(0, None),
        units=pyunits.kg,
        doc="Total permeate salt produced over all time periods",
    )

    mp.permeate_concentration = Var(
        initialize=0.5,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.liter,
        doc="Concentration of salt in the permeate",
    )

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
        initialize=3600,
        domain=NonNegativeReals,
        bounds=(10, None),
        units=pyunits.s,
        doc="Total cycle time including flushing",
    )

    # mp.total_flushing_time = Var(
    #     initialize=3600,
    #     domain=NonNegativeReals,
    #     bounds=(10, None),
    #     units=pyunits.s,
    #     doc="Total flushing time",
    # )

    mp.cycle_time_ratio = Var(
        initialize=1,
        domain=NonNegativeReals,
        # bounds=(0, 1.0001),
        units=pyunits.dimensionless,
        doc="Ratio of total cycle time to filtration time",
    )
    mp.total_flush_volume = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.m**3,
        doc="Total flushing volume over all time periods",
    )

    mp.filtration_ramp_rate = Var(
        # mp.filtration_set,
        initialize=0.25,
        bounds=(0.01, 50),
        units=pyunits.bar / pyunits.min,
        doc="Pressure ramp rate during filtration steps",
    )


def add_multiperiod_constraints(mp, cc_configuration=None):
    """
    Add constraints to the multiperiod model.
    """

    # Get all filtration time blocks
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]
    bf = blks[mp.TIME.last()]

    # Define time step variables for use in constraints based on operation mode
    dt_filt = b0.fs.dead_volume.accumulation_time[0]
    if bf.fs.operation_mode == "flushing":
        dt_flush = bf.fs.flushing.flushing_time
    if bf.fs.operation_mode == "flushing_with_filtration":
        dt_flush = mp.flushing.flushing_time / mp.flushing_points

    # RO membrane area should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane area equality through all filtration periods",
    )
    def ro_membrane_area_constraint(b, t):
        if t == b.TIME.first() or blks[t].fs.operation_mode == "flushing":
            return Constraint.Skip
        return blks[t].fs.RO.area == b0.fs.RO.area

    # RO membrane length should be the same across all time periods - except flushing
    @mp.Constraint(
        mp.TIME,
        doc="RO membrane length equality through all filtration periods",
    )
    def ro_membrane_length_constraint(b, t):
        if t == b.TIME.first() or blks[t].fs.operation_mode == "flushing":
            return Constraint.Skip
        return blks[t].fs.RO.length == b0.fs.RO.length

    # Dead volume should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Dead volume equality through all filtration periods",
    )
    def equal_dead_volume_constraint(b, t):
        if t == b.TIME.first():
            return Constraint.Skip
        elif (
            blks[t].fs.operation_mode == "flushing" and blks[t].fs.ro_model_with_hold_up
        ):
            return (
                blks[t].fs.dead_volume.volume[0, "Liq"]
                == b0.fs.dead_volume.volume[0, "Liq"] + b0.fs.RO.feed_side.volume
            )
        else:
            return (
                blks[t].fs.dead_volume.volume[0, "Liq"]
                == b0.fs.dead_volume.volume[0, "Liq"]
            )

    @mp.Constraint(
        mp.TIME, doc="Delta dead volume equality through all filtration periods"
    )
    def equal_delta_dead_volume_constraint(b, t):
        return (
            blks[t].fs.dead_volume.volume[0, "Liq"]
            == blks[t].fs.dead_volume.delta_state.volume[0, "Liq"]
        )

    # Recycle rate should be the same across all time periods
    @mp.Constraint(
        mp.TIME,
        doc="Recycle rate equality through all filtration periods",
    )
    def equal_recycle_rate(b, t):
        if t == b.TIME.first() or blks[t].fs.operation_mode == "flushing":
            return Constraint.Skip
        elif blks[t].fs.operation_mode == "filtration":
            return (
                blks[t].fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
                == b0.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
            )
        elif blks[t].fs.operation_mode == "flushing_with_filtration":
            return (
                blks[t].fs.conduit_feed.properties[0].flow_vol_phase["Liq"]
                == b0.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]
            )

    @mp.Constraint(doc="Global dead volume constraint")
    def global_dead_volume_constraint(b):
        if b0.fs.ro_model_with_hold_up:
            return (
                b0.fs.dead_volume.volume[0, "Liq"]
                == b.dead_volume_to_area_ratio
                * b0.fs.RO.area
                * b.dead_volume_to_area_multiplier
                * b.pipe_to_module_ratio
            )
        else:
            return b0.fs.dead_volume.volume[
                0, "Liq"
            ] == b0.fs.RO.area * b.dead_volume_to_area_ratio * b.dead_volume_to_area_multiplier * (
                1 + b.pipe_to_module_ratio
            )

    calculate_variable_from_constraint(
        b0.fs.dead_volume.volume[0, "Liq"], mp.global_dead_volume_constraint
    )

    if b0.fs.ro_model_with_hold_up:

        @mp.Constraint(doc="Global dead volume constraint")
        def global_ro_volume_constraint(b):
            return (
                b0.fs.RO.feed_side.volume
                == b0.fs.RO.area
                * b.dead_volume_to_area_ratio
                * b.dead_volume_to_area_multiplier
            )

        calculate_variable_from_constraint(
            b0.fs.RO.feed_side.volume, mp.global_ro_volume_constraint
        )

        @mp.Constraint(
            mp.TIME,
            doc="Dead volume equality through all filtration periods",
        )
        def equal_ro_volume_constraint(b, t):
            if t == b.TIME.first() or blks[t].fs.operation_mode == "flushing":
                return Constraint.Skip
            else:
                return blks[t].fs.RO.feed_side.volume == b0.fs.RO.feed_side.volume

    # Density at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Density equality between end of flushing and start of filtration"
    )
    def cycle_end_density_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            == bf.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"]
        )

    # Mass fraction at the start of cycle should be the same as end of flushing
    # (Initial condition and after flushing)
    @mp.Constraint(
        doc="Mass fraction equality between end of flushing and start of filtration"
    )
    def cycle_end_mass_frac_constraint(b):
        return (
            b0.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
            == bf.fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp[
                "Liq", "NaCl"
            ]
        )

    ### linking ro start stat to final flushed state here
    if b0.fs.ro_model_with_hold_up:
        if bf.fs.operation_mode == "flushing":

            @mp.Constraint(
                b0.fs.RO.difference_elements,
                doc="Density equality between end of flushing and start of filtration",
            )
            def ro_cycle_end_density_constraint(b, i):
                return (
                    b0.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"]
                    == bf.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase[
                        "Liq"
                    ]
                )

            # Mass fraction at the start of cycle should be the same as end of flushing
            # (Initial condition and after flushing)
            @mp.Constraint(
                b0.fs.RO.difference_elements,
                doc="Mass fraction equality between end of flushing and start of filtration",
            )
            def ro_cycle_end_mass_frac_constraint(b, i):
                return (
                    b0.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                        0, i, "Liq", "NaCl"
                    ]
                    == bf.fs.dead_volume.dead_volume.properties_out[
                        0
                    ].mass_frac_phase_comp["Liq", "NaCl"]
                )

    if bf.fs.operation_mode == "flushing_with_filtration":

        @mp.Constraint(
            b0.fs.RO.difference_elements,
            doc="Density equality between end of flushing and start of filtration",
        )
        def ro_cycle_end_density_constraint(b, i):
            return (
                b0.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"]
                == bf.fs.RO.feed_side.properties[0, i].dens_mass_phase["Liq"]
            )

        # Mass fraction at the start of cycle should be the same as end of flushing
        # (Initial condition and after flushing)
        @mp.Constraint(
            b0.fs.RO.difference_elements,
            doc="Mass fraction equality between end of flushing and start of filtration",
        )
        def ro_cycle_end_mass_frac_constraint(b, i):
            return (
                b0.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                    0, i, "Liq", "NaCl"
                ]
                == bf.fs.RO.feed_side.properties[0, i].mass_frac_phase_comp[
                    "Liq", "NaCl"
                ]
            )

    @mp.Expression(mp.TIME, doc="Operation time for each time period")
    def operation_time_points(b, t):
        times = []
        for m in blks[: t + 1]:
            if m.fs.operation_mode == "filtration":
                times.append(dt_filt)
            else:
                times.append(dt_flush)
        return sum(times)

    @mp.Expression(doc="Total flushing time")
    def total_flushing_time(b):
        return sum(dt_flush for _ in b.flushing_set)

    @mp.Constraint(doc="Total filtration time constraint")
    def total_filtration_time_constraint(b):
        return b.total_filtration_time == sum(dt_filt for _ in b.filtration_set)

    @mp.Expression(doc="Change in operating pressure across filtration steps")
    def dp_filtration(b):
        m0 = blks[b.filtration_set.first()]
        # m0 = blks[b.TIME.last()]
        mf = blks[b.filtration_set.last()]
        return (
            mf.fs.P1.control_volume.properties_out[0].pressure
            - m0.fs.P1.control_volume.properties_out[0].pressure
        )

    @mp.Constraint(doc="Filtration ramp rate constraint")
    def filtration_ramp_rate_constraint(b, t):
        return b.filtration_ramp_rate == pyunits.convert(
            b.dp_filtration / b.total_filtration_time,
            to_units=pyunits.bar / pyunits.minute,
        )

    @mp.Expression(mp.filtration_set, doc="Ramp rate for filtration steps")
    def ramp_rate_filtration(b, t):
        p_now = blks[t].fs.P1.control_volume.properties_out[0].pressure
        if t == b.filtration_set.first():
            p_last = (
                blks[b.flushing_set.last()]
                .fs.P1.control_volume.properties_out[0]
                .pressure
            )
        else:
            p_last = blks[t - 1].fs.P1.control_volume.properties_out[0].pressure
        dp = pyunits.convert(p_now - p_last, to_units=pyunits.bar)
        return pyunits.convert(dp / dt_filt, to_units=pyunits.bar / pyunits.minute)

    @mp.Expression(mp.flushing_set, doc="Ramp rate for flushing steps")
    def ramp_rate_flushing(b, t):
        p_now = blks[t].fs.P1.control_volume.properties_out[0].pressure
        if t == b.flushing_set.first():
            p_last = (
                blks[b.filtration_set.last()]
                .fs.P1.control_volume.properties_out[0]
                .pressure
            )
        else:
            p_last = blks[t - 1].fs.P1.control_volume.properties_out[0].pressure
        dp = pyunits.convert(p_now - p_last, to_units=pyunits.bar)
        return pyunits.convert(dp / dt_flush, to_units=pyunits.bar / pyunits.minute)

    @mp.Expression(mp.TIME, doc="Ramp rate for each time period")
    def ramp_rate(b, t):
        if t in b.filtration_set:
            return b.ramp_rate_filtration[t]
        if t in b.flushing_set:
            return b.ramp_rate_flushing[t]

    @mp.Constraint(doc="Total cycle time constraint")
    def total_cycle_time_constraint(b):
        return b.total_cycle_time == b.total_filtration_time + b.total_flushing_time

    @mp.Constraint(doc="Cycle time ratio constraint")
    def cycle_time_ratio_constraint(b):
        return b.cycle_time_ratio == (b.total_filtration_time / b.total_cycle_time)

    @mp.Constraint(doc="Concentration in the recycle loop")
    def recycle_loop_concentration_constraint(b):
        if bf.fs.operation_mode == "flushing":
            return (
                b.recycle_loop_concentration
                == bf.fs.flushing.pre_flushing_concentration
            )
        else:
            return (
                b.recycle_loop_concentration == mp.flushing.pre_flushing_concentration
            )

    # Total permeate water volume
    @mp.Constraint(doc="Total permeate produced over all time periods")
    def total_permeate_vol_constraint(b):
        total_permeate_vol = []
        for m in blks:
            if m.fs.operation_mode == "filtration":
                total_permeate_vol.append(
                    m.fs.product.properties[0].flow_vol_phase["Liq"] * dt_filt
                )
            elif m.fs.operation_mode == "flushing_with_filtration":
                total_permeate_vol.append(
                    m.fs.product.properties[0].flow_vol_phase["Liq"] * dt_flush
                )
        return b.total_permeate_vol == sum(total_permeate_vol)

    # Total permeate salt mass
    @mp.Constraint(doc="Total mass of salt in permeate over all time periods")
    def total_permeate_salt_constraint(b):
        total_permeate_salt = []
        for m in blks:
            if m.fs.operation_mode == "filtration":
                total_permeate_salt.append(
                    m.fs.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
                    * dt_filt
                )
            elif m.fs.operation_mode == "flushing_with_filtration":
                total_permeate_salt.append(
                    m.fs.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
                    * dt_flush
                )
        return b.total_permeate_salt == sum(total_permeate_salt)

    @mp.Constraint(doc="Permeate concentration constraint")
    def permeate_concentration_constraint(b):
        return b.permeate_concentration == pyunits.convert(
            b.total_permeate_salt / b.total_permeate_vol,
            to_units=pyunits.g / pyunits.liter,
        )

    @mp.Constraint(doc="Upper bound on permeate concentration")
    def max_permeate_concentration_constraint(b):
        return b.permeate_concentration <= 0.5

    @mp.Constraint(doc="Total flush volume over all time periods")
    def total_flush_volume_constraint(b):
        total_flush_volume = []
        for m in blks:
            if m.fs.operation_mode == "flushing":
                total_flush_volume.append(m.fs.dead_volume.volume[0, "Liq"])
            elif m.fs.operation_mode == "flushing_with_filtration":
                total_flush_volume.append(
                    m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"] * dt_flush
                )
        return b.total_flush_volume == sum(total_flush_volume)

    if mp.find_component("conduit") is not None:
        mp.conduit.volume_eq = Constraint(
            expr=mp.conduit.volume == mp.total_flush_volume
        )
        mp.conduit.volume.unfix()

    @mp.Constraint(doc="Average product flow rate over all time periods")
    def eq_avg_product_flow_rate(mp):
        return mp.avg_product_flow_rate == pyunits.convert(
            # mp.total_permeate_vol / sum(total_operating_time),
            mp.total_permeate_vol / mp.total_cycle_time,
            to_units=pyunits.m**3 / pyunits.s,
        )

    # Total feed
    @mp.Constraint(doc="Total feed volume over all time periods")
    def total_feed_vol_constraint(b):
        total_feed_vol = []
        for m in blks:
            if m.fs.operation_mode == "filtration":
                total_feed_vol.append(
                    m.fs.raw_feed.properties[0].flow_vol_phase["Liq"] * dt_filt
                )
            elif m.fs.operation_mode == "flushing_with_filtration":
                total_feed_vol.append(
                    m.fs.raw_feed.properties[0].flow_vol_phase["Liq"] * dt_flush
                )
                total_feed_vol.append(
                    m.fs.conduit_feed.properties[0].flow_vol_phase["Liq"] * dt_flush
                )
            elif m.fs.operation_mode == "flushing":
                total_feed_vol.append(m.fs.dead_volume.volume[0, "Liq"])
        return b.total_feed_vol == sum(total_feed_vol)

    @mp.Expression(doc="Avg feed water flow rate over all time periods")
    def avg_feed_flow_rate(b):
        return pyunits.convert(
            b.total_feed_vol / b.total_cycle_time,
            to_units=pyunits.m**3 / pyunits.s,
        )

    # Overall water recovery
    @mp.Constraint(doc="Overall water recovery for system")
    def overall_water_recovery_constraint(b):
        return b.total_permeate_vol == b.overall_recovery * b.total_feed_vol
        # return b.avg_product_flow_rate == b.overall_recovery * b.avg_feed_flow_rate

    deactivate_mp_constraints(mp)


def deactivate_mp_constraints(mp):
    """
    Deactivate constraints on the multiperiod model
    """
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]
    mp.ro_membrane_area_constraint.deactivate()
    mp.ro_membrane_length_constraint.deactivate()
    mp.equal_dead_volume_constraint.deactivate()
    mp.equal_delta_dead_volume_constraint.deactivate()
    mp.equal_recycle_rate.deactivate()
    mp.global_dead_volume_constraint.deactivate()
    
    if b0.fs.ro_model_with_hold_up:
        mp.global_ro_volume_constraint.deactivate()
        mp.equal_ro_volume_constraint.deactivate()

    mp.filtration_ramp_rate_constraint.deactivate()
    mp.max_permeate_concentration_constraint.deactivate()


def fix_overall_water_recovery(mp, overall_water_recovery):
    mp.overall_recovery.fix(overall_water_recovery)
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]
    if mp.flushing_points > 1:
        flushing_block = blks[mp.time_points]
    else:
        flushing_block = blks[-1]
    print("flushing block:", flushing_block.name, flushing_block.fs.operation_mode)
    # Fixed for accumulation time for initialization
    build_acc_time_links = False
    for t, m in enumerate(blks):
        m.fs.dead_volume.accumulation_time.unfix()
        if m.fs.find_component("RO.feed_side.accumulation_time") is not None:
            m.fs.RO.feed_side.accumulation_time.unfix()
            build_acc_time_links = True
        # m.fs.dead_volume.accumulation_time.setlb(1)
        # m.fs.dead_volume.accumulation_time.setub(400)
        iscale.set_scaling_factor(m.fs.dead_volume.accumulation_time, 1e-2)

    # Equal accumulation time across all filtration periods
    @mp.Constraint(mp.TIME, doc="Equal accumulation time across all filtration periods")
    def accumulation_time_cons(b, t):
        if (
            t == b.TIME.first()
            or blks[t].fs.operation_mode == "flushing"
            or flushing_block.name == blks[t].name
        ):
            return Constraint.Skip
        elif blks[t].fs.operation_mode == "filtration":
            return blks[t].fs.dead_volume.accumulation_time[0] == (
                b0.fs.dead_volume.accumulation_time[0]
            )
        elif blks[t].fs.operation_mode == "flushing_with_filtration":
            return blks[t].fs.dead_volume.accumulation_time[0] == (
                flushing_block.fs.dead_volume.accumulation_time[0]
            )

    for c in mp.accumulation_time_cons.values():
        iscale.constraint_scaling_transform(c, 1)
    if build_acc_time_links:

        @mp.Constraint(mp.TIME)
        def ro_accumulation_time_cons(b, t):
            if blks[t].fs.operation_mode == "flushing":
                return Constraint.Skip
            elif blks[t].fs.operation_mode == "filtration":
                return blks[t].fs.RO.feed_side.accumulation_time[0] == (
                    b0.fs.dead_volume.accumulation_time[0]
                )
            elif blks[t].fs.operation_mode == "flushing_with_filtration":
                return blks[t].fs.RO.feed_side.accumulation_time[0] == (
                    flushing_block.fs.dead_volume.accumulation_time[0]
                )

        for c in mp.ro_accumulation_time_cons.values():
            iscale.constraint_scaling_transform(c, 1)
