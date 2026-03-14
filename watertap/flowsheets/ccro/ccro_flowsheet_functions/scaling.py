import idaes.core.util.scaling as iscale


def overscale_ro(ro, props):
    h2o_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", "H2O")]
    NaCl_scale = props._default_scaling_factors["flow_mass_phase_comp", ("Liq", "NaCl")]
    scales = {"H2O": h2o_scale, "NaCl": NaCl_scale}

    iscale.set_scaling_factor(ro.area, 1e-2)
    iscale.constraint_scaling_transform(ro.eq_area, 1 / 100)
    iscale.set_scaling_factor(ro.width, 1)
    iscale.set_scaling_factor(ro.length, 1)
    stage = ro
    for e in stage.feed_side.velocity:
        iscale.set_scaling_factor(stage.feed_side.velocity[e], 1)

    iscale.constraint_scaling_transform(stage.eq_area, 1 / 10)

    for temp_stream in [
        stage.eq_permeate_isothermal,
        stage.feed_side.eq_equal_temp_interface,
        stage.feed_side.eq_feed_isothermal,
        stage.eq_permeate_outlet_isothermal,
    ]:
        for e in temp_stream:
            iscale.constraint_scaling_transform(temp_stream[e], 1e-2)
    for pressure_stream in [
        stage.eq_permeate_outlet_isobaric,
        stage.feed_side.eq_equal_pressure_interface,
    ]:
        for e in pressure_stream:
            iscale.constraint_scaling_transform(pressure_stream[e], 1e-5)
    for e in stage.eq_pressure_drop:
        iscale.constraint_scaling_transform(stage.eq_pressure_drop[e], 1e-4)
    for e in stage.feed_side.eq_K:
        iscale.constraint_scaling_transform(stage.feed_side.eq_K[e], 1e4)

    for e in stage.feed_side.eq_N_Sh_comp:
        iscale.constraint_scaling_transform(stage.feed_side.eq_N_Sh_comp[e], 1e-2)
    for e in stage.feed_side.eq_N_Re:
        iscale.constraint_scaling_transform(stage.feed_side.eq_N_Re[e], 1e2)

    for e in stage.feed_side.eq_friction_factor:
        iscale.constraint_scaling_transform(stage.feed_side.eq_friction_factor[e], 1e-2)
    for e in stage.feed_side.eq_dP_dx:
        iscale.constraint_scaling_transform(stage.feed_side.eq_dP_dx[e], 1e-3)

    for e in stage.feed_side.eq_equal_flow_vol_interface:
        iscale.constraint_scaling_transform(
            stage.feed_side.eq_equal_flow_vol_interface[e], 1e1
        )

    for e in stage.eq_mass_transfer_term:
        sf = scales[e[-1]]
        iscale.constraint_scaling_transform(stage.eq_mass_transfer_term[e], sf * 10)
    for e in stage.feed_side.mass_transfer_term:
        sf = scales[e[-1]]
        iscale.set_scaling_factor(stage.feed_side.mass_transfer_term[e], sf * 10)
    for e in stage.eq_mass_flux_equal_mass_transfer:
        sf = scales[e[-1]]
        iscale.constraint_scaling_transform(
            stage.eq_mass_flux_equal_mass_transfer[e], sf
        )
    for e in stage.eq_connect_mass_transfer:
        if e[-1] == "H2O":
            sf = sf * 10
        if e[-1] == "NaCl":
            sf = sf * 100
        sf = scales[e[-1]]
        iscale.constraint_scaling_transform(stage.eq_connect_mass_transfer[e], sf)
    for e in stage.eq_recovery_mass_phase_comp:
        sf = scales[e[-1]]
        if e[-1] == "H2O":
            sf = sf * 10
        if e[-1] == "NaCl":
            sf = sf * 100
        iscale.set_scaling_factor(stage.eq_recovery_mass_phase_comp[e], sf)
    for e in stage.eq_permeate_production:
        sf = scales[e[-1]]
        if e[-1] == "H2O":
            sf = sf * 10
        if e[-1] == "NaCl":
            sf = sf * 100
        iscale.constraint_scaling_transform(stage.eq_permeate_production[e], sf)
    for e in stage.eq_flux_mass:
        sf = scales[e[-1]]
        if e[-1] == "H2O":
            sf = sf * 10
        if e[-1] == "NaCl":
            sf = sf * 100
        iscale.constraint_scaling_transform(stage.eq_flux_mass[e], sf)


def scale_flushing_system(m=None):
    """
    Scale flushing model configuration
    """
    flow_mass = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    flow_nacl = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value
    print(f"Flow H2O: {flow_mass}, Flow NaCl: {flow_nacl}")
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_mass, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_nacl, index=("Liq", "NaCl")
    )

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-1)
    overscale_ro(m.fs.RO, m.fs.properties)

    iscale.set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 1)
    iscale.set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 1)
    iscale.constraint_scaling_transform(m.fs.post_flushing_conc_constraint, 1)
    iscale.constraint_scaling_transform(m.fs.pre_flushing_conc_constraint, 1)
    iscale.set_scaling_factor(m.fs.flushing.pre_flushing_concentration, 1)
    iscale.set_scaling_factor(m.fs.flushing.post_flushing_concentration, 1)
    iscale.calculate_scaling_factors(m)


def scale_filtration_system(m):
    """
    Scale filtration model configuration
    """

    flow_mass = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    flow_nacl = m.fs.raw_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value
    print(f"Flow H2O: {flow_mass}, Flow NaCl: {flow_nacl}")
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_mass, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / flow_nacl, index=("Liq", "NaCl")
    )

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-1)
    overscale_ro(m.fs.RO, m.fs.properties)

    iscale.set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 1)
    iscale.set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 1)

    iscale.calculate_scaling_factors(m)


def scale_multiperiod_model(mp):

    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]
    flow_vol = b0.fs.raw_feed.properties[0].flow_vol_phase["Liq"].value
    perm_conc = b0.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value

    ### VARIABLES ###
    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)
    iscale.set_scaling_factor(mp.total_cycle_time, 1e-3)
    iscale.set_scaling_factor(mp.total_filtration_time, 1e-3)
    iscale.set_scaling_factor(mp.filtration_ramp_rate, 1)
    iscale.set_scaling_factor(mp.permeate_concentration, 10)
    iscale.set_scaling_factor(mp.total_feed_vol, 1 / (flow_vol * 3600))
    iscale.set_scaling_factor(mp.total_permeate_vol, 10)
    iscale.set_scaling_factor(mp.total_permeate_salt, 10)
    iscale.set_scaling_factor(mp.avg_product_flow_rate, 1 / flow_vol)

    ### CONSTRAINTS ###
    iscale.constraint_scaling_transform(
        mp.total_feed_vol_constraint, 1 / (flow_vol * 3600)
    )
    iscale.constraint_scaling_transform(mp.eq_avg_product_flow_rate, 1 / flow_vol)
    iscale.constraint_scaling_transform(mp.total_permeate_vol_constraint, 0.1)
    iscale.constraint_scaling_transform(mp.total_permeate_salt_constraint, 0.1)
    iscale.constraint_scaling_transform(mp.permeate_concentration_constraint, 0.1)

    iscale.constraint_scaling_transform(mp.total_filtration_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.total_cycle_time_constraint, 1e-3)
    iscale.constraint_scaling_transform(mp.recycle_loop_concentration_constraint, 1e-1)
    iscale.constraint_scaling_transform(mp.global_dead_volume_constraint, 1e2)
    for c in mp.equal_recycle_rate.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_ro_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.cycle_end_density_constraint.values():
        iscale.constraint_scaling_transform(c, 1)
    for c in mp.ro_cycle_end_density_constraint.values():
        iscale.constraint_scaling_transform(c, 1)
    for c in mp.ro_cycle_end_mass_frac_constraint.values():
        iscale.constraint_scaling_transform(c, 1)
    for c in mp.ro_membrane_area_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.ro_membrane_length_constraint.values():
        iscale.constraint_scaling_transform(c, 1e-1)
    for c in mp.equal_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
    for c in mp.equal_delta_dead_volume_constraint.values():
        iscale.constraint_scaling_transform(c, 1e2)
