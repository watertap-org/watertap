import idaes.core.util.scaling as iscale


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
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)

    iscale.set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 1)
    iscale.set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 1)

    iscale.calculate_scaling_factors(m)


def scale_multiperiod_model(mp):

    iscale.set_scaling_factor(mp.dead_volume_to_area_ratio, 1e2)
    iscale.set_scaling_factor(mp.dead_volume_to_area_multiplier, 1)
    iscale.set_scaling_factor(mp.total_cycle_time, 1e-3)
    iscale.set_scaling_factor(mp.total_filtration_time, 1e-3)
    blks = list(mp.get_active_process_blocks())
    b0 = blks[mp.TIME.first()]
    flow_vol = b0.fs.raw_feed.properties[0].flow_vol_phase["Liq"].value
    iscale.set_scaling_factor(mp.total_feed, 1 / (flow_vol * 3600))
    iscale.constraint_scaling_transform(
        mp.total_permeate_constraint, 1 / (flow_vol * 3600)
    )
    iscale.set_scaling_factor(mp.total_permeate, 1 / (flow_vol * 3600))
    iscale.constraint_scaling_transform(mp.total_feed_constraint, 1 / (flow_vol * 3600))

    iscale.set_scaling_factor(mp.avg_product_flow_rate, 1 / (flow_vol))
    iscale.constraint_scaling_transform(mp.eq_avg_product_flow_rate, 1 / (flow_vol))

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
