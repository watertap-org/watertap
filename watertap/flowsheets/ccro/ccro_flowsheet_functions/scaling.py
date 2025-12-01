import idaes.core.util.scaling as iscale


def scale_flushing_system(m=None):
    """
    Scale flushing model configuration
    """

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)

    iscale.set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 10)
    iscale.set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 10)
    iscale.constraint_scaling_transform(m.fs.post_flushing_conc_constraint, 1e-1)
    iscale.constraint_scaling_transform(m.fs.pre_flushing_conc_constraint, 1e-1)
    iscale.set_scaling_factor(m.fs.flushing.pre_flushing_concentration, 1e-1)
    iscale.set_scaling_factor(m.fs.flushing.post_flushing_concentration, 1e-1)

    iscale.calculate_scaling_factors(m)


def scale_filtration_system(m):
    """
    Scale filtration model configuration
    """

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)

    iscale.set_scaling_factor(m.fs.dead_volume.dead_volume.mass_frac_phase_comp, 10)
    iscale.set_scaling_factor(m.fs.dead_volume.delta_state.mass_frac_phase_comp, 10)

    iscale.calculate_scaling_factors(m)


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
