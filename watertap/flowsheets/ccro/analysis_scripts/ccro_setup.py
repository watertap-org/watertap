from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale


def build(
    time_steps=10,
    n_flushing_points=5,
    use_hold_up=True,
    feed_tds=5,
    fixed_setup=False,
    A_comp=1.5,  # LMH/bar; default for SWRO
    B_comp=0.1,  # LMH; default for SWRO
    min_cycle_time_hr=10 / 60,  # 10 minutes
    max_cycle_time_hr=1,  # 1 hour
    cross_flow=2,
    osmotic_overpressure=2,
    recycle_flow_bounds=(1, 100),
    overall_water_recovery=0.5,
    accumulation_time=5,
    flushing_efficiency=0.25,
    use_interval_initializer=True,
    **kwargs,
):

    if feed_tds <= 10:
        # Use brackish water RO parameters if salinity <=10 g/L
        A_comp = 5
        B_comp = 0.5

    if feed_tds > 70:
        use_interval_initializer = False

    if feed_tds >= 50:
        min_cycle_time_hr = 5 / 60  # 5 minutes

    if feed_tds >= 75:
        min_cycle_time_hr = 1 / 60  # 1 minute

    cc_config = CCROConfiguration()
    cc_config["A_comp"] = A_comp * (
        pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar
    )
    cc_config["B_comp"] = B_comp * (pyunits.liter / pyunits.m**2 / pyunits.hour)
    cc_config["raw_feed_conc"] = feed_tds * pyunits.g / pyunits.L
    cc_config["osmotic_overpressure"] = osmotic_overpressure * pyunits.dimensionless
    cc_config["accumulation_time"] = accumulation_time * pyunits.second
    cc_config["flushing_efficiency"] = flushing_efficiency * pyunits.dimensionless

    mp = CCRO.create_ccro_multiperiod(
        n_time_points=time_steps,
        n_flushing_points=n_flushing_points,
        include_costing=True,
        use_ro_with_hold_up=use_hold_up,
        cc_configuration=cc_config,
        use_interval_initializer=use_interval_initializer,
    )
    # if not fixed_setup:
    CCRO.setup_optimization(
        mp,
        overall_water_recovery=overall_water_recovery,
        min_cycle_time_hr=min_cycle_time_hr,
        max_cycle_time_hr=max_cycle_time_hr,
        recycle_flow_bounds=recycle_flow_bounds,
    )
    # else:
    #     CCRO.setup_optimization(
    #         mp,
    #         overall_water_recovery=0.5,
    #         max_cycle_time_hr=5,
    #         recycle_flow_bounds=(1, 50),
    #         # use_ro_with_hold_up=True,
    #     )
    # CCRO.fix_optimization_dofs(
    #     mp,
    #     overal_water_recovery=0.5,
    #     add_water_recovery_objective=True,
    #     membrane_area=200 * pyunits.meter**2,
    #     membrane_length=1 * pyunits.meter,
    #     recycle_rate=10 * pyunits.L / pyunits.s,
    #     flushing_efficiency=0.8,
    # )

    # CCRO.print_results_table(mp)
    blks = list(mp.get_active_process_blocks())
    mp.overall_recovery.fix()  # Unfixed with times fixed should get 1 DOF!
    first_block = blks[0]
    first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix()
    # mp.cost_objective.deactivate()
    # mp.cost_product_objective.activate()
    solve_model(mp)

    mp.ramp_rate.display()

    # first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
    # mp.max_permeate_concentration_constraint.activate()
    mp.permeate_concentration.setub(0.5)

    solve_model(mp)

    first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
    print(
        "CCRO initalized, solving for optimal design point now with all constraints active and DOF free"
    )

    return mp


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

    for i in iscale.list_unscaled_variables(m):
        print("Var with no scale:", i)
    for i in iscale.list_unscaled_constraints(m):
        print("Constraint with no scale:", i)
    for var, scale in iscale.badly_scaled_var_generator(m):
        print(
            "Badly scaled variable:",
            var.name,
            var.value,
            iscale.get_scaling_factor(var),
        )
    return cond_number


def build_for_flush_eff(overall_recovery=0.5, **kwargs):

    mp = build(**kwargs)
    mp.overall_recovery.fix(overall_recovery)
    solve_model(mp)

    return mp


def solve_model(mp, **kwargs):
    results = CCRO.solve(mp)
    CCRO.print_results_table(mp)
    return results


if __name__ == "__main__":

    # mp = build(
    #     feed_tds=35,
    #     A_comp=1.5,
    #     B_comp=0.1,
    #     osmotic_overpressure=2,
    #     overall_water_recovery=0.4,
    # )

    mp = build(
        feed_tds=35,
        A_comp=1.5,
        B_comp=0.1,
        overall_water_recovery=0.45,
        time_steps=20,
        n_flushing_points=5,
        # osmotic_overpressure=3,
        # overall_water_recovery=0.2,
        use_interval_initializer=True,
        min_cycle_time_hr=0.16667,
    )
    mp.overall_recovery.fix(0.45)
    solve_model(mp)
    mp.ramp_rate.display()
