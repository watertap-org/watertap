from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits
from idaes.core.util.model_statistics import degrees_of_freedom


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
    recycle_flow_bounds=(1, 50),
    overall_water_recovery=0.5,
    accumulation_time=5,
    flushing_efficiency=0.25,
    **kwargs,
):
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
    )
    if not fixed_setup:
        CCRO.setup_optimization(
            mp,
            overall_water_recovery=overall_water_recovery,
            min_cycle_time_hr=min_cycle_time_hr,
            max_cycle_time_hr=max_cycle_time_hr,
            recycle_flow_bounds=recycle_flow_bounds,
        )
    else:
        CCRO.setup_optimization(
            mp,
            overall_water_recovery=0.5,
            max_cycle_time_hr=5,
            recycle_flow_bounds=(1, 50),
            # use_ro_with_hold_up=True,
        )
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
    solve_model(mp)
    first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].unfix()
    # first_block.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix(
    #     cross_flow * pyunits.L / pyunits.s
    # )
    return mp


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

    # mp = build(feed_tds=35, A_comp=1.5, B_comp=0.1)
    # mp.overall_recovery.fix(0.5)
    # results = solve_model(mp)
    # mp.flushing.flushing_efficiency.fix(0.4)
    # results = solve_model(mp)

    mp = build_for_flush_eff(
        overall_recovery=0.5,
        feed_tds=35,
        A_comp=1.5,
        B_comp=0.1,
        flushing_efficiency=0.6,
    )
    mp.flushing.flushing_efficiency.display()
    mp.permeate_concentration.display()
    # for x in [0.25, 0.4, 0.5, 0.6, 0.75, 0.9]:
    #     mp.flushing.flushing_efficiency.fix(x)
    #     results = solve_model(mp)
