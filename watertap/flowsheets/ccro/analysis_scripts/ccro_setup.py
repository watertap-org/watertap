from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits
from idaes.core.util.model_statistics import degrees_of_freedom


def build(
    time_steps=20,
    n_flushing_points=5,
    use_hold_up=True,
    feed_tds=5,
    fixed_setup=False,
    A_comp=1.5,  # LMH/bar; default for SWRO
    B_comp=0.1,  # LMH; default for SWRO
    min_cycle_time_hr=10 / 60,  # 1 minute
    max_cycle_time_hr=1,  # 5 hours
    cross_flow=2,
    osmotic_overpressure=2,
    recycle_flow_bounds=(1, 50),
    overall_water_recovery=0.5,
    accumulation_time=5,
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


def solve_model(mp, **kwargs):
    results = CCRO.solve(mp)
    CCRO.print_results_table(mp)
    return results


if __name__ == "__main__":

    mp = build(feed_tds=35, A_comp=1.5, B_comp=0.1)
    mp.overall_recovery.fix(0.5)
    results = solve_model(mp)
    # mp.costing.aggregate_flow_electricity.pprint()
    # print(mp.total_utilization())

    # mp = build(
    #     time_steps=10,
    #     feed_tds=35,
    #     A_comp=1.5,
    #     B_comp=0.1,
    #     osmotic_overpressure=3,
    #     accumulation_time=5,
    #     overall_water_recovery=0.5,
    # )
    # mp.filtration_ramp_rate.unfix()
    # mp.overall_recovery.fix(0.5)
    # results = solve_model(mp)
    
    # CCRO.print_results_table(mp_bw)
    # mp_bw.cycle_time_ratio.display()
    # mp_bw.total_cycle_time.display()
    # mp_bw.total_flushing_time.display()
    # mp_bw.total_filtration_time.display()
    mp.filtration_ramp_rate.display()    
    mp.ramp_rate.display()
    mp.filtration_ramp_rate.setub(1)
    results = solve_model(mp)
    mp.filtration_ramp_rate.display()    
    mp.ramp_rate.display()
    mp.dp.display()
    mp.filtration_set.display()
    print(mp.filtration_set.first(), mp.filtration_set.last())
    blks = list(mp.get_active_process_blocks())
    m0 = blks[mp.filtration_set.first()]
    mf = blks[mp.filtration_set.last()]
    from pyomo.environ import value, units as pyunits
    print(value(pyunits.convert(m0.fs.P1.control_volume.properties_out[0].pressure, to_units=pyunits.bar)))
    print(value(pyunits.convert(mf.fs.P1.control_volume.properties_out[0].pressure, to_units=pyunits.bar)))
    print(value(pyunits.convert(blks[mp.flushing_set.last()].fs.P1.control_volume.properties_out[0].pressure, to_units=pyunits.bar)))
    # # mp.total_utilization.display()

    # mp_bw = build(
    #     time_steps=10,
    #     feed_tds=5,
    #     A_comp=5,
    #     B_comp=0.5,
    #     osmotic_overpressure=1,
    #     accumulation_time=20,
    #     overall_water_recovery=0.25,
    # )
    # mp_bw.overall_recovery.fix(0.5)
    # results = solve_model(mp_bw)
    # # CCRO.print_results_table(mp_bw)
    # # mp_bw.cycle_time_ratio.display()
    # # mp_bw.total_cycle_time.display()
    # # mp_bw.total_flushing_time.display()
    # # mp_bw.total_filtration_time.display()
    # mp_bw.filtration_ramp_rate.display()
    # mp_bw.ramp_rate.display()

    # # mp_bw.overall_recovery.unfix()
    # mp_bw.filtration_ramp_rate_constraint.activate()
    # mp_bw.flushing.flushing_efficiency.unfix()
    # results = solve_model(mp_bw)

    # mp_bw.total_filtration_time.display()
    # mp_bw.filtration_ramp_rate.display()
    # mp_bw.ramp_rate.display()

    # mp_sw = build(feed_tds=35, A_comp=1.5, B_comp=0.1)
    # mp_sw.overall_recovery.fix(0.5)
    # results = solve_model(mp_sw)

    # mp_pw = build(feed_tds=75, A_comp=1.5, B_comp=0.1)
    # mp_pw.overall_recovery.fix(0.5)
    # results = solve_model(mp_pw)

    # CCRO.print_results_table(mp_sw)
    # CCRO.print_results_table(mp_pw)
    # for r in [0.55, 0.56, 0.57, 0.58, 0.59, 0.6]:
    #     mp.overall_recovery.fix(r)
    #     results = solve_model(mp)
