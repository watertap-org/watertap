from pyomo.environ import value, assert_optimal_termination, check_optimal_termination
from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale

from watertap.kurby.whatever import build_with_fixed_recovery


def build(
    time_steps=10,
    n_flushing_points=5,
    use_hold_up=True,
    feed_tds=5,
    A_comp=1.5,  # LMH/bar; default for SWRO
    B_comp=0.1,  # LMH; default for SWRO
    total_cycle_time_lb=10,  # minutes
    total_cycle_time_ub=60,  # minutes
    osmotic_overpressure=2,
    overall_water_recovery=0.5,
    accumulation_time=5,  # seconds
    flushing_efficiency=0.25,
    flushing_time_lb=10,  # seconds
    use_interval_initializer=True,
    recycle_flowrate=10,  # L/s; default in cc_config is 10
    recycle_flowrate_lb=1,  # L/s
    recycle_flowrate_ub=100,  # L/s
    cycle_time_ratio_lb=0.8,
    cycle_time_ratio_ub=0.999999,
    permeate_concentration_lb=0.001,  # g/L
    permeate_concentration_ub=0.5,  # g/L
    flushing_efficiency_lb=0.1,
    flushing_efficiency_ub=0.9999,
    use_perm_conc_target=True,  # activate max perm conc constraint
    rejection_lb=0.985,
    rejection_ub=1,
    use_rejection_target=False,  # activate min rejection constraint
    high_pressure_membrane_cost=False,
    **kwargs,
):

    if feed_tds <= 15:
        # Use brackish water RO parameters if salinity <=15 g/L
        A_comp = 5
        B_comp = 0.5

    if feed_tds > 70:
        use_interval_initializer = False

    if feed_tds >= 50:
        total_cycle_time_lb = 1

    cc_config = CCROConfiguration()
    cc_config["A_comp"] = A_comp * (
        pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar
    )
    cc_config["B_comp"] = B_comp * (pyunits.liter / pyunits.m**2 / pyunits.hour)
    cc_config["raw_feed_conc"] = feed_tds * pyunits.g / pyunits.L
    cc_config["osmotic_overpressure"] = osmotic_overpressure * pyunits.dimensionless
    cc_config["accumulation_time"] = accumulation_time * pyunits.second
    cc_config["flushing_efficiency"] = flushing_efficiency * pyunits.dimensionless
    cc_config["recycle_flowrate"] = recycle_flowrate * pyunits.L / pyunits.s

    if recycle_flowrate_lb > value(cc_config["recycle_flowrate"]):
        cc_config["recycle_flowrate"] = recycle_flowrate_lb * pyunits.L / pyunits.s

    cc_config["max_permeate_concentration"] = permeate_concentration_ub * (
        pyunits.g / pyunits.L
    )
    cc_config["min_overall_rejection"] = rejection_lb * pyunits.dimensionless

    mp = CCRO.create_ccro_multiperiod(
        n_time_points=time_steps,
        n_flushing_points=n_flushing_points,
        include_costing=True,
        use_ro_with_hold_up=use_hold_up,
        cc_configuration=cc_config,
        use_interval_initializer=use_interval_initializer,
        high_pressure_membrane_cost=high_pressure_membrane_cost,
    )

    CCRO.setup_optimization(
        mp,
        overall_water_recovery=overall_water_recovery,
        total_cycle_time_lb=total_cycle_time_lb,
        total_cycle_time_ub=total_cycle_time_ub,
        recycle_flowrate_lb=recycle_flowrate_lb,
        recycle_flowrate_ub=recycle_flowrate_ub,
        flushing_time_lb=flushing_time_lb,
        use_perm_conc_target=use_perm_conc_target,
        use_rejection_target=use_rejection_target,
        flushing_efficiency_lb=flushing_efficiency_lb,
        flushing_efficiency_ub=flushing_efficiency_ub,
    )

    mp.overall_recovery.fix()  # Unfixed with times fixed should get 1 DOF!
    mp.recycle_flowrate.fix()
    # mp.cost_objective.deactivate()
    # mp.cost_product_objective.activate()
    solve_model(mp)

    mp.cycle_time_ratio.setlb(cycle_time_ratio_lb)
    mp.cycle_time_ratio.setub(cycle_time_ratio_ub)

    solve_model(mp)

    mp.recycle_flowrate.unfix()

    print(
        "CCRO initalized, solving for optimal design point now with all constraints active and DOF free"
    )

    return mp


def build_with_fixed_recovery(recovery=0.5, **kwargs):

    mp = build(**kwargs)
    mp.overall_recovery.fix(recovery)
    solve_model(mp)

    return mp


def solve_model(mp, **kwargs):
    results = CCRO.solve(mp)
    CCRO.print_results_table(mp)
    return results


if __name__ == "__main__":

    run_kwargs = {
        "feed_tds": 5,
        "time_steps": 20,
        "overall_water_recovery": 0.5,
        "recovery": 0.75,
        "osmotic_overpressure": 2,
        "total_cycle_time_lb": 10,  # minutes
        "total_cycle_time_ub": 60,  # minutes
        "cycle_time_ratio_lb": 0.8,
        "cycle_time_ratio_ub": 0.999999,
        "flushing_time_lb": 10,
        "rejection_lb": 0.985,
        "rejection_ub": 1,
        # "use_high_pressure_membrane_cost": True,
        "use_high_pressure_membrane_cost": False,
        "use_perm_conc_target": True,
        # "use_perm_conc_target": False,
        # "use_rejection_target": True,
        "use_rejection_target": False,
    }

    mp = build_with_fixed_recovery(**run_kwargs)
