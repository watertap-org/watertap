from pyomo.environ import value, assert_optimal_termination
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
    min_flushing_time=10,
    use_interval_initializer=True,
    recycle_flowrate=10,  # default in cc_config is 10
    cycle_time_ratio_bounds=(0.8, 0.99),
    permeate_concentration_bounds=(0.001, 0.5),
    use_perm_conc_target=True,
    rejection_bounds=(0.9, 1),
    use_rejection_target=False,
    high_pressure_membrane_cost=False,
    **kwargs,
):

    if feed_tds <= 15:
        # Use brackish water RO parameters if salinity <=15 g/L
        A_comp = 5
        B_comp = 0.5

    if feed_tds > 70:
        use_interval_initializer = False

    # if feed_tds >= 50:
    #     min_cycle_time_hr = 5 / 60  # 5 minutes

    if feed_tds >= 50:
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
    cc_config["recycle_flowrate"] = recycle_flowrate * pyunits.L / pyunits.s

    if recycle_flow_bounds[0] > value(cc_config["recycle_flowrate"]):
        cc_config["recycle_flowrate"] = recycle_flow_bounds[0] * pyunits.L / pyunits.s

    cc_config["max_permeate_concentration"] = permeate_concentration_bounds[1] * (
        pyunits.g / pyunits.L
    )
    cc_config["min_overall_rejection"] = rejection_bounds[0] * pyunits.dimensionless

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
        min_cycle_time_hr=min_cycle_time_hr,
        max_cycle_time_hr=max_cycle_time_hr,
        recycle_flow_bounds=recycle_flow_bounds,
        min_flushing_time=min_flushing_time,
        use_perm_conc_target=use_perm_conc_target,
        use_rejection_target=use_rejection_target,
    )

    mp.overall_recovery.fix()  # Unfixed with times fixed should get 1 DOF!
    mp.recycle_flowrate.fix()
    # mp.cost_objective.deactivate()
    # mp.cost_product_objective.activate()
    solve_model(mp)

    mp.cycle_time_ratio.setlb(cycle_time_ratio_bounds[0])
    mp.cycle_time_ratio.setub(cycle_time_ratio_bounds[1])

    solve_model(mp)

    mp.permeate_concentration.setlb(permeate_concentration_bounds[0])
    mp.permeate_concentration.setub(permeate_concentration_bounds[1])
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

    base_case_kwargs = {
    "BW": {
        "feed_tds": 5,
        "overall_water_recovery": 0.5,
        "recovery": 0.9,
    },
    "SW": {
        "feed_tds": 35,
        "overall_water_recovery": 0.5,
        "recovery": 0.5,
    },
    "PW": {
        "feed_tds": 75,
        "overall_water_recovery": 0.4,
        "recovery": 0.4,
        "high_pressure_membrane_cost": True,
    },
}
    
    kwargs = base_case_kwargs["SW"]
    mp = build_with_fixed_recovery(**kwargs)
