from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits


def build(time_steps=11, use_hold_up=True, feed_tds=5, fixed_setup=False, **kwargs):
    cc_config = CCROConfiguration()
    cc_config["raw_feed_conc"] = feed_tds * pyunits.g / pyunits.L
    mp = CCRO.create_ccro_multiperiod(
        n_time_points=time_steps,
        include_costing=True,
        use_ro_with_hold_up=use_hold_up,
        cc_configuration=cc_config,
    )
    if fixed_setup == False:
        CCRO.setup_optimization(
            mp,
            overall_water_recovery=0.5,
            max_cycle_time_hr=2,
            recycle_flow_bounds=(1, 50),
        )
    else:
        CCRO.setup_optimization(
            mp,
            overall_water_recovery=0.5,
            max_cycle_time_hr=10,
            recycle_flow_bounds=(1, 50),
        )
        CCRO.fix_optimization_dofs(
            mp,
            overal_water_recovery=0.5,
            add_water_recovery_objective=True,
            membrane_area=200 * pyunits.meter**2,
            membrane_length=1 * pyunits.meter,
            recycle_rate=10 * pyunits.L / pyunits.s,
            flushing_efficiency=0.8,
        )

    results = CCRO.solve(mp)

    return mp


def solve_model(mp, **kwargs):
    results = CCRO.solve(mp)
    return results


if __name__ == "__main__":
    mp = build()
    mp.overall_recovery.fix(0.85)
    results = CCRO.solve(mp)
