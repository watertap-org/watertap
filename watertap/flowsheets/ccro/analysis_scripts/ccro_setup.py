from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import units as pyunits


def build(time_steps=11, feed_tds=5, **kwargs):
    cc_config = CCROConfiguration()
    cc_config["raw_feed_conc"] = feed_tds * pyunits.g / pyunits.L
    mp = CCRO.create_ccro_multiperiod(
        n_time_points=time_steps, include_costing=True, cc_configuration=cc_config
    )
    CCRO.setup_optimization(
        mp,
        overall_water_recovery=0.8,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
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
