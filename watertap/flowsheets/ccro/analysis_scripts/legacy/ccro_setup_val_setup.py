from watertap.flowsheets.ccro import CCRO

from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration
from pyomo.environ import Var, units as pyunits


def build(use_ro_with_hold_up=False, time_steps=100, **kwargs):
    cc_config = CCRO.validation_seedling_configs()
    mp = CCRO.create_ccro_multiperiod(
        n_time_points=time_steps,
        include_costing=False,
        use_ro_with_hold_up=use_ro_with_hold_up,
        cc_configuration=cc_config,
    )
    CCRO.setup_optimization(
        mp,
        overall_water_recovery=0.8,
        max_cycle_time_hr=1,
        recycle_flow_bounds=(0.1, 100),
    )
    CCRO.fix_optimization_dofs(
        mp,
        accumulation_time=21 * 60 / (time_steps - 1),
        target_pressure=18 * pyunits.bar,
        add_initial_pressure_objective=True,
    )
    results = CCRO.solve(mp)
    mp.dummy_var = Var(initialize=1.0)
    mp.dummy_var.fix()
    return mp


def solve_model(mp, **kwargs):
    results = CCRO.solve(mp)
    return results


if __name__ == "__main__":
    mp = build()
    mp.overall_recovery.fix(0.85)
    results = CCRO.solve(mp)
