from watertap.flowsheets.ccro import PS_RO


from pyomo.environ import Var, units as pyunits


def build(use_hold_up=True, time_steps=11, data_point_skips=40, **kwargs):

    cc_config = PS_RO.validation_configs()

    mp = PS_RO.create_ccro_multiperiod(
        n_time_points=time_steps,
        cc_configuration=cc_config,
        use_ro_with_hold_up=use_hold_up,
    )
    PS_RO.load_validation_data_into_model(
        mp,
        file_path="../validation_data/sine_700-900psi_60s_period.csv",
        data_point_skips=data_point_skips,
    )
    mp.dummy_var = Var(initialize=1.0)
    mp.dummy_var.fix()
    return mp


def solve_model(mp, **kwargs):
    results = PS_RO.solve(mp, use_ipoptv2=False)
    return results


if __name__ == "__main__":
    mp = build()
    results = PS_RO.solve(mp)
