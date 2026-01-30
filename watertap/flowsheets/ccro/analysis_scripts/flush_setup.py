from watertap.flowsheets.ccro import PS_RO


from pyomo.environ import Var, units as pyunits


from watertap.flowsheets.ccro.utils.cc_configuration import CCROConfiguration


def build(use_hold_up=True, time_steps=11, add_flushing=False, **kwargs):

    cc_config = CCROConfiguration()  # validation_configs()
    cc_config["raw_feed_conc"] = 35 * pyunits.g / pyunits.L
    # cc_config["membrane_area"] = 100 * pyunits.m**2
    cc_config["raw_feed_flowrate"] = 6 * pyunits.L / pyunits.s  # validation_configs()

    mp = PS_RO.create_ccro_multiperiod(
        n_time_points=time_steps,
        cc_configuration=cc_config,
        use_ro_with_hold_up=use_hold_up,
    )
    PS_RO.create_flush_setup(mp)
    PS_RO.solve(mp, use_ipoptv2=False)
    if add_flushing:
        PS_RO.add_flushing_unit(mp, cc_configuration=cc_config)
    mp.dummy_var = Var(initialize=1.0)
    mp.dummy_var.fix()
    return mp


def solve_model(mp, **kwargs):
    results = PS_RO.solve(mp, use_ipoptv2=False)
    return results


if __name__ == "__main__":
    mp = build()
    results = PS_RO.solve(mp)
