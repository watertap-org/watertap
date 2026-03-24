import pytest
from pyomo.environ import value, assert_optimal_termination
from watertap.flowsheets.ccro.analysis_scripts import ccro_setup as ccro_setup

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

# 3.24.2026 -KAS
# LCOWs for BW, SW, PW: [0.2209698655431355, 0.58445275299291, 0.8553982136011674]
base_case_lcows = {
    "BW": 0.22096,
    "SW": 0.58445,
    "PW": 0.85540,
}


def assert_build_with_fixed_recovery(**kwargs):
    """
    build_with_fixed_recovery from ccro_setup but will assert_optimal_termination
    """
    recovery = kwargs.get("recovery", None)
    assert recovery is not None, "Recovery must be specified in kwargs"

    mp = ccro_setup.build(**kwargs)
    mp.overall_recovery.fix(recovery)
    results = ccro_setup.solve_model(mp)
    assert_optimal_termination(results)
    return mp


@pytest.mark.parametrize("water", base_case_kwargs.keys())
def test_base_cases(water):
    """
    Test conditions for base case against changes in model formulation.
    Unless an argument is listed, the default is supposed to result in an optimal solve
    for the given water type.
    """

    kwargs = base_case_kwargs[water]
    lcow = base_case_lcows[water]
    mp = assert_build_with_fixed_recovery(**kwargs)
    assert value(mp.costing.LCOW) == pytest.approx(lcow, rel=1e-4)
