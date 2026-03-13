#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest

from pyomo.environ import value

from watertap.flowsheets.crystallization.crystallizer_live_steam_with_condenser_chiller import (
    main,
)


@pytest.mark.component
def test_crystallizer_live_steam_with_condenser_chiller_main():
    m = main()
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.22166
    assert (
        pytest.approx(value(m.fs.crystallizer.inlet.temperature[0]), rel=1e-3) == 323.15
    )
    assert (
        pytest.approx(value(m.fs.crystallizer.temperature_operating), rel=1e-3)
        == 319.58
    )
    assert (
        pytest.approx(
            value(m.fs.mixer.recycle.flow_mass_phase_comp[0, "Liq", "H2O"]), rel=1e-3
        )
        == 2279.92
    )
    assert (
        pytest.approx(
            value(m.fs.heater.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]),
            rel=1e-3,
        )
        == 18.14
    )
    assert pytest.approx(value(m.fs.heater.area), rel=1e-3) == 215.10
