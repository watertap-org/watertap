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

from watertap.flowsheets.crystallization.crystallizer_TVC import main
from pyomo.environ import value


@pytest.mark.component
def test_crystallizer_TVC_main():
    m = main()
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.41046
    assert (
        pytest.approx(value(m.fs.crystallizer.temperature_operating), rel=1e-3)
        == 319.576
    )
    assert (
        pytest.approx(value(m.fs.heater.cold_side_inlet.temperature[0]), rel=1e-3)
        == 319.15
    )
    assert (
        pytest.approx(value(m.fs.condenser.cold_side_outlet.temperature[0]), rel=1e-3)
        == 491.68
    )
    assert (
        pytest.approx(
            value(m.fs.condenser.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]),
            rel=1e-3,
        )
        == 20
    )
    assert pytest.approx(value(m.fs.SteamEjector.compression_ratio), rel=1e-3) == 2
    assert (
        pytest.approx(
            value(m.fs.separator.recycle.flow_mass_phase_comp[0, "Liq", "H2O"]),
            rel=1e-3,
        )
        == 2279.91
    )
    assert (
        pytest.approx(
            value(m.fs.separator.purge.flow_mass_phase_comp[0, "Liq", "H2O"]), rel=1e-3
        )
        == 23.2430
    )
    assert pytest.approx(value(m.fs.heater.area), rel=1e-3) == 2941.86
