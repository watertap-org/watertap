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

from watertap.flowsheets.crystallization.crystallizer_MVR import main
from pyomo.environ import value


@pytest.mark.component
def test_crystallizer_MVR_main():
    m = main()
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.7863
    assert (
        pytest.approx(value(m.fs.crystallizer.inlet.temperature[0]), rel=1e-3)
        == 323.150
    )
    assert (
        pytest.approx(value(m.fs.crystallizer.temperature_operating), rel=1e-3)
        == 319.576
    )
    assert pytest.approx(value(m.fs.compressor.pressure_ratio), rel=1e-3) == 1.203411
    assert (
        pytest.approx(
            value(m.fs.separator.recycle.flow_mass_phase_comp[0, "Liq", "H2O"]),
            rel=1e-3,
        )
        == 2279.918
    )
    assert (
        pytest.approx(
            value(m.fs.distillate.flow_mass_phase_comp[0, "Liq", "H2O"]), rel=1e-3
        )
        == 15.68956
    )
