#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from watertap.examples.flowsheets.nf_dspmde.nf import main


@pytest.mark.component
def test_main():
    m = main()
    test_dict = {
        "lcow": [m.fs.costing.LCOW.value, 0.1518574055112589],
        "pressure": [m.fs.NF.pump.outlet.pressure[0].value / 1e5, 6.168174681542402],
        "area": [m.fs.NF.nfUnit.area.value, 470.174600076009],
        "recovery": [
            m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"].value * 100,
            91.67052799785642,
        ],
    }
    for key, (solve, testval) in test_dict.items():
        assert round(solve, 4) == round(testval, 4)
