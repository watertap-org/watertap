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
from watertap.examples.flowsheets.nf_dspmde.nf_with_bypass import main


@pytest.mark.component
def test_main():
    m = main()
    test_dict = {
        "lcow": [m.fs.costing.LCOW.value, 0.13769517724604913],
        "pressure": [m.fs.NF.pump.outlet.pressure[0].value / 1e5, 6.702808007030784],
        "area": [m.fs.NF.nfUnit.area.value, 419.77504107728066],
        "recovery": [
            m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"].value * 100,
            94.99999897892856,
        ],
        "bypass": [
            m.fs.by_pass_splitter.split_fraction[0, "bypass"].value * 100,
            10.58832033562513,
        ],
    }
    for key, (solve, testval) in test_dict.items():
        assert round(solve, 4) == round(testval, 4)
