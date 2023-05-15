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
        "lcow": [m.fs.costing.LCOW.value, 0.12178636975214299],
        "pressure": [m.fs.NF.pump.outlet.pressure[0].value / 1e5, 6.364618667926368],
        "area": [m.fs.NF.nfUnit.area.value, 174.97127919329446],
        "recovery": [
            m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"].value * 100,
            53.60319873234126,
        ],
        "bypass": [
            m.fs.by_pass_splitter.split_fraction[0, "bypass"].value * 100,
            4.127411072051739,
        ],
    }
    for key, (solve, testval) in test_dict.items():
        assert round(solve, 4) == round(testval, 4)
