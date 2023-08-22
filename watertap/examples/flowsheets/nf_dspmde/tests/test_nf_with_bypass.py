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
from pyomo.environ import value
from watertap.examples.flowsheets.nf_dspmde.nf_with_bypass import main


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_main():
    m = main()
    test_dict = {
        "lcow": [m.fs.costing.LCOW, 0.1376],
        "pressure": [m.fs.NF.pump.outlet.pressure[0] / 1e5, 6.702781588429595],
        "area": [m.fs.NF.nfUnit.area, 419.7559215919783],
        "recovery": [
            m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"] * 100,
            94.9999989788293,
        ],
        "bypass": [
            m.fs.by_pass_splitter.split_fraction[0, "bypass"] * 100,
            10.592027658361907,
        ],
    }
    for model_result, testval in test_dict.values():
        assert pytest.approx(testval, rel=1e-3) == value(model_result)
