###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import os
import pytest
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.suboxic_activated_sludge_process import (
    multi_sweep,
)

sweep_list = []
for case_num in [1, 2, 3, 4]:
    sweep_list.append(case_num)


@pytest.mark.parametrize("case_num", sweep_list)
@pytest.mark.integration
def test_multi_sweep(case_num, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    nx = 2
    global_results, sweep_params, m = multi_sweep.run_analysis(
        case_num, nx, interpolate_nan_outputs=False
    )
    os.chdir(cwd)
