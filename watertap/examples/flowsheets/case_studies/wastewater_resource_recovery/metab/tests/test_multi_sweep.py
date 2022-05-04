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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.multi_sweep import (
    main,
)

pytest_parameterize_list = []
for case_num in [1, 2, 3, 4, 5, 6, 7]:
    pytest_parameterize_list.append(case_num)


@pytest.mark.parametrize("case_num", pytest_parameterize_list)
@pytest.mark.integration
def test_multi_sweep(case_num, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    nx = 1
    global_results, sweep_params = main(case_num, nx, interpolate_nan_outputs=False)
    os.chdir(cwd)
