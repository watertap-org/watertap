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

import numpy as np
import pandas as pd

from watertap.examples.flowsheets.lsrro.multi_sweep import run_case

from .gha_divider import get_test_cases_subset

_this_file_path = os.path.dirname(os.path.abspath(__file__))

_test_cases = list(range(1, 6))

# comment out this line if you want to run the entire baseline
_test_cases = get_test_cases_subset(_test_cases)


@pytest.mark.parametrize("test_case_index", _test_cases)
@pytest.mark.integration
def test_against_multisweep(test_case_index, tmp_path):
    csv_file_name = f"{test_case_index}_stage_results_LSRRO.csv"
    csv_test_file_name = os.path.join(tmp_path, csv_file_name)
    csv_baseline_file_name = os.path.join(
        _this_file_path, "parameter_sweep_baselines", csv_file_name
    )
    run_case(test_case_index, 2, output_filename=csv_test_file_name)

    baseline = pd.read_csv(csv_baseline_file_name).astype(float).T.to_dict()
    test = pd.read_csv(csv_test_file_name).astype(float).T.to_dict()

    assert len(baseline) == len(test)

    for k in test:
        assert pytest.approx(baseline[k], nan_ok=True, rel=1e-02, abs=1e-07) == test[k]
