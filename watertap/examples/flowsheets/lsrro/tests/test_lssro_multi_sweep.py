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

import math
import os
import pytest

import pandas as pd

from watertap.examples.flowsheets.lsrro.multi_sweep import run_case

_this_file_path = os.path.dirname(os.path.abspath(__file__))

# NOTE: we used to test up to 5 stages, but those are
#       excluded by the rule below, so no point in running
_test_cases = list(range(1, 3))


@pytest.mark.parametrize("number_of_stages", _test_cases)
@pytest.mark.integration
def test_against_multisweep(number_of_stages, tmp_path):
    csv_file_name = f"{number_of_stages}_stage_results_LSRRO.csv"
    csv_test_file_name = os.path.join("./", csv_file_name)
    csv_baseline_file_name = os.path.join(
        _this_file_path, "parameter_sweep_baselines", csv_file_name
    )
    run_case(number_of_stages, 2, output_filename=csv_test_file_name)

    baseline = pd.read_csv(csv_baseline_file_name).astype(float).T.to_dict()
    test = pd.read_csv(csv_test_file_name).astype(float).T.to_dict()

    assert len(baseline) == len(test)

    for k, base in baseline.items():
        # Don't test those cases which have too many stages
        for s in range(1, number_of_stages + 1):
            if math.isclose(base[f"Membrane Area-Stage {s}"], 1.0, abs_tol=1e-4):
                print(
                    f"Skipping feed concentration {base['# Feed Concentration']}, recovery {base['Volumetric Recovery Rate']}; stage {s} has membrane area of 1"
                )
                break
        else:  # no break
            assert pytest.approx(base, nan_ok=True, rel=1e-02, abs=1e-07) == test[k]
