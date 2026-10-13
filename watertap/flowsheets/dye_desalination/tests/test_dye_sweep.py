#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest
import os
import tempfile
from watertap.flowsheets.dye_desalination.dye_sweep import (
    run_analysis,
)

_PR_REPRESENTATIVE_CASES = {1, 6, 11, 12, 16, 19}


def _case_params(case_numbers):
    return [
        pytest.param(
            case_num,
            marks=(
                pytest.mark.volatile
                if case_num in _PR_REPRESENTATIVE_CASES
                else pytest.mark.nightly
            ),
            id=f"case_{case_num}",
        )
        for case_num in case_numbers
    ]


@pytest.mark.parametrize("case_num", _case_params(range(1, 12)))
@pytest.mark.integration
def test_dye_sweep(case_num, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    # test every other sweep with RO
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.close()
    results, params, model = run_analysis(
        case_num,
        nx=1,
        interpolate_nan_outputs=False,
        output_filename=temp.name,
    )
    os.remove(temp.name)
    os.chdir(cwd)
    return


@pytest.mark.parametrize("case_num", _case_params(range(12, 20)))
@pytest.mark.integration
def test_dye_sweep2(case_num, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    # test only with RO
    run_analysis(case_num, nx=1, interpolate_nan_outputs=False)
    os.chdir(cwd)
    return


@pytest.mark.integration
def test_out_of_range_cases():
    # test a case number that is not included in the dye sweep file
    with pytest.raises(Exception):
        run_analysis(30, 1, interpolate_nan_outputs=False)
    return
