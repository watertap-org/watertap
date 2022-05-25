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

# TODO: Load and spot-check parameter sweep analysis corresponding to
#       s*_stage/results_LSRRO.csv's for stages 1-8

# TODO: Maybe we want to randomly generate which of the various test
#       scenarios we run, so we get coverage in expectation.

import os
import sys
import platform
import pytest

import pandas as pd

from watertap.examples.flowsheets.lsrro.multi_sweep import run_case

_this_file_path = os.path.dirname(os.path.abspath(__file__))

_test_cases = list(range(1, 9))

# comment out these lines if you want to run the entire baseline
_supported_systems = ["Linux", "Windows"]
_number_of_python_versions = 3
_python_version_index = sys.version_info.minor % _number_of_python_versions

_this_platform = platform.system()
try:
    _platform_index = _supported_systems.index(_this_platform)
except ValueError:
    _platform_index = None

if _platform_index is None:
    _test_cases = []
else:
    _test_cases = [
        test_case
        for idx, test_case in enumerate(_test_cases)
        if (idx % len(_supported_systems) == _platform_index)
        and (idx % _number_of_python_versions == _python_version_index)
    ]
# END code to comment out


@pytest.fixture(scope="module")
def cleanup_parameter_sweep_files():
    yield
    ## TODO: remove any files that could
    ##       have been created by the
    ##       parameter sweep


@pytest.mark.parametrize("test_case_index", _test_cases)
def test_against_multisweep(test_case_index, cleanup_parameter_sweep_files):
    run_case(test_case_index, 2)
