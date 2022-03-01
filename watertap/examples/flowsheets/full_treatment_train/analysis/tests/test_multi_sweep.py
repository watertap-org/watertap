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
import pytest
import os
from pyomo.environ import value
from watertap.examples.flowsheets.full_treatment_train.analysis import (flowsheet_NF,
                                                               flowsheet_NF_no_bypass,
                                                               flowsheet_single_stage,
                                                               flowsheet_two_stage,
                                                               flowsheet_NF_two_stage,
                                                               flowsheet_softening,
                                                               flowsheet_softening_two_stage)
from watertap.examples.flowsheets.full_treatment_train.analysis.multi_sweep import *
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep

@pytest.mark.component
def test_multi_sweep():
    # Start MPI communicator
    comm, rank, num_procs = _init_mpi()
    nx = 2
    RO_type = '0D'

    fail_counter = 0
    failed_cases = []
    for case_num in [1,2,3,4,8,9]: # range(1,10):
        try:
            # raise ValueError()
            global_results, sweep_params = run_analysis(case_num, nx, RO_type)
        except:
            fail_counter += 1
            failed_cases.append(case_num)

    if fail_counter > 0:
        pytest.fail(f"multi_sweep.py failed for cases {failed_cases} ")
