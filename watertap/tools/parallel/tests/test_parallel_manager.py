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

from watertap.tools.parallel.parallel_manager import return_arguments_as_list
from watertap.tools.parallel.concurrent_futures_parallel_manager import (
    ConcurrentFuturesParallelManager,
)
from watertap.tools.parallel.single_process_parallel_manager import (
    SingleProcessParallelManager,
)
from watertap.tools.parameter_sweep import ParameterSweep


# top-level so it's picklable - used as the common arg rebuilding function in
# test_rebuilding_common_sweep_args
def build_common_args(base=100):
    return [base]


class MockParamSweep(ParameterSweep):
    """
    Mock implementation of the ParameterSweep class, used for testing.
    """

    SWEEP_COMMON_PARAMS = [100]
    SWEEP_PARAMS = [1, 2, 3, 4, 5, 6]
    EXPECTED_SWEEP_RESULTS = [101, 104, 109, 116, 125, 136]

    def _do_param_sweep(self, common_param, local_parameters):
        """
        Simple implementation intended to be run for each local worker. Goal is to
        remove the bookkeeping complexity and dependencies of the real class but
        still return a result that can be verified.
        """
        return [common_param + param**2 for param in local_parameters]


class TestParallelManager:
    @pytest.mark.component
    def test_single_process_scatter_gather(self):

        param_sweep_instance = MockParamSweep(
            parallel_manager_class=SingleProcessParallelManager
        )

        execution_results = param_sweep_instance.run_scatter_gather(
            MockParamSweep.SWEEP_COMMON_PARAMS,
            None,
            None,
            MockParamSweep.SWEEP_PARAMS,
        )

        # there should be no fan-out; all results should come back from a single node
        assert len(execution_results) == 1
        assert execution_results[0].parameters == MockParamSweep.SWEEP_PARAMS
        assert execution_results[0].results == MockParamSweep.EXPECTED_SWEEP_RESULTS

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1, 2, 3, 4, 8, 16])
    def test_concurrent_futures_scatter_gather(self, number_of_subprocesses):
        param_sweep_instance = MockParamSweep(
            parallel_manager_class=ConcurrentFuturesParallelManager,
            number_of_subprocesses=number_of_subprocesses,
        )

        execution_results = param_sweep_instance.run_scatter_gather(
            MockParamSweep.SWEEP_COMMON_PARAMS,
            None,
            None,
            MockParamSweep.SWEEP_PARAMS,
        )

        # the parallel manager shouldn't kick off more subprocesses than can do work
        number_of_subprocesses_used = min(
            number_of_subprocesses, len(MockParamSweep.SWEEP_PARAMS)
        )

        # there should be a set of local results for each subprocess
        assert len(execution_results) == number_of_subprocesses_used

        # verify that each subprocess did actual work
        for i in range(len(execution_results)):
            assert len(execution_results[i].parameters) > 0
            assert len(execution_results[i].parameters) == len(
                execution_results[i].results
            )

        # verify that the full list of results matches the expected
        all_results = [r for result in execution_results for r in result.results]
        assert all_results == MockParamSweep.EXPECTED_SWEEP_RESULTS

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1, 2])
    @pytest.mark.parametrize(
        "parallel_manager_class",
        [SingleProcessParallelManager, ConcurrentFuturesParallelManager],
    )
    def test_rebuilding_common_sweep_args(
        self, number_of_subprocesses, parallel_manager_class
    ):

        # don't run unnecessary tests
        if (
            number_of_subprocesses > 1
            and parallel_manager_class is SingleProcessParallelManager
        ):
            return

        # number_of_subprocesses will be ignored when the class is the SingleProcessParallelManager
        param_sweep_instance = MockParamSweep(
            parallel_manager_class=parallel_manager_class,
            number_of_subprocesses=number_of_subprocesses,
        )

        execution_results = param_sweep_instance.run_scatter_gather(
            MockParamSweep.SWEEP_COMMON_PARAMS,
            build_common_args,
            {"base": 20},
            MockParamSweep.SWEEP_PARAMS,
        )

        # verify that the full list of results matches the expected
        all_results = [r for result in execution_results for r in result.results]
        expected_results = [21, 24, 29, 36, 45, 56]
        assert all_results == expected_results
