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

from pyomo.common.dependencies import attempt_import

from watertap.tools.parallel.concurrent_futures_parallel_manager import (
    ConcurrentFuturesParallelManager,
)
from watertap.tools.parallel.single_process_parallel_manager import (
    SingleProcessParallelManager,
)
from watertap.tools.parallel.multiprocessing_parallel_manager import (
    MultiprocessingParallelManager,
)

from watertap.tools.parallel.parallel_manager_factory import create_parallel_manager

ray, ray_available = attempt_import("ray", defer_check=False)
if ray_available:
    from watertap.tools.parallel.ray_io_parallel_manager import (
        RayIoParallelManager,
    )


def do_build(base=10):
    return [base, -1]


def do_execute(local_parameters, base, multiplier):
    return [(base + p) * multiplier for p in local_parameters]


class TestParallelManager:
    @pytest.mark.component
    def test_single_process(self):
        all_parameters = [1, 2, 3, 4]

        parallel_manager = SingleProcessParallelManager()
        parallel_manager.scatter(do_build, dict(), do_execute, all_parameters)
        execution_results = parallel_manager.gather()

        # there should be no fan-out; all results should come back from a single node
        assert len(execution_results) == 1
        assert execution_results[0].parameters == all_parameters
        assert execution_results[0].results == [-11, -12, -13, -14]

    @pytest.mark.component
    def test_single_process_with_build_kwargs(self):
        all_parameters = [1, 2, 3, 4]

        parallel_manager = SingleProcessParallelManager()
        parallel_manager.scatter(do_build, {"base": 100}, do_execute, all_parameters)
        execution_results = parallel_manager.gather()

        # there should be no fan-out; all results should come back from a single node
        assert len(execution_results) == 1
        assert execution_results[0].parameters == all_parameters
        assert execution_results[0].results == [-101, -102, -103, -104]

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1, 4, 8])
    def test_multiple_subprocesses(self, number_of_subprocesses):
        all_parameters = [1, 2, 3, 4, 5]

        parallel_manager = ConcurrentFuturesParallelManager(number_of_subprocesses)

        # the parallel manager shouldn't kick off more subprocesses than can do work
        number_of_subprocesses_used = min(number_of_subprocesses, len(all_parameters))

        parallel_manager.scatter(do_build, {"base": 100}, do_execute, all_parameters)
        execution_results = parallel_manager.gather()

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
        assert all_results == [-101, -102, -103, -104, -105]

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1, 2])
    def test_multiple_subprocesses_with_multiprocessing(self, number_of_subprocesses):
        all_parameters = [1, 2, 3, 4, 5]

        parallel_manager = MultiprocessingParallelManager(number_of_subprocesses)
        # the parallel manager shouldn't kick off more subprocesses than can do work
        number_of_subprocesses_used = min(number_of_subprocesses, len(all_parameters))

        parallel_manager.scatter(do_build, {"base": 100}, do_execute, all_parameters)
        execution_results = parallel_manager.gather()
        # print(execution_results)
        # running in async, so there should be as many results as parametrs
        assert len(execution_results) == len(all_parameters)

        # verify that each subprocess did actual work
        for i in range(len(execution_results)):
            assert len(execution_results[i].parameters) > 0
            assert len(execution_results[i].parameters) == len(
                execution_results[i].results
            )

        # verify that the full list of results matches the expected
        all_results = [r for result in execution_results for r in result.results]
        assert all_results == [-101, -102, -103, -104, -105]

    @pytest.mark.component
    @pytest.mark.parametrize("number_of_subprocesses", [1])
    def test_rayio_multiprocessing(self, number_of_subprocesses):
        all_parameters = [1, 2, 3, 4, 5]
        if ray_available:
            parallel_manager = RayIoParallelManager(number_of_subprocesses)
            # the parallel manager shouldn't kick off more subprocesses than can do work
            number_of_subprocesses_used = min(
                number_of_subprocesses, len(all_parameters)
            )

            parallel_manager.scatter(
                do_build, {"base": 100}, do_execute, all_parameters
            )
            execution_results = parallel_manager.gather()
            # running in async, so there should be as many results as parametrs
            assert len(execution_results) == len(all_parameters)

            # verify that each subprocess did actual work
            for i in range(len(execution_results)):
                assert len(execution_results[i].parameters) > 0
                assert len(execution_results[i].parameters) == len(
                    execution_results[i].results
                )

            # verify that the full list of results matches the expected
            all_results = [r for result in execution_results for r in result.results]
            assert all_results == [-101, -102, -103, -104, -105]
        else:
            pass

    @pytest.mark.component
    def test_create_parallel_manager(
        self,
    ):
        pm = create_parallel_manager(
            number_of_subprocesses=2,
            parallel_back_end="ConcurrentFutures",
        )
        assert isinstance(pm, ConcurrentFuturesParallelManager)
        pm = create_parallel_manager(
            number_of_subprocesses=2,
            parallel_back_end="MultiProcessing",
        )
        assert isinstance(pm, MultiprocessingParallelManager)
        if ray_available:
            pm = create_parallel_manager(
                number_of_subprocesses=2,
                parallel_back_end="RayIo",
            )

            assert isinstance(pm, RayIoParallelManager)
        else:
            with pytest.raises(ModuleNotFoundError):
                pm = create_parallel_manager(
                    number_of_subprocesses=2,
                    parallel_back_end="RayIo",
                )
