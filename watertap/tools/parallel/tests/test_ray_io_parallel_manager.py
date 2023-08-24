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

from watertap.tools.parallel.ray_io_parallel_manager import (
    RayIoParallelManager,
)


class TestRayIoParallelManager:
    @pytest.mark.component
    def test_combine_data_with_peers(self):
        parallel_manager = RayIoParallelManager(number_of_subprocesses=1)
        assert parallel_manager.combine_data_with_peers("test_data") == ["test_data"]

    @pytest.mark.component
    def test_gather_arrays_to_root(self):
        parallel_manager = RayIoParallelManager(number_of_subprocesses=1)

        sendbuf = [1, 2, 3]
        recvbuf_spec = [[0, 0, 0], [1, 2]]
        parallel_manager.gather_arrays_to_root(sendbuf, recvbuf_spec)
        assert recvbuf_spec[0] == [1, 2, 3]

    @pytest.mark.component
    def test_sum_values_and_sync(self):
        parallel_manager = RayIoParallelManager(number_of_subprocesses=1)

        sendbuf = [[1, 2, 3]]
        recvbuf = [[]]
        parallel_manager.sum_values_and_sync(sendbuf, recvbuf)
        assert recvbuf[0] == [1, 2, 3]
