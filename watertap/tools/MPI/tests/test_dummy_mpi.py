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
import numpy as np
from watertap.tools.MPI.dummy_mpi import DummyCOMM


@pytest.mark.integration
def test_dummy_mpi():
    comm = DummyCOMM
    rank = comm.Get_rank()
    size = comm.Get_size()

    n_arr = 3
    dummy_array1 = np.arange(n_arr)
    comm.Bcast(dummy_array1)
    returned_array = comm.bcast(dummy_array1)
    returned_list = comm.allgather(dummy_array1)
    recvbuf = (np.zeros(n_arr), [n_arr])
    comm.Gatherv(sendbuf=dummy_array1, recvbuf=recvbuf, root=0)

    assert comm is DummyCOMM
    assert rank == 0
    assert size == 1
    assert returned_array is dummy_array1
    assert returned_list == [dummy_array1]
    assert (recvbuf[0] == dummy_array1).all
