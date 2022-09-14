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
import numpy as np


class DummyCOMM:

    rank = 0
    size = 1

    @staticmethod
    def Get_rank():
        return 0

    @staticmethod
    def Get_size():
        return 1

    @staticmethod
    def Barrier():
        pass

    @staticmethod
    def Bcast(array, root=0):
        pass

    @staticmethod
    def bcast(array, root=0):
        return array

    @staticmethod
    def allgather(value):
        return [value]

    @staticmethod
    def Gatherv(sendbuf=None, recvbuf=(None, None), root=0):
        assert len(recvbuf) == 2
        assert isinstance(recvbuf, tuple)
        receive_arr = recvbuf[0]
        receive_sizes = recvbuf[1]

        assert isinstance(receive_arr, (np.ndarray, list))
        assert len(receive_arr) == sum(receive_sizes)
        receive_arr[:] = sendbuf[:]
