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

class DummyMPI:
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
    def allgather(value):
        return [value]

    @staticmethod
    def Gatherv(sendbuf, recvbuf, root=0):
        receive_arr = recvbuf[0]
        receive_sizes = recvbuf[1]
        assert sum(receive_sizes) == len(receive_arr)
        assert sendbuf.size == receive_arr.size
        receive_arr[:] = sendbuf[:]

if __name__ == '__main__':

    dummy_comm = DummyMPI()
    rank = dummy_comm.Get_rank()
    size = dummy_comm.Get_size()
    dummy_comm.Bcast(np.zeros(2))
    returned_list = dummy_comm.allgather(np.zeros(2))

    print("dummy_comm = ", dummy_comm)
    print("rank = ", rank)
    print("size = ", size)
    print("returned_list = ", returned_list)
