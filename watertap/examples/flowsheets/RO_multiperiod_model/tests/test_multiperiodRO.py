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

from idaes.core.solvers import get_solver
import watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO as mpro

solver = get_solver()


class TestMPmodel:
    @pytest.fixture(scope="class")
    def mp_model(self):
        m = mpro.create_multiperiod_swro_model(n_time_points=2)
        return m

    @pytest.mark.unit
    def test_model_creation(self, mp_model):
        m = mp_model

        # test time blocks exist
        assert len(m.pyomo_model.blocks) == 2
        assert m.pyomo_model.active

    def test_model_properties(self, mp_model):
        m = mp_model

        # test sub-blocks
        assert hasattr(m.pyomo_model.blocks[0], "process")
        assert hasattr(m.pyomo_model.blocks[1].process.ro_mp.fs, "properties")
        assert hasattr(m.pyomo_model.blocks[1].process.ro_mp.fs, "costing")
        assert hasattr(m.pyomo_model.blocks[1].process.ro_mp.fs, "feed")
        assert hasattr(m.pyomo_model.blocks[1].process.ro_mp.fs, "RO")
