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
import matplotlib.pyplot as plt
from pyomo.environ import value

from idaes.core.solvers import get_solver
import idaes.core.util.model_statistics as mod_stats
import watertap.examples.flowsheets.RO_multiperiod_model.pricetaker_multiperiod_RO as pt

solver = get_solver()


class TestMPflowsheet:
    @pytest.fixture(scope="class")
    def frame(self):
        base_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        data_path = "dagget_CA_LMP_hourly_2015.csv"
        sample_path = os.path.join(base_path, data_path)
        sample_data = pt._get_lmp(2, sample_path)

        m = pt.build_flowsheet(2)
        m, t_blocks = pt.set_objective(m, sample_data)

        return m, t_blocks, sample_data

    @pytest.mark.unit
    def test_set_objective(self, frame):
        m, t_blocks, lmp = frame

        # check if block-level objectives are properly set
        assert mod_stats.number_activated_objectives(m) == 1
        assert m.obj.is_expression_type()
        assert t_blocks[0].weighted_LCOW.is_expression_type()
        assert t_blocks[0].water_prod.is_expression_type()
        assert t_blocks[0].energy_consumption.is_expression_type()

        # check if there are 2 blocks (1 for each timestep)
        assert len(t_blocks) == 2

        # verify initial pressure is fixed
        assert t_blocks[0].ro_mp.previous_pressure.is_fixed()
        assert pytest.approx(55e5, rel=1e-3) == value(
            t_blocks[0].ro_mp.previous_pressure
        )

    @pytest.mark.unit
    def test_solve(self, frame):
        m, t_blocks, lmp = frame
        m, results = pt.solve(m)
        assert results.solver.status == "ok"
        assert results.solver.termination_condition == "optimal"

    @pytest.mark.unit
    def test_visualize(self, frame):
        m, t_blocks, lmp = frame
        pt.visualize_results(m, t_blocks, lmp)

        assert plt.fignum_exists(1)
        plt.close("all")

        assert pytest.approx(0.33453, rel=1e-3) == m.obj()

    @pytest.mark.component
    def test_full_dataset(self):
        # bring in full set of sample data
        base_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        data_path = "dagget_CA_LMP_hourly_2015.csv"
        sample_path = os.path.join(base_path, data_path)

        sample_data = pt._get_lmp(-1, sample_path)
        assert len(sample_data) == 8759
