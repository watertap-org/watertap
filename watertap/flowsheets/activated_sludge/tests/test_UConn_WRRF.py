#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for UConn WRRF flowsheet example.

"""

__author__ = "Adam Atia"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.activated_sludge import UConn_WRRF as uconn 
from watertap.flowsheets.activated_sludge.UConn_WRRF import ASMModel



class TestUConnFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m = uconn.build_flowsheet(asm_model=ASMModel.asm3)

        uconn.set_operating_conditions(m, asm_model=ASMModel.asm3)
        uconn.scale_flowsheet(m)

        uconn.initialize_flowsheet(m)

        uconn.scale_flowsheet(m)

        res = uconn.solve_flowsheet(m)

        m.results = res

        return m

    @pytest.mark.integration
    def test_structure(self, model):
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 0
        assert_optimal_termination(model.results)

    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(0.20904, rel=1e-4)
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-5
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_S"]) == pytest.approx(
            8.89e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            4.39e-3, rel=1e-3
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            1.88e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_BH"]) == pytest.approx(
            9.78e-3, rel=1e-3
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_BA"]) == pytest.approx(
            5.73e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_P"]) == pytest.approx(
            1.73e-3, rel=1e-2
        )
        # S_O is slightly off, but probably due to differences in injection rates
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O"]) == pytest.approx(
            4.49e-4, rel=1e-2
        )
        # Slightly off in last significant digit
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            10.14e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH"]) == pytest.approx(
            1.86e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_ND"]) == pytest.approx(
            6.88e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_ND"]) == pytest.approx(
            1.35e-5, rel=1e-2
        )
        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(4.13e-3, rel=1e-2)

        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(4.457e-3, rel=1e-4)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-5
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_S"]) == pytest.approx(
            8.89e-4, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            2247e-3, rel=1e-3
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            96.8e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_BH"]) == pytest.approx(
            5004e-3, rel=1e-3
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_BA"]) == pytest.approx(
            292e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_P"]) == pytest.approx(
            884e-3, rel=1e-2
        )
        # S_O is slightly off, but probably due to differences in injection rates
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O"]) == pytest.approx(
            4.49e-4, rel=1e-2
        )
        # Slightly off in last significant digit
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            10.14e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH"]) == pytest.approx(
            1.86e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_ND"]) == pytest.approx(
            6.88e-4, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_ND"]) == pytest.approx(
            6.92e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.alkalinity[0]) == pytest.approx(4.13e-3, rel=1e-2)
