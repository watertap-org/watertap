#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
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

from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.activated_sludge.UConn_WRRF import (
    build_flowsheet,
    set_operating_conditions,
    scale_flowsheet,
    initialize_flowsheet,
    solve_flowsheet,
)
from watertap.flowsheets.activated_sludge.UConn_WRRF import ASMModel


class TestUConnFlowsheetASM3:
    @pytest.fixture(scope="class")
    @classmethod
    def model(self):
        m = build_flowsheet(asm_model=ASMModel.asm3)

        set_operating_conditions(m, asm_model=ASMModel.asm3)
        scale_flowsheet(m)

        initialize_flowsheet(m)

        return m

    @pytest.mark.integration
    def test_structure(self, model):
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 0
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.integration
    def test_solve(self, model):
        res = solve_flowsheet(model)
        assert_optimal_termination(res)
        dt = DiagnosticsToolbox(model)
        warnings, _ = dt._collect_numerical_warnings()
        assert len(warnings) == 1
        assert "WARNING: 1 Variable at or outside bounds (tol=0.0E+00)" in warnings

    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(
            1.06747685185185, rel=1e-4
        )
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(288.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-5
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0.02888, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            0.0036856, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NOX"]) == pytest.approx(
            0.005157, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O"]) == pytest.approx(
            9.1727e-5, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_S"]) == pytest.approx(
            0.00025491, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_A"]) == pytest.approx(
            0.13149, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.6305, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.4627, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            0.2136, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.6305, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_STO"]) == pytest.approx(
            0.3125, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.0305, rel=1e-3
        )

        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(0.00495, rel=1e-2)
