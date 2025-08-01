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
Tests for modified ASM2D flowsheet example.

Based on flowsheet from:

Flores-Alsina X., Gernaey K.V. and Jeppsson, U. "Benchmarking biological
nutrient removal in wastewater treatment plants: influence of mathematical model
assumptions", 2012, Wat. Sci. Tech., Vol. 65 No. 8, pp. 1496-1505
"""

# Some more information about this module
__author__ = "Chenyu Wang"

import pytest

from pyomo.environ import assert_optimal_termination, value
from watertap.flowsheets.activated_sludge.modified_ASM2D_flowsheet import (
    build_flowsheet,
)
from idaes.core.util import DiagnosticsToolbox


# lbianchi-lbl 2024-12-12: adding this as part of watertap-org/watertap#1540
# this was only observed for macOS (EXPERIMENTAL) CI env which is slated for removal soon
@pytest.mark.filterwarnings(
    "ignore:divide by zero encountered in divide:RuntimeWarning"
)
class TestASM2DFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m, res = build_flowsheet()

        m.results = res

        return m

    @pytest.mark.integration
    def test_structure(self, model):
        dt = DiagnosticsToolbox(model)
        assert dt.get_dulmage_mendelsohn_partition() == ([], [], [], [])

    @pytest.mark.component
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model, variable_bounds_violation_tolerance=1e-9)
        dt.assert_no_numerical_warnings()
        assert_optimal_termination(model.results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(0.20921, rel=1e-4)
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_A"]) == pytest.approx(
            2.653e-5, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.7466e-4, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.0e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            7.886e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            4.425e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            1.88e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.677e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            2.7906e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            1.0481e-1, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_K"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            1.02129e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            5.489e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.75925e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            1.0456e-3, rel=1e-4
        )

        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(4.2808e-3, rel=1e-4)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_A"]) == pytest.approx(
            2.65295e-5, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.74656e-4, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.0e-2, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            7.8862e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            4.4252e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            1.88e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.6767e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            2.7906e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            0.10481, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_K"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            4.513e-2, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_H"]) == pytest.approx(
            2.4477, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            0.7921, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            0.04601, rel=1e-4
        )
