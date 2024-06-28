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
Tests for ASM2D flowsheet example.

Based on flowsheet from:

Flores-Alsina X., Gernaey K.V. and Jeppsson, U. "Benchmarking biological
nutrient removal in wastewater treatment plants: influence of mathematical model
assumptions", 2012, Wat. Sci. Tech., Vol. 65 No. 8, pp. 1496-1505
"""

# Some more information about this module
__author__ = "Andrew Lee, Alejandro Garciadiego"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.activated_sludge.ASM2D_flowsheet import (
    build_flowsheet,
)


class TestASM2DFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m, res = build_flowsheet()

        m.results = res

        return m

    @pytest.mark.integration
    def test_structure(self, model):
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 0
        assert_optimal_termination(model.results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(0.20921, rel=1e-4)
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_A"]) == pytest.approx(
            5.4560e-5, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.5736e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            13.837e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            8.98349e-6, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.3924e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            7.8744e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            2.7958e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            1.3277e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            36.331e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            21.1808e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            2.111e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            8.981e-5, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            0.66913e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            1.15897e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            76.489e-3, rel=1e-4
        )
        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(
            7.243e-3, rel=1e-4
        )
        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(4.2808e-3, rel=1e-4)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_A"]) == pytest.approx(
            5.4560e-5, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.5736e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            13.837e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            8.9834e-6, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.3924e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            7.8744e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            2.7958e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            58.6744e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.6200, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            953.686e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            93.795e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            4.0913e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            30.522e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            50.995e-3, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.4078, rel=1e-4
        )
        assert value(model.fs.Sludge.alkalinity[0]) == pytest.approx(
            7.2430e-3, rel=1e-4
        )
