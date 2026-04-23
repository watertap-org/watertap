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
            4.75097e-5, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.102e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.0e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            1.418e-2, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            1.3076e-5, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.1e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.4735e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            9.4822e-4, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            1.326e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            3.356e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            2.0912e-2, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            5.260e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            2.168e-4, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            1.6414e-3, abs=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            1.2306e-3, rel=1e-4
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            7.993e-2, rel=1e-4
        )
        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(
            7.336e-3, rel=1e-4
        )
        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(4.2808e-3, rel=1e-4)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_A"]) == pytest.approx(
            4.751e-5, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.102e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.0e-2, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            1.418e-2, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            1.3076e-5, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.100e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.473e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            9.4822e-4, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            5.859e-2, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.4966, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            0.941598, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            0.2337, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            9.873e-3, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            7.488e-2, abs=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            5.415e-2, rel=1e-4
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.5611, rel=1e-4
        )
        assert value(model.fs.Sludge.alkalinity[0]) == pytest.approx(7.336e-3, rel=1e-4)
