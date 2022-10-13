#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for ASM2D flowsheet example.
"""

# Some more information about this module
__author__ = "Alejandro Garciadiego, Andrew Lee"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.activated_sludge.ASM2D_flowsheet_noPHA import (
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

    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(0.1564, rel=1e-4)
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_A"]) == pytest.approx(
            1.09e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.43e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.00e-2, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            1.5e-2, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            1.51e-2, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            1.06e-8, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.49e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.25e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            1.27e-9, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_H"]) == pytest.approx(
            7.53e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.98e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            9.99e-11, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            1.16e-9, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            5.18e-8, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            4.35e-8, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            4.0e-8, rel=1e2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            2.73e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            1.32e-3, rel=1e-2
        )
        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(7.86e-3, rel=1e-2)

        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(0.05708, rel=1e-2)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(298.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_A"]) == pytest.approx(
            1.09e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_F"]) == pytest.approx(
            2.43e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            3.00e-2, rel=1e-3
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            1.5e-2, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            1.51e-2, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            1.06e-08, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            4.49e-4, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.25e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            1.03e-7, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_H"]) == pytest.approx(
            0.385, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            0.101, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            1.36e-9, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            1.2e-8, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            1.95e-6, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            3.64e-6, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            1.60e-6, rel=1e2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            0.140, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            0.677, rel=1e-2
        )
        assert value(model.fs.Sludge.alkalinity[0]) == pytest.approx(7.86e-3, rel=1e-2)
