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
Tests for electroNP flowsheet example.
"""

# Some more information about this module
__author__ = "Chenyu Wang"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.electroNP.BSM2_electroNP_no_bioP import (
    main,
)


class TestElectroNPFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m, res = m, results = main(has_electroNP=True)

        m.results = res

        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_structure(self, model):
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 0
        assert_optimal_termination(model.results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_results(self, model):
        # Treated water
        assert value(model.fs.electroNP.treated.flow_vol[0]) == pytest.approx(
            0.00285, rel=1e-4
        )
        assert value(model.fs.electroNP.treated.temperature[0]) == pytest.approx(
            308.15, rel=1e-4
        )
        assert value(model.fs.electroNP.treated.pressure[0]) == pytest.approx(
            101325, rel=1e-4
        )
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_A"]
        ) == pytest.approx(8.48827, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_F"]
        ) == pytest.approx(22.1026, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_I"]
        ) == pytest.approx(0.05804, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_N2"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_NH4"]
        ) == pytest.approx(1.42581, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_NO3"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_O2"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "S_PO4"]
        ) == pytest.approx(3.5846, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_AUT"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_H"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_I"]
        ) == pytest.approx(0.30187, rel=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_PAO"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_PHA"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_PP"]
        ) == pytest.approx(0, abs=1e-4)
        assert value(
            model.fs.electroNP.treated.conc_mass_comp[0, "X_S"]
        ) == pytest.approx(0.076503, rel=1e-4)
