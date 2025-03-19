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
Tests for full Water Resource Recovery Facility with ASM2d and ADM1
(WRRF; a.k.a., wastewater treatment plant) flowsheet example with ASM1 and ADM1.
The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.
"""

# Some more information about this module
__author__ = "Chenyu Wang"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.anaerobic_digester.modified_ADM1_flowsheet_with_translators import (
    main,
)
from watertap.core.solvers import get_solver

solver = get_solver()


@pytest.mark.requires_idaes_solver
class TestADM1BioPFalse:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=False)
        m.results = res
        return m

    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.003, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            8.4935, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            22.163, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.05745, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(2.0167, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(65.479, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            1.0654, rel=1e-2
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.71312, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(1.2023, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            14.875, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            3.3158, rel=1e-3
        )


@pytest.mark.requires_idaes_solver
class TestADM1BioPTrue:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        m.results = res
        return m

    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.003, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            9.1528, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            23.908, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.05745, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(1.9658, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(4.0233, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            1.2776, rel=1e-2
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.92538, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.63549, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            14.881, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.85963, rel=1e-3
        )
