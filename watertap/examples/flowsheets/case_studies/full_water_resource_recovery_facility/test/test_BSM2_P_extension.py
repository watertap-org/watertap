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
import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2_P_extension import (
    main,
    solve,
)


# class TestFullFlowsheet:
#     @pytest.mark.requires_idaes_solver
#     @pytest.fixture(scope="class")
#     def system_frame(self):
#         m, res = main(bio_P=False)
#         m.results = res
#         return m
#
#     @pytest.mark.requires_idaes_solver
#     @pytest.mark.integration
#     def test_structure(self, system_frame):
#         assert_units_consistent(system_frame)
#         assert degrees_of_freedom(system_frame) == 0
#         assert_optimal_termination(system_frame.results)
#
#     @pytest.mark.requires_idaes_solver
#     @pytest.mark.component
#     def test_solve(self, system_frame):
#         m = system_frame
#
#         assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
#             0.2422, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
#             5.21698e-07, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
#             0.00029232, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
#             0.057450, rel=1e-3
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
#         ) == pytest.approx(0.058906, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
#         ) == pytest.approx(0.00012576, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
#         ) == pytest.approx(0.0097215, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
#         ) == pytest.approx(0.00772288, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
#         ) == pytest.approx(0.495041, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
#             0.37011, rel=1e-2
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
#         ) == pytest.approx(0.01912677, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
#         ) == pytest.approx(0.130732, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
#         ) == pytest.approx(0.00053442, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
#             0.012031, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
#             0.0116256, rel=1e-3
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
#         ) == pytest.approx(0.00805345, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
#         ) == pytest.approx(2.632e-06, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
#         ) == pytest.approx(0.0027335, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
#             0.00018726, rel=1e-3
#         )


class TestFullFlowsheetBioPTrue:
    @pytest.mark.requires_idaes_solver
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        m.results = res
        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.2422, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            5.524735e-07, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.000290589, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.0587986, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.00014743, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.0077378, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.0076636, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.0045396, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.37011, rel=1e-2
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.020824, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.159079, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.00051909, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.013848, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.012709, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.0114519, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(5.0524e-06, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.0038119, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00022427, rel=1e-3
        )
