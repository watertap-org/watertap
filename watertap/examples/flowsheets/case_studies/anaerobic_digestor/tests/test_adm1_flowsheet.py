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
Tests for ADM1 flowsheet example.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

"""

# Some more information about this module
__author__ = "Alejandro Garciadiego"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.anaerobic_digestor.ADM1_flowsheet import (
    build_flowsheet,
)


class TestADM1Flowsheet:
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

        assert value(model.fs.R1.liquid_outlet.flow_vol[0]) == pytest.approx(
            0.0019675, rel=1e-2
        )
        assert value(model.fs.R1.liquid_outlet.temperature[0]) == pytest.approx(
            308.15, rel=1e-2
        )
        assert value(model.fs.R1.liquid_outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-2
        )
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_su"]
        ) == pytest.approx(1.193e-2, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_aa"]
        ) == pytest.approx(5.31e-3, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_fa"]
        ) == pytest.approx(9.862e-2, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_va"]
        ) == pytest.approx(1.16e-2, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_bu"]
        ) == pytest.approx(0.0132, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_pro"]
        ) == pytest.approx(0.01578, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_ac"]
        ) == pytest.approx(0.17319, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_h2"]
        ) == pytest.approx(2.35e-7, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_ch4"]
        ) == pytest.approx(0.0541, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_IC"]
        ) == pytest.approx(0.15436 * 12, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_IN"]
        ) == pytest.approx(0.13022 * 14, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "S_I"]
        ) == pytest.approx(0.32869, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_c"]
        ) == pytest.approx(0.3086, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_ch"]
        ) == pytest.approx(0.02794, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_pr"]
        ) == pytest.approx(0.10257, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_li"]
        ) == pytest.approx(0.02948, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_su"]
        ) == pytest.approx(0.42016, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_aa"]
        ) == pytest.approx(1.1791, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_fa"]
        ) == pytest.approx(0.2430, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_c4"]
        ) == pytest.approx(0.4319, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_pro"]
        ) == pytest.approx(0.13730, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_ac"]
        ) == pytest.approx(0.76056, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_h2"]
        ) == pytest.approx(0.3170, rel=1e-2)
        assert value(
            model.fs.R1.liquid_outlet.conc_mass_comp[0, "X_I"]
        ) == pytest.approx(25.617, rel=1e-2)
        assert value(model.fs.R1.liquid_outlet.anions[0]) == pytest.approx(
            2e-2, rel=1e-2
        )
        assert value(model.fs.R1.liquid_outlet.cations[0]) == pytest.approx(
            4e-2, rel=1e-2
        )
