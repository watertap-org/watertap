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
Tests for the METAB + ASM1 flowsheet example.
"""

# Some more information about this module
__author__ = "Marcus Holly"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.flowsheets.METAB.METAB_BSM1_flowsheet import (
    main,
)


class TestMETABASM1Flowsheet:
    @pytest.fixture(scope="class")
    @classmethod
    def model(self):
        m, res = main()

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
        assert value(model.fs.Treated.flow_vol[0]) == pytest.approx(1.0902e-4, rel=1e-4)
        assert value(model.fs.Treated.temperature[0]) == pytest.approx(308.15, rel=1e-4)
        assert value(model.fs.Treated.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Treated.conc_mass_comp[0, "S_I"]) == pytest.approx(
            8.6312e-2, rel=1e-5
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_S"]) == pytest.approx(
            3.9661e-2, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_I"]) == pytest.approx(
            9.5593e-3, rel=1e-3
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_S"]) == pytest.approx(
            6.421e-2, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_BH"]) == pytest.approx(
            1.2043e-11, abs=1e-6
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_BA"]) == pytest.approx(
            8.877e-12, abs=1e-6
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_P"]) == pytest.approx(
            8.6168e-12, abs=1e-6
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_O"]) == pytest.approx(
            4.49e-4, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            1.1787e-10, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_NH"]) == pytest.approx(
            2.0578e-1, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "S_ND"]) == pytest.approx(
            5.2525e-3, rel=1e-2
        )
        assert value(model.fs.Treated.conc_mass_comp[0, "X_ND"]) == pytest.approx(
            4.9807e-3, rel=1e-2
        )
        assert value(model.fs.Treated.alkalinity[0]) == pytest.approx(
            4.6668e-2, rel=1e-2
        )

        # Sludge stream
        assert value(model.fs.Sludge.flow_vol[0]) == pytest.approx(2.3245e-6, rel=1e-4)
        assert value(model.fs.Sludge.temperature[0]) == pytest.approx(308.15, rel=1e-4)
        assert value(model.fs.Sludge.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_I"]) == pytest.approx(
            8.6312e-2, rel=1e-5
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_S"]) == pytest.approx(
            3.9661e-2, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_I"]) == pytest.approx(
            4.8936, rel=1e-3
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_S"]) == pytest.approx(
            32.8708, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_BH"]) == pytest.approx(
            3.4696e-8, rel=1e-3
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_BA"]) == pytest.approx(
            3.31293e-8, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_P"]) == pytest.approx(
            3.300396e-8, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_O"]) == pytest.approx(
            4.49016e-4, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            2.87103e-8, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_NH"]) == pytest.approx(
            2.05778e-1, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "S_ND"]) == pytest.approx(
            5.2525e-3, rel=1e-2
        )
        assert value(model.fs.Sludge.conc_mass_comp[0, "X_ND"]) == pytest.approx(
            2.54976, rel=1e-2
        )
        assert value(model.fs.Sludge.alkalinity[0]) == pytest.approx(
            4.6668e-2, rel=1e-2
        )
