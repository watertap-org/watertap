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
import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Block,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.unit_models.crystallizer import Crystallization
from watertap.flowsheets.crystallization.sim_simple_crystallizer import main as cryst_ex
from watertap.property_models.unit_specific.cryst_prop_pack import NaClParameterBlock
from watertap.costing import CrystallizerCostType


class TestCrystallizerBuild:
    @pytest.fixture(scope="class")
    def crystallizer(self):
        m = cryst_ex()
        return m

    @pytest.mark.unit
    def test_model_properties(self, crystallizer):
        m = crystallizer

        # Test basic build
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties, NaClParameterBlock)
        assert isinstance(m.fs.costing, Block)
        assert isinstance(m.fs.crystallizer, Crystallization)

        # Test port
        assert isinstance(m.fs.crystallizer.inlet, Port)
        assert isinstance(m.fs.crystallizer.vapor, Port)
        assert isinstance(m.fs.crystallizer.solids, Port)
        assert isinstance(m.fs.crystallizer.outlet, Port)

        # Test consting
        assert isinstance(m.fs.crystallizer.costing, Block)
        assert isinstance(m.fs.crystallizer.costing.capital_cost, Var)

        var_str_list = [
            "total_capital_cost",
            "total_operating_cost",
        ]
        for var_str in var_str_list:
            var = getattr(m.fs.costing, var_str)
            assert isinstance(var, Var)

        # Test the crystallizer properties
        # test configrations
        assert len(m.fs.crystallizer.config) == 4
        assert not m.fs.crystallizer.config.dynamic
        assert not m.fs.crystallizer.config.has_holdup
        assert m.fs.crystallizer.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list
        assert "NaCl" in m.fs.properties.component_list
        assert (
            m.fs.crystallizer.costing.config["costing_method_arguments"]["cost_type"]
            == CrystallizerCostType.mass_basis
        )

    # @pytest.mark.component
    def test_final_simulation_conditions(self, crystallizer):
        m = crystallizer
        assert degrees_of_freedom(m) == 0
        assert value(m.fs.crystallizer.pressure_operating) == pytest.approx(
            11992, rel=1e-3
        )
        assert value(
            m.fs.crystallizer.product_volumetric_solids_fraction
        ) == pytest.approx(0.0874, rel=1e-3)
        assert value(m.fs.crystallizer.dens_mass_magma) == pytest.approx(
            184.93, rel=1e-3
        )
        assert value(
            m.fs.crystallizer.vapor.flow_mass_phase_comp[0.0, "Vap", "H2O"]
        ) == pytest.approx(21.505, rel=1e-3)
        assert value(
            m.fs.crystallizer.solids.flow_mass_phase_comp[0.0, "Sol", "NaCl"]
        ) == pytest.approx(4.0815, rel=1e-3)
