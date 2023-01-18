###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Objective,
    assert_optimal_termination,
    Block,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product
import idaes.logger as idaeslogger
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

import watertap.examples.flowsheets.ion_exchange.ion_exchange_demo as ixf

import math

__author__ = "Kurban Sitterley"

target_ion = "Ca_2+"
ions = [target_ion]
mass_frac = 1e-4
feed_mass_frac = {target_ion: mass_frac}
solver = get_solver()


class TestIXDemo:
    @pytest.fixture(scope="class")
    def ix_0D(self):

        m = ixf.ix_build(ions)
        return m

    @pytest.mark.unit
    def test_build_model(self, ix_0D):
        m = ix_0D

        # Test basic build
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties, MCASParameterBlock)
        assert isinstance(m.fs.costing, Block)
        assert isinstance(m.fs.feed, Feed)
        assert isinstance(m.fs.ion_exchange, IonExchange0D)
        assert isinstance(m.fs.product, Product)
        assert isinstance(m.fs.regen, Product)

        # Test port
        assert isinstance(m.fs.feed.outlet, Port)
        assert isinstance(m.fs.ion_exchange.inlet, Port)
        assert isinstance(m.fs.ion_exchange.outlet, Port)
        assert isinstance(m.fs.ion_exchange.regen, Port)
        assert isinstance(m.fs.product.inlet, Port)
        assert isinstance(m.fs.regen.inlet, Port)

        # # Test consting
        assert isinstance(m.fs.ion_exchange.costing, Block)
        assert isinstance(m.fs.ion_exchange.costing.capital_cost, Var)
        assert isinstance(m.fs.ion_exchange.costing.fixed_operating_cost, Var)

        var_str_list = [
            "total_capital_cost",
            "maintenance_labor_chemical_operating_cost",
            "total_operating_cost",
        ]
        for var_str in var_str_list:
            var = getattr(m.fs.costing, var_str)
            assert isinstance(var, Var)

        # Test arcs
        arc_dict = {
            m.fs.feed_to_ix: (m.fs.feed.outlet, m.fs.ion_exchange.inlet),
            m.fs.ix_to_product: (m.fs.ion_exchange.outlet, m.fs.product.inlet),
            m.fs.ix_to_regen: (
                m.fs.ion_exchange.regen,
                m.fs.regen.inlet,
            ),
        }
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

        # test configrations
        assert len(m.fs.ion_exchange.config) == 7
        assert not m.fs.ion_exchange.config.dynamic
        assert not m.fs.ion_exchange.config.has_holdup
        assert m.fs.ion_exchange.config.target_ion == "Ca_2+"
        assert m.fs.ion_exchange.ion_exchange_type == IonExchangeType.cation
        assert m.fs.ion_exchange.config.regenerant == RegenerantChem.NaCl

        assert m.fs.ion_exchange.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_specific_operating_conditions(self, ix_0D):
        m = ix_0D
        ixf.set_operating_conditions(m, feed_mass_frac=feed_mass_frac)
        ixf.initialize_system(m)
        assert degrees_of_freedom(m) == 0

        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)
        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.05, abs=1e-6
        )
        assert value(m.fs.product.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.049901566, abs=1e-6
        )
        assert value(m.fs.regen.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.016666, abs=1e-6
        )
        assert value(
            sum(
                m.fs.feed.properties[0].conc_equiv_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(5, rel=1e-3)
        assert value(
            sum(
                m.fs.product.properties[0].conc_equiv_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(0.01057086, rel=1e-3)
        assert value(
            sum(
                m.fs.regen.properties[0].conc_equiv_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(471.99797, rel=1e-3)

        assert value(m.fs.ion_exchange.number_columns) == pytest.approx(4, rel=1e-6)
        assert value(m.fs.ion_exchange.bed_depth) == pytest.approx(1.7, rel=1e-6)
        assert value(m.fs.ion_exchange.t_breakthru) == pytest.approx(
            56759.75759, rel=1e-6
        )
        assert value(m.fs.ion_exchange.partition_ratio) == pytest.approx(
            235.99899, rel=1e-6
        )

        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            0.057352, abs=0.001
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.219, abs=0.01)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_optimization(self, ix_0D):
        m = ix_0D
        ixf.optimize_system(m)
        isinstance(m.fs.obj, Objective)
        assert m.fs.obj.expr == m.fs.costing.LCOW
        assert degrees_of_freedom(m) == 2

        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.05, abs=1e-6
        )
        assert value(m.fs.product.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.049901566, abs=1e-6
        )
        assert value(m.fs.regen.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.016666, abs=1e-6
        )
        ix = m.fs.ion_exchange
        num_col = math.ceil(
            ix.number_columns()
        )  # To eliminate fractional number of columns
        bed_depth = ix.bed_depth()
        ix.bed_depth.fix(bed_depth)
        ix.number_columns.fix(num_col)
        results = solver.solve(m)
        assert_optimal_termination(results)
        assert degrees_of_freedom(m) == 0
        assert value(m.fs.ion_exchange.number_columns) == 5
        assert value(m.fs.ion_exchange.bed_depth) == pytest.approx(1.24567, abs=1e-4)
        assert value(m.fs.ion_exchange.dimensionless_time) == pytest.approx(
            1.10981, abs=1e-4
        )
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            0.032097, abs=0.0001
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.160, abs=0.001)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.unit
    def test_main_fun(self):
        m = ixf.main()

        assert degrees_of_freedom(m) == 0
        assert value(m.fs.ion_exchange.number_columns) == 5
        assert value(m.fs.ion_exchange.bed_depth) == pytest.approx(1.24567, abs=1e-4)
        assert value(m.fs.ion_exchange.dimensionless_time) == pytest.approx(
            1.10981, abs=1e-4
        )
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            0.032097, abs=0.0001
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.160, abs=0.001)
