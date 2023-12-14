#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
    Objective,
    assert_optimal_termination,
    Block,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
    IsothermType,
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

        # Test port
        assert isinstance(m.fs.feed.outlet, Port)
        assert isinstance(m.fs.ion_exchange.inlet, Port)
        assert isinstance(m.fs.ion_exchange.outlet, Port)
        assert isinstance(m.fs.product.inlet, Port)

        # # Test consting
        assert isinstance(m.fs.ion_exchange.costing, Block)
        assert isinstance(m.fs.ion_exchange.costing.capital_cost, Var)
        assert isinstance(m.fs.ion_exchange.costing.fixed_operating_cost, Var)

        var_str_list = [
            "total_capital_cost",
            "total_operating_cost",
        ]
        for var_str in var_str_list:
            var = getattr(m.fs.costing, var_str)
            assert isinstance(var, Var)

        # Test arcs
        arc_dict = {
            m.fs.feed_to_ix: (m.fs.feed.outlet, m.fs.ion_exchange.inlet),
            m.fs.ix_to_product: (m.fs.ion_exchange.outlet, m.fs.product.inlet),
        }
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

        # test configrations
        assert len(m.fs.ion_exchange.config) == 11
        assert not m.fs.ion_exchange.config.dynamic
        assert not m.fs.ion_exchange.config.has_holdup
        assert m.fs.ion_exchange.config.target_ion == "Ca_2+"
        assert m.fs.ion_exchange.ion_exchange_type == IonExchangeType.cation
        assert m.fs.ion_exchange.config.regenerant == RegenerantChem.NaCl
        assert m.fs.ion_exchange.config.isotherm == IsothermType.langmuir

        assert m.fs.ion_exchange.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        assert_units_consistent(m)

    @pytest.mark.component
    def test_specific_operating_conditions(self, ix_0D):

        m = ix_0D
        ixf.set_operating_conditions(m)
        ixf.initialize_system(m)
        assert degrees_of_freedom(m) == 0

        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)
        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.05, rel=1e-3
        )
        assert value(m.fs.product.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.049995010, rel=1e-3
        )
        assert value(
            sum(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(0.1, rel=1e-3)
        assert value(
            sum(
                m.fs.product.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(0.000211438, rel=1e-3)

        results_dict = {
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_surf_per_vol": 4285.714,
            "c_norm": {"Ca_2+": 0.47307},
            "bed_vol_tot": 11.999,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488,
            "col_diam": 1.498964,
            "col_height_to_diam_ratio": 2.3274,
            "number_columns": 4,
            "t_breakthru": 56759.7575,
            "ebct": 240.0,
            "vel_bed": 0.00708,
            "service_flow_rate": 15,
            "N_Re": 4.9583,
            "N_Sc": {"Ca_2+": 1086.956},
            "N_Sh": {"Ca_2+": 26.296358},
            "N_Pe_particle": 0.1078,
            "N_Pe_bed": 261.8677,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.6857,
            "resin_unused_capacity": 1.3142,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 7079.9696},
            "num_transfer_units": 35.54838,
            "dimensionless_time": 1,
            "partition_ratio": 235.998,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.4560e-05},
            "pressure_drop": 9.4501,
            "bed_vol": 3,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 1.3574,
            "regen_tank_vol": 29.9999,
            "bw_flow": 0.00980392,
            "bed_expansion_frac": 0.463950,
            "rinse_flow": 0.04999,
            "t_cycle": 60359.75756,
            "bw_pump_power": 0.798485720,
            "rinse_pump_power": 4.072277,
            "bed_expansion_h": 0.788715,
            "main_pump_power": 4.072277,
            "col_vol_per": 6.15655,
            "col_vol_tot": 24.6262,
            "t_contact": 120.0,
            "vel_inter": 0.01416,
            "bv_calc": 236.4989,
            "lh": 0.0,
            "separation_factor": {"Ca_2+": 1.428571},
            "rate_coeff": {"Ca_2+": 0.000211597},
            "HTU": {"Ca_2+": 0.047822},
        }
        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 810841.861208,
            "aggregate_fixed_operating_cost": 4099.25715,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 10.300,
            "aggregate_flow_NaCl": 2352713.226,
            "aggregate_flow_costs": {
                "electricity": 6320.57182,
                "NaCl": 214194.7689,
            },
            "total_capital_cost": 810841.8612,
            "total_operating_cost": 226888.3196,
            "aggregate_direct_capital_cost": 405420.9305,
            "maintenance_labor_chemical_operating_cost": 24325.25583,
            "total_fixed_operating_cost": 28424.5129,
            "total_variable_operating_cost": 198463.8066,
            "total_annualized_cost": 307972.50579,
            "annual_water_production": 1419950.29103,
            "LCOW": 0.21688,
            "specific_energy_consumption": 0.057231,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.component
    def test_optimization(self, ix_0D):
        m = ix_0D
        ixf.optimize_system(m)
        isinstance(m.fs.obj, Objective)
        assert m.fs.obj.expr == m.fs.costing.LCOW
        assert degrees_of_freedom(m) == 2

        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.05, rel=1e-3
        )
        assert value(m.fs.product.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            0.04999624, rel=1e-3
        )
        assert value(
            sum(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(0.1, rel=1e-3)
        assert value(
            sum(
                m.fs.product.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(0.025, rel=1e-3)
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

        results_dict = {
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_surf_per_vol": 4285.71428,
            "c_norm": {"Ca_2+": 0.989774},
            "bed_vol_tot": 12,
            "bed_depth": 1.63866,
            "bed_porosity": 0.5,
            "col_height": 3.398,
            "col_diam": 1.3655,
            "col_height_to_diam_ratio": 2.4890,
            "number_columns": 5,
            "t_breakthru": 133431.6785,
            "ebct": 240.0,
            "vel_bed": 0.0068277,
            "service_flow_rate": 15,
            "N_Re": 4.77944,
            "N_Sc": {"Ca_2+": 1086.9565},
            "N_Sh": {"Ca_2+": 25.96987},
            "N_Pe_particle": 0.105942,
            "N_Pe_bed": 248.007270,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 2.97846,
            "resin_unused_capacity": 0.0215398,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 12509.53260},
            "num_transfer_units": 35.10703,
            "dimensionless_time": 1.33210,
            "partition_ratio": 416.98442,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.41318e-05},
            "pressure_drop": 8.78589,
            "bed_vol": 2.39999,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 1.262012,
            "regen_tank_vol": 29.9999,
            "bw_flow": 0.010170,
            "bed_expansion_frac": 0.46395,
            "rinse_flow": 0.05,
            "t_cycle": 137031.6785,
            "bw_pump_power": 0.77014,
            "rinse_pump_power": 3.7860,
            "bed_expansion_h": 0.7602,
            "main_pump_power": 3.7860,
            "col_vol_per": 4.9780,
            "col_vol_tot": 24.8904,
            "t_contact": 120.0,
            "vel_inter": 0.01365,
            "bv_calc": 555.96532,
            "lh": 11.6591,
            "separation_factor": {"Ca_2+": 1.428571},
            "rate_coeff": {"Ca_2+": 0.00020897},
            "HTU": {"Ca_2+": 0.046676},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 847088.1118,
            "aggregate_fixed_operating_cost": 3935.286,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 9.60422,
            "aggregate_flow_NaCl": 994870.919,
            "aggregate_flow_costs": {
                "electricity": 5893.3472,
                "NaCl": 90574.6370,
            },
            "total_capital_cost": 847088.111,
            "total_operating_cost": 116169.116137,
            "aggregate_direct_capital_cost": 423544.0559,
            "maintenance_labor_chemical_operating_cost": 25412.6433,
            "total_fixed_operating_cost": 29347.93021,
            "total_variable_operating_cost": 86821.18591,
            "total_annualized_cost": 200877.9273,
            "annual_water_production": 1419985.4904,
            "LCOW": 0.14146,
            "specific_energy_consumption": 0.05336,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.unit
    def test_main_fun(self):
        m = ixf.main()

        assert degrees_of_freedom(m) == 0
        results_dict = {
            "c_norm": {"Ca_2+": 0.989774},
            "bed_vol_tot": 12,
            "bed_depth": 1.63866,
            "number_columns": 5,
            "t_breakthru": 133431.678,
            "ebct": 240,
            "resin_eq_capacity": 2.9784,
            "mass_removed": {"Ca_2+": 12509.532},
            "num_transfer_units": 35.1070,
            "dimensionless_time": 1.3321,
            "partition_ratio": 416.984,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.413e-05},
            "bv_calc": 555.96,
            "lh": 11.6590,
            "separation_factor": {"Ca_2+": 1.42857},
            "rate_coeff": {"Ca_2+": 0.00020897},
            "HTU": {"Ca_2+": 0.04667},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 847088.11,
            "aggregate_flow_electricity": 9.604,
            "aggregate_flow_NaCl": 994870.919,
            "total_capital_cost": 847088.1118,
            "total_operating_cost": 116169.116,
            "total_annualized_cost": 200877.927,
            "annual_water_production": 1419985.4904,
            "LCOW": 0.14145,
            "specific_energy_consumption": 0.05336,
        }
        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
