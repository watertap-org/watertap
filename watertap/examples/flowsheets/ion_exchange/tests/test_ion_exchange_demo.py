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
            "c_norm": {"Ca_2+": 0.473078},
            "bed_vol_tot": 11.999,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 1.498964,
            "col_height_to_diam_ratio": 2.327417,
            "number_columns": 4.0,
            "t_breakthru": 56759.757,
            "ebct": 240.0,
            "vel_bed": 0.007083333,
            "service_flow_rate": 15.0,
            "N_Re": 4.958333,
            "N_Sc": {"Ca_2+": 1086.956},
            "N_Sh": {"Ca_2+": 26.296},
            "N_Pe_particle": 0.107827,
            "N_Pe_bed": 261.867,
            "resin_max_capacity": 3.0,
            "resin_eq_capacity": 1.685707,
            "resin_unused_capacity": 1.314292,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 7079.969},
            "num_transfer_units": 35.548,
            "dimensionless_time": 1.0,
            "partition_ratio": 235.998,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.456e-05},
            "pressure_drop": 9.450141,
            "bed_vol": 2.999999,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "t_cycle": 60359.757,
            "regen_pump_power": 0.040480054,
            "regen_tank_vol": 29.999,
            "bw_flow": 0.009803921,
            "bed_expansion_frac": 0.46395,
            "rinse_flow": 0.049999999,
            "bw_pump_power": 0.007937265,
            "rinse_pump_power": 0.080960109,
            "bed_expansion_h": 0.788715,
            "main_pump_power": 3.829396,
            "col_vol_per": 6.156555,
            "col_vol_tot": 24.626,
            "t_contact": 120.0,
            "vel_inter": 0.014166666,
            "bv_calc": 236.498,
            "lh": 0.0,
            "separation_factor": {"Ca_2+": 1.428571},
            "rate_coeff": {"Ca_2+": 0.000211597},
            "HTU": {"Ca_2+": 0.047822147},
        }
        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 810841.861,
            "aggregate_fixed_operating_cost": 4099.257,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 3.958774,
            "aggregate_flow_NaCl": 2352713.226,
            "aggregate_flow_costs": {"electricity": 2429.183, "NaCl": 214194.768},
            "total_capital_cost": 810841.861,
            "total_operating_cost": 223386.069,
            "capital_recovery_factor": 0.1,
            "HCl_cost": 0.459459,
            "NaOH_cost": 1.966666,
            "MeOH_cost": 3.395,
            "NaCl_cost": 0.09,
            "aggregate_direct_capital_cost": 405420.93,
            "maintenance_labor_chemical_operating_cost": 24325.255,
            "total_fixed_operating_cost": 28424.512,
            "total_variable_operating_cost": 194961.556,
            "total_annualized_cost": 304470.255,
            "annual_water_production": 1419950.291,
            "LCOW": 0.214423,
            "specific_energy_consumption": 0.021995385,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.component
    @pytest.mark.requires_idaes_solver
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
            "resin_surf_per_vol": 4285.714,
            "c_norm": {"Ca_2+": 0.99398},
            "bed_vol_tot": 11.999,
            "bed_depth": 2.211,
            "bed_porosity": 0.5,
            "col_height": 4.237165,
            "col_diam": 1.17555,
            "col_height_to_diam_ratio": 3.604409,
            "number_columns": 5.0,
            "t_breakthru": 133829.353,
            "ebct": 239.999,
            "vel_bed": 0.009213559,
            "service_flow_rate": 15.0,
            "N_Re": 6.449491,
            "N_Sc": {"Ca_2+": 1086.956},
            "N_Sh": {"Ca_2+": 28.755},
            "N_Pe_particle": 0.122332,
            "N_Pe_bed": 386.44,
            "resin_max_capacity": 3.0,
            "resin_eq_capacity": 2.987337,
            "resin_unused_capacity": 0.01266295,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 12546.815},
            "num_transfer_units": 38.872,
            "dimensionless_time": 1.332104,
            "partition_ratio": 418.227,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.7792e-05},
            "pressure_drop": 16.049,
            "bed_vol": 2.399999,
            "t_rinse": 1199.999,
            "t_waste": 3600.0,
            "t_cycle": 137429.353,
            "regen_pump_power": 0.030195009,
            "regen_tank_vol": 29.999,
            "bw_flow": 0.0075372,
            "bed_expansion_frac": 0.46395,
            "rinse_flow": 0.049999999,
            "bw_pump_power": 0.004551716,
            "rinse_pump_power": 0.060390018,
            "bed_expansion_h": 1.025911,
            "main_pump_power": 6.734964,
            "col_vol_per": 4.598836,
            "col_vol_tot": 22.994,
            "t_contact": 119.999,
            "vel_inter": 0.018427119,
            "bv_calc": 557.622,
            "lh": 12.909,
            "separation_factor": {"Ca_2+": 1.428571},
            "rate_coeff": {"Ca_2+": 0.000231384},
            "HTU": {"Ca_2+": 0.05688455},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 827291.05,
            "aggregate_fixed_operating_cost": 3935.28,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 6.8301,
            "aggregate_flow_NaCl": 991992.09,
            "aggregate_flow_costs": {"electricity": 4191.08, "NaCl": 90312.54},
            "total_capital_cost": 827291.05,
            "total_operating_cost": 113807.28,
            "aggregate_direct_capital_cost": 413645.52,
            "maintenance_labor_chemical_operating_cost": 24818.73,
            "total_fixed_operating_cost": 28754.0186,
            "total_variable_operating_cost": 85053.26,
            "total_annualized_cost": 196536.39,
            "annual_water_production": 1419985.49,
            "LCOW": 0.1384073,
            "specific_energy_consumption": 0.0379478,
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
            "resin_surf_per_vol": 4285.7142,
            "c_norm": {"Ca_2+": 0.993980},
            "bed_vol_tot": 11.9999,
            "bed_depth": 2.2112,
            "bed_porosity": 0.5,
            "col_height": 4.2371,
            "col_diam": 1.1755,
            "col_height_to_diam_ratio": 3.6044,
            "number_columns": 5.0,
            "t_breakthru": 133829.3,
            "ebct": 239.9999,
            "vel_bed": 0.00921355996,
            "service_flow_rate": 15.0,
            "N_Re": 6.44949197,
            "N_Sc": {"Ca_2+": 1086.9565},
            "N_Sh": {"Ca_2+": 28.7554},
            "N_Pe_particle": 0.122332,
            "N_Pe_bed": 386.4407,
            "resin_max_capacity": 3.0,
            "resin_eq_capacity": 2.98733704,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 12546.8156},
            "num_transfer_units": 38.8726,
            "dimensionless_time": 1.33210447,
            "partition_ratio": 418.2271,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.7792874e-05},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 827291.06,
            "aggregate_fixed_operating_cost": 3935.28,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 6.8301,
            "aggregate_flow_NaCl": 991992.0893,
            "aggregate_flow_costs": {"electricity": 4191.0865, "NaCl": 90312.54},
            "total_capital_cost": 827291.058,
            "total_operating_cost": 113807.2859,
            "aggregate_direct_capital_cost": 413645.529,
            "maintenance_labor_chemical_operating_cost": 24818.73,
            "total_fixed_operating_cost": 28754.0186,
            "total_variable_operating_cost": 85053.26,
            "total_annualized_cost": 196536.39,
            "annual_water_production": 1419985.49,
            "LCOW": 0.1384073,
            "specific_energy_consumption": 0.03794785,
        }
        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
