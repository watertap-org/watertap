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
            "resin_surf_per_vol": 4285.714285714285,
            "c_norm": {"Ca_2+": 0.4730788932825343},
            "bed_vol_tot": 11.999999999999998,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 1.4989640803696807,
            "col_height_to_diam_ratio": 2.327417344877002,
            "number_columns": 4,
            "t_breakthru": 56759.757564984655,
            "ebct": 240.0,
            "vel_bed": 0.007083333333333333,
            "service_flow_rate": 15,
            "N_Re": 4.958333333333333,
            "N_Sc": {"Ca_2+": 1086.9565217391305},
            "N_Sh": {"Ca_2+": 26.296358158587932},
            "N_Pe_particle": 0.1078279006415783,
            "N_Pe_bed": 261.86775870097586,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.685707070386448,
            "resin_unused_capacity": 1.3142929296135517,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 7079.96969562308},
            "num_transfer_units": 35.54838294744621,
            "dimensionless_time": 1,
            "partition_ratio": 235.99898985410272,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.4560927865572706e-05},
            "pressure_drop": 9.450141899999998,
            "bed_vol": 2.9999999999999996,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 1.357425724718769,
            "regen_tank_vol": 29.999999999999996,
            "bw_flow": 0.009803921568627449,
            "bed_expansion_frac": 0.46395000000000003,
            "rinse_flow": 0.04999999999999999,
            "t_cycle": 60359.757564984655,
            "bw_pump_power": 0.7984857204228053,
            "rinse_pump_power": 4.072277174156307,
            "bed_expansion_h": 0.788715,
            "main_pump_power": 4.072277174156307,
            "col_vol_per": 6.15655588235294,
            "col_vol_tot": 24.62622352941176,
            "t_contact": 120.0,
            "vel_inter": 0.014166666666666666,
            "bv_calc": 236.49898985410272,
            "lh": 0.0,
            "separation_factor": {"Ca_2+": 1.4285714285714286},
            "rate_coeff": {"Ca_2+": 0.00021159751754432275},
            "HTU": {"Ca_2+": 0.04782214714276131},
        }
        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 810841.861208609,
            "aggregate_fixed_operating_cost": 4099.257151281435,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 10.30046579345419,
            "aggregate_flow_NaCl": 2352713.2269726195,
            "aggregate_flow_costs": {
                "electricity": 6320.5718201793625,
                "NaCl": 214194.76894808252,
            },
            "total_capital_cost": 810841.861208609,
            "total_operating_cost": 226888.31967897541,
            "aggregate_direct_capital_cost": 405420.9306043045,
            "maintenance_labor_chemical_operating_cost": 24325.255836258268,
            "total_fixed_operating_cost": 28424.512987539703,
            "total_variable_operating_cost": 198463.8066914357,
            "total_annualized_cost": 307972.5057998363,
            "annual_water_production": 1419950.2910321492,
            "LCOW": 0.21688963884501467,
            "specific_energy_consumption": 0.05723052091619846,
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
            "resin_surf_per_vol": 4285.714285714285,
            "c_norm": {"Ca_2+": 0.9897743909187179},
            "bed_vol_tot": 11.999999999999998,
            "bed_depth": 1.6386683694597846,
            "bed_porosity": 0.5,
            "col_height": 3.3989285594706518,
            "col_diam": 1.3655737026706347,
            "col_height_to_diam_ratio": 2.489011430379571,
            "number_columns": 5,
            "t_breakthru": 133431.67855270393,
            "ebct": 240.0,
            "vel_bed": 0.006827784872749103,
            "service_flow_rate": 15,
            "N_Re": 4.779449410924372,
            "N_Sc": {"Ca_2+": 1086.9565217391305},
            "N_Sh": {"Ca_2+": 25.969879651128718},
            "N_Pe_particle": 0.1059427840517767,
            "N_Pe_bed": 248.00727028307858,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 2.978460143483587,
            "resin_unused_capacity": 0.02153985651641334,
            "langmuir": {"Ca_2+": 0.7},
            "mass_removed": {"Ca_2+": 12509.532602631058},
            "num_transfer_units": 35.10703730797482,
            "dimensionless_time": 1.332100914432498,
            "partition_ratio": 416.9844200877022,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.413184182719774e-05},
            "pressure_drop": 8.785890093905437,
            "bed_vol": 2.3999999999999995,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 1.2620120792068785,
            "regen_tank_vol": 29.999999999999996,
            "bw_flow": 0.010170860057646149,
            "bed_expansion_frac": 0.46395000000000003,
            "rinse_flow": 0.04999999999999999,
            "t_cycle": 137031.67855270393,
            "bw_pump_power": 0.7701448949203327,
            "rinse_pump_power": 3.7860362376206353,
            "bed_expansion_h": 0.7602601900108671,
            "main_pump_power": 3.7860362376206367,
            "col_vol_per": 4.978083848301045,
            "col_vol_tot": 24.890419241505228,
            "t_contact": 120.0,
            "vel_inter": 0.013655569745498206,
            "bv_calc": 555.9653273029331,
            "lh": 11.65907919299426,
            "separation_factor": {"Ca_2+": 1.4285714285714286},
            "rate_coeff": {"Ca_2+": 0.00020897046016651682},
            "HTU": {"Ca_2+": 0.04667635024524126},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 847088.1118126865,
            "aggregate_fixed_operating_cost": 3935.286865230177,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 9.604229449368482,
            "aggregate_flow_NaCl": 994870.9191908962,
            "aggregate_flow_costs": {
                "electricity": 5893.347274721491,
                "NaCl": 90574.63707273173,
            },
            "total_capital_cost": 847088.1118126865,
            "total_operating_cost": 116169.11613231867,
            "aggregate_direct_capital_cost": 423544.0559063433,
            "maintenance_labor_chemical_operating_cost": 25412.643354380594,
            "total_fixed_operating_cost": 29347.93021961077,
            "total_variable_operating_cost": 86821.18591270791,
            "total_annualized_cost": 200877.92731358734,
            "annual_water_production": 1419985.4904372608,
            "LCOW": 0.14146477458141515,
            "specific_energy_consumption": 0.05336083243675618,
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
            "c_norm": {"Ca_2+": 0.9897743909187181},
            "bed_vol_tot": 11.999999999999998,
            "bed_depth": 1.638668369459785,
            "number_columns": 5,
            "t_breakthru": 133431.67855270393,
            "ebct": 239.99999999999994,
            "resin_eq_capacity": 2.978460143483587,
            "mass_removed": {"Ca_2+": 12509.53260263106},
            "num_transfer_units": 35.107037307974814,
            "dimensionless_time": 1.3321009144324982,
            "partition_ratio": 416.98442008770223,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.413184182719774e-05},
            "bv_calc": 555.9653273029331,
            "lh": 11.659079192994266,
            "separation_factor": {"Ca_2+": 1.4285714285714286},
            "rate_coeff": {"Ca_2+": 0.00020897046016651682},
            "HTU": {"Ca_2+": 0.04667635024524128},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ion_exchange, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

        sys_cost_results = {
            "aggregate_capital_cost": 847088.1118126865,
            "aggregate_flow_electricity": 9.60422944936849,
            "aggregate_flow_NaCl": 994870.9191908962,
            "total_capital_cost": 847088.1118126865,
            "total_operating_cost": 116169.11613231867,
            "total_annualized_cost": 200877.92731358734,
            "annual_water_production": 1419985.4904372608,
            "LCOW": 0.14146477458141515,
            "specific_energy_consumption": 0.05336083243675621,
        }
        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
