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
    Param,
    assert_optimal_termination,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    EnergyBalanceType,
    MomentumBalanceType,
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
)
import idaes.logger as idaeslog

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
    IsothermType,
)
from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestIonExchangeLangmuir:
    @pytest.fixture(scope="class")
    def IX_lang(self):
        target_ion = "Ca_2+"
        ion_props = {
            "solute_list": [target_ion],
            "diffusivity_data": {("Liq", target_ion): 9.2e-10},
            "mw_data": {"H2O": 0.018, target_ion: 0.04},
            "charge": {target_ion: 2},
        }
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**ion_props)
        ix_config = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
        }
        ix = m.fs.ix = IonExchange0D(**ix_config)
        ix.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 0.1,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
        ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
        ix.regeneration_stream[0].flow_mass_phase_comp[...]

        ix.service_flow_rate.fix(15)
        ix.langmuir[target_ion].fix(0.9)
        ix.resin_max_capacity.fix(3)
        ix.bed_depth.fix(1.7)
        ix.dimensionless_time.fix()
        ix.number_columns.fix(8)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()

        return m

    @pytest.mark.unit
    def test_config(self, IX_lang):
        m = IX_lang

        assert len(m.fs.ix.config) == 11

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert not m.fs.ix.config.hazardous_waste
        assert m.fs.ix.config.regenerant is RegenerantChem.NaCl
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.cation
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.langmuir
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_default_build(self, IX_lang):
        m = IX_lang
        ix = m.fs.ix
        # test ports and variables
        port_lst = ["inlet", "outlet", "regen"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(ix, port_str)
            port = getattr(ix, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects
        ix_params = [
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "bw_rate",
            "distributor_h",
            "number_columns_redund",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "rinse_bv",
            "service_to_regen_flow_ratio",
            "t_bw",
            "t_regen",
            "underdrain_h",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "N_Pe_bed",
            "N_Pe_particle",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "bed_depth",
            "bed_porosity",
            "bed_vol_tot",
            "c_norm",
            "col_diam",
            "col_height",
            "col_height_to_diam_ratio",
            "dimensionless_time",
            "ebct",
            "fluid_mass_transfer_coeff",
            "langmuir",
            "mass_removed",
            "num_transfer_units",
            "number_columns",
            "partition_ratio",
            "resin_bulk_dens",
            "resin_diam",
            "resin_eq_capacity",
            "resin_max_capacity",
            "resin_surf_per_vol",
            "resin_unused_capacity",
            "service_flow_rate",
            "t_breakthru",
            "vel_bed",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 70
        assert number_total_constraints(m) == 44
        assert number_unused_variables(m) == 10

    @pytest.mark.unit
    def test_dof(self, IX_lang):
        m = IX_lang
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_lang):
        m = IX_lang
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 10, index=("Liq", "Ca_2+")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, IX_lang):
        m = IX_lang
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, IX_lang):
        m = IX_lang
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_mass_balance(self, IX_lang):
        m = IX_lang

        ix = m.fs.ix
        target = ix.config.target_ion
        pf = ix.process_flow
        prop_in = pf.properties_in[0]
        prop_out = pf.properties_out[0]
        regen = ix.regeneration_stream[0]

        assert value(prop_in.flow_mass_phase_comp["Liq", target]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", target])
            + value(regen.flow_mass_phase_comp["Liq", target]),
            rel=1e-3,
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", "H2O"]),
            rel=1e-3,
        )

        assert -1 * value(pf.mass_transfer_term[0, "Liq", target]) == pytest.approx(
            value(regen.flow_mol_phase_comp["Liq", target]), rel=1e-3
        )

    @pytest.mark.component
    def test_solution(self, IX_lang):
        m = IX_lang

        # results for all Var and Expressions on unit model
        results_dict = {
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "resin_surf_per_vol": 4285.714,
            "c_norm": {"Ca_2+": 0.4919},
            "bed_vol_tot": 120.0,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488,
            "col_diam": 3.351,
            "col_height_to_diam_ratio": 1.0408,
            "number_columns": 8,
            "t_breakthru": 52360.644,
            "ebct": 240.0,
            "vel_bed": 0.007083,
            "service_flow_rate": 15,
            "N_Re": 4.958,
            "N_Sc": {"Ca_2+": 1086.9565},
            "N_Sh": {"Ca_2+": 26.296},
            "N_Pe_particle": 0.10782,
            "N_Pe_bed": 261.8677,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.554,
            "resin_unused_capacity": 1.445,
            "langmuir": {"Ca_2+": 0.9},
            "mass_removed": {"Ca_2+": 65300.8052},
            "num_transfer_units": 35.5483,
            "dimensionless_time": 1,
            "partition_ratio": 217.669,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.45609e-05},
            "pressure_drop": 9.450,
            "bed_vol": 15.0,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 0.43662,
            "regen_tank_vol": 300.0,
            "bw_flow": 0.09803,
            "bed_expansion_frac": 0.46395,
            "rinse_flow": 0.5,
            "t_cycle": 55960.6441,
            "bw_pump_power": 0.08561221,
            "rinse_pump_power": 0.873244,
            "bed_expansion_h": 0.788,
            "main_pump_power": 38.103,
            "col_vol_per": 30.7827,
            "col_vol_tot": 246.2622,
            "t_contact": 120.0,
            "vel_inter": 0.01416,
            "bv_calc": 218.169,
            "lh": 0.0,
            "separation_factor": {"Ca_2+": 1.11111},
            "rate_coeff": {"Ca_2+": 0.00021159},
            "HTU": {"Ca_2+": 0.04782},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ix, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

    @pytest.mark.component
    def test_costing(self, IX_lang):
        m = IX_lang
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        check_dof(m, fail_flag=True)

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 3993072.469,
            "aggregate_fixed_operating_cost": 36893.314,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 39.4985,
            "aggregate_flow_NaCl": 22838957.969,
            "aggregate_flow_costs": {
                "electricity": 24237.080,
                "NaCl": 2079295.202,
            },
            "total_capital_cost": 3993072.4698,
            "total_operating_cost": 2049864.542,
            "LCOW": 0.1724829,
            "SEC": 0.021945,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 3993072.4698,
            "fixed_operating_cost": 36893.314,
            "capital_cost_vessel": 101131.881,
            "capital_cost_resin": 81985.1430,
            "capital_cost_regen_tank": 215778.261,
            "capital_cost_backwash_tank": 132704.7555,
            "operating_cost_hazardous": 0,
            "flow_mass_regen_soln": 22838957.969,
            "total_pumping_power": 39.4985,
            "backwash_tank_vol": 174042.7639,
            "regeneration_tank_vol": 79251.61570,
            "direct_capital_cost": 1996536.234,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.ix.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


class TestIonExchangeFreundlich:
    @pytest.fixture(scope="class")
    def IX_fr(self):

        target_ion = "Cl_-"

        ion_props = {
            "solute_list": [target_ion],
            "diffusivity_data": {("Liq", target_ion): 1e-9},
            "mw_data": {"H2O": 0.018, target_ion: 35.45e-3},
            "charge": {target_ion: -1},
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**ion_props)
        ix_config = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
            "isotherm": "freundlich",
            "regenerant": "NaOH",
            "hazardous_waste": True,
        }
        m.fs.ix = ix = IonExchange0D(**ix_config)

        ix.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 1e-6,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
        ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
        ix.regeneration_stream[0].flow_mass_phase_comp[...]

        ix.freundlich_n.fix(1.2)
        ix.bv_50.fix(20000)
        ix.bv.fix(18000)
        ix.resin_bulk_dens.fix(0.72)
        ix.bed_porosity.fix()
        ix.vel_bed.fix(6.15e-3)
        ix.resin_diam.fix(6.75e-4)
        ix.number_columns.fix(16)
        ix.service_flow_rate.fix(15)
        ix.c_norm.fix(0.25)

        return m

    @pytest.mark.unit
    def test_config(self, IX_fr):
        m = IX_fr

        assert len(m.fs.ix.config) == 11

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert m.fs.ix.config.hazardous_waste
        assert m.fs.ix.config.regenerant is RegenerantChem.NaOH
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.anion
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.freundlich
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_fr_build(self, IX_fr):
        m = IX_fr
        ix = m.fs.ix
        # test ports and variables

        port_lst = ["inlet", "outlet", "regen"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(ix, port_str)
            port = getattr(ix, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects
        ix_params = [
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "bw_rate",
            "c_trap_min",
            "distributor_h",
            "number_columns_redund",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "rinse_bv",
            "service_to_regen_flow_ratio",
            "t_bw",
            "t_regen",
            "underdrain_h",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "N_Pe_bed",
            "N_Pe_particle",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "bed_depth",
            "bed_porosity",
            "bed_vol_tot",
            "bv",
            "bv_50",
            "c_norm",
            "c_norm_avg",
            "c_traps",
            "col_diam",
            "col_height",
            "col_height_to_diam_ratio",
            "ebct",
            "freundlich_n",
            "mass_transfer_coeff",
            "number_columns",
            "resin_bulk_dens",
            "resin_diam",
            "resin_surf_per_vol",
            "service_flow_rate",
            "t_breakthru",
            "tb_traps",
            "traps",
            "vel_bed",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 80
        assert number_total_constraints(m) == 51
        assert number_unused_variables(m) == 11

    @pytest.mark.unit
    def test_dof(self, IX_fr):
        m = IX_fr
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_fr):
        m = IX_fr
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e6, index=("Liq", "Cl_-")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, IX_fr):
        m = IX_fr
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, IX_fr):
        m = IX_fr
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_mass_balance(self, IX_fr):
        m = IX_fr

        ix = m.fs.ix
        target = ix.config.target_ion
        pf = ix.process_flow
        prop_in = pf.properties_in[0]
        prop_out = pf.properties_out[0]
        regen = ix.regeneration_stream[0]

        assert value(prop_in.flow_mass_phase_comp["Liq", target]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", target])
            + value(regen.flow_mass_phase_comp["Liq", target]),
            rel=1e-3,
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", "H2O"]),
            rel=1e-3,
        )

        assert -1 * value(pf.mass_transfer_term[0, "Liq", target]) == pytest.approx(
            value(regen.flow_mol_phase_comp["Liq", target]), rel=1e-3
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, IX_fr):
        m = IX_fr
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )
        ix.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 5116213.693,
            "aggregate_fixed_operating_cost": 323306.067,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 30.805,
            "aggregate_flow_NaOH": 279183.597,
            "aggregate_flow_costs": {
                "electricity": 18902.5053,
                "NaOH": 555415.521,
            },
            "total_capital_cost": 5116213.693,
            "total_operating_cost": 993678.703,
            "aggregate_direct_capital_cost": 2558106.846,
            "maintenance_labor_chemical_operating_cost": 153486.410,
            "total_fixed_operating_cost": 476792.478,
            "total_variable_operating_cost": 516886.223,
            "total_annualized_cost": 1505300.072,
            "LCOW": 0.1060001,
            "SEC": 0.017113,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 5116213.6933,
            "fixed_operating_cost": 323306.0679,
            "capital_cost_vessel": 75000.3204,
            "capital_cost_resin": 54924.687,
            "capital_cost_regen_tank": 215778.261,
            "capital_cost_backwash_tank": 133603.4529,
            "operating_cost_hazardous": 276620.0837,
            "flow_mass_regen_soln": 279183.597,
            "total_pumping_power": 30.8049,
            "backwash_tank_vol": 176401.0669,
            "regeneration_tank_vol": 79251.615,
            "direct_capital_cost": 2558106.846,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.ix.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


class TestIonExchangeInert:
    @pytest.fixture(scope="class")
    def IX_inert(self):

        target_ion = "Cl_-"
        inert = "Ca_2+"

        ion_props = {
            "solute_list": [target_ion, inert],
            "diffusivity_data": {("Liq", target_ion): 1e-9, ("Liq", inert): 9.2e-10},
            "mw_data": {"H2O": 0.018, target_ion: 35.45e-3, inert: 40e-3},
            "charge": {target_ion: -1, inert: 2},
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**ion_props)
        ix_config = {
            "property_package": m.fs.properties,
            "target_ion": target_ion,
            "regenerant": "single_use",
            "isotherm": "freundlich",
        }
        m.fs.ix = ix = IonExchange0D(**ix_config)

        ix.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 1e-6,
                ("conc_mass_phase_comp", ("Liq", inert)): 1e-7,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
        ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
        ix.regeneration_stream[0].flow_mass_phase_comp[...]

        ix.freundlich_n.fix(1.2)
        ix.bv_50.fix(20000)
        ix.bv.fix(18000)
        ix.resin_bulk_dens.fix(0.72)
        ix.bed_porosity.fix()
        ix.vel_bed.fix(6.15e-3)
        ix.resin_diam.fix(6.75e-4)
        ix.number_columns.fix(16)
        ix.service_flow_rate.fix(15)
        ix.c_norm.fix(0.25)

        return m

    @pytest.mark.unit
    def test_config(self, IX_inert):
        m = IX_inert

        assert len(m.fs.ix.config) == 11

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert not m.fs.ix.config.hazardous_waste
        assert m.fs.ix.config.regenerant is RegenerantChem.single_use
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.anion
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.freundlich
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_inert_build(self, IX_inert):
        m = IX_inert
        ix = m.fs.ix
        # test ports and variables

        port_lst = ["inlet", "outlet", "regen"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(ix, port_str)
            port = getattr(ix, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects
        ix_params = [
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "bw_rate",
            "c_trap_min",
            "distributor_h",
            "number_columns_redund",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "rinse_bv",
            "service_to_regen_flow_ratio",
            "t_bw",
            "t_regen",
            "underdrain_h",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "N_Pe_bed",
            "N_Pe_particle",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "bed_depth",
            "bed_porosity",
            "bed_vol_tot",
            "bv",
            "bv_50",
            "c_norm",
            "c_norm_avg",
            "c_traps",
            "col_diam",
            "col_height",
            "col_height_to_diam_ratio",
            "ebct",
            "freundlich_n",
            "mass_transfer_coeff",
            "number_columns",
            "resin_bulk_dens",
            "resin_diam",
            "resin_surf_per_vol",
            "service_flow_rate",
            "t_breakthru",
            "tb_traps",
            "traps",
            "vel_bed",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 90
        assert number_total_constraints(m) == 57
        assert number_unused_variables(m) == 12

    @pytest.mark.unit
    def test_dof(self, IX_inert):
        m = IX_inert
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_inert):
        m = IX_inert
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e5, index=("Liq", "Cl_-")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, IX_inert):
        m = IX_inert
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, IX_inert):
        m = IX_inert
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_mass_balance(self, IX_inert):
        m = IX_inert

        ix = m.fs.ix
        target = ix.config.target_ion
        inert = "Ca_2+"
        pf = ix.process_flow
        prop_in = pf.properties_in[0]
        prop_out = pf.properties_out[0]
        regen = ix.regeneration_stream[0]

        assert value(prop_in.flow_mass_phase_comp["Liq", target]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", target])
            + value(regen.flow_mass_phase_comp["Liq", target]),
            rel=1e-3,
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", "H2O"]),
            rel=1e-3,
        )

        assert -1 * value(pf.mass_transfer_term[0, "Liq", target]) == pytest.approx(
            value(regen.flow_mol_phase_comp["Liq", target]), rel=1e-3
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", inert]) == pytest.approx(
            value(prop_out.flow_mass_phase_comp["Liq", inert]),
            rel=1e-3,
        )

        assert value(pf.mass_transfer_term[0, "Liq", inert]) == 0

    @pytest.mark.component
    def test_solution(self, IX_inert):
        m = IX_inert

        results_dict = {
            "resin_diam": 0.000675,
            "resin_bulk_dens": 0.72,
            "resin_surf_per_vol": 4444.44,
            "c_norm": {"Cl_-": 0.25},
            "bed_vol_tot": 120,
            "bed_depth": 1.476,
            "bed_porosity": 0.5,
            "col_height": 3.1607902,
            "col_diam": 2.5435,
            "col_height_to_diam_ratio": 1.2426,
            "number_columns": 16,
            "t_breakthru": 4320000.0,
            "ebct": 240.0,
            "vel_bed": 0.00615,
            "service_flow_rate": 15,
            "N_Re": 4.15125,
            "N_Sc": {"Cl_-": 1000},
            "N_Sh": {"Cl_-": 24.0830},
            "N_Pe_particle": 0.0990,
            "N_Pe_bed": 216.5102,
            "c_traps": {
                0: 0,
                1: 0.01,
                2: 0.07,
                3: 0.13,
                4: 0.19,
                5: 0.25,
            },
            "tb_traps": {
                0: 0,
                1: 3344557.580,
                2: 3825939.149,
                3: 4034117.386,
                4: 4188551.111,
                5: 4320000.0,
            },
            "traps": {
                1: 0.00387,
                2: 0.00445724,
                3: 0.00481894,
                4: 0.00571976,
                5: 0.00669,
            },
            "c_norm_avg": {"Cl_-": 0.02556},
            "freundlich_n": 1.2,
            "mass_transfer_coeff": 0.1593,
            "bv": 18000,
            "bv_50": 20000,
            "pressure_drop": 7.1513,
            "bed_vol": 7.5,
            "t_rinse": 1200.0,
            "t_waste": 1800.0,
            "bw_flow": 0.1129,
            "bed_expansion_frac": 0.4639,
            "rinse_flow": 0.5,
            "t_cycle": 4321800.0,
            "bw_pump_power": 0.0009661,
            "rinse_pump_power": 0.0085566,
            "bed_expansion_h": 0.68479,
            "main_pump_power": 30.816,
            "col_vol_per": 16.060,
            "col_vol_tot": 256.97,
            "t_contact": 120.0,
            "vel_inter": 0.0123,
            "c_breakthru": {"Cl_-": 2.50e-07},
        }

        for v, r in results_dict.items():
            ixv = getattr(m.fs.ix, v)
            if ixv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(ixv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(ixv)

    @pytest.mark.component
    def test_costing(self, IX_inert):

        m = IX_inert
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )
        check_dof(m, fail_flag=True)

        ix.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 4684657.16,
            "aggregate_fixed_operating_cost": 6419597.45,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 30.813,
            "aggregate_flow_costs": {"electricity": 18907.75},
            "total_capital_cost": 4684657.16,
            "total_operating_cost": 6577154.15,
            "aggregate_direct_capital_cost": 2342328.58,
            "maintenance_labor_chemical_operating_cost": 140539.72,
            "total_fixed_operating_cost": 6560137.16,
            "total_variable_operating_cost": 17016.98,
            "total_annualized_cost": 7045619.86,
            "LCOW": 0.49613827,
            "SEC": 0.017118,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 4684657.16,
            "fixed_operating_cost": 6419597.45,
            "capital_cost_vessel": 75000.32,
            "capital_cost_resin": 54924.68,
            "capital_cost_regen_tank": 0,
            "capital_cost_backwash_tank": 133603.45,
            "operating_cost_hazardous": 0,
            "flow_mass_regen_soln": 0,
            "total_pumping_power": 30.813,
            "flow_vol_resin": 876.599,
            "single_use_resin_replacement_cost": 6419597.454,
            "direct_capital_cost": 2342328.58,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.ix.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
