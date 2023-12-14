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
            "resin_surf_per_vol": 4285.714285714285,
            "c_norm": {"Ca_2+": 0.4919290557789296},
            "bed_vol_tot": 120.00000000000001,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 3.351785579537064,
            "col_height_to_diam_ratio": 1.0408526790314099,
            "number_columns": 8,
            "t_breakthru": 52360.64416318655,
            "ebct": 240.0,
            "vel_bed": 0.007083333333333333,
            "service_flow_rate": 15,
            "N_Re": 4.958333333333333,
            "N_Sc": {"Ca_2+": 1086.9565217391307},
            "N_Sh": {"Ca_2+": 26.29635815858793},
            "N_Pe_particle": 0.10782790064157834,
            "N_Pe_bed": 261.8677587009759,
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.5547810762853225,
            "resin_unused_capacity": 1.4452189237146778,
            "langmuir": {"Ca_2+": 0.9},
            "mass_removed": {"Ca_2+": 65300.80520398353},
            "num_transfer_units": 35.54838294744621,
            "dimensionless_time": 1,
            "partition_ratio": 217.66935067994396,
            "fluid_mass_transfer_coeff": {"Ca_2+": 3.4560927865572706e-05},
            "pressure_drop": 9.450141899999998,
            "bed_vol": 15.000000000000002,
            "t_rinse": 1200.0,
            "t_waste": 3600.0,
            "regen_pump_power": 13.574257247187694,
            "regen_tank_vol": 300.00000000000006,
            "bw_flow": 0.09803921568627452,
            "bed_expansion_frac": 0.46395000000000003,
            "rinse_flow": 0.5,
            "t_cycle": 55960.64416318655,
            "bw_pump_power": 7.984857204228056,
            "rinse_pump_power": 40.72277174156308,
            "bed_expansion_h": 0.788715,
            "main_pump_power": 40.722771741563086,
            "col_vol_per": 30.782779411764707,
            "col_vol_tot": 246.26223529411766,
            "t_contact": 120.0,
            "vel_inter": 0.014166666666666666,
            "bv_calc": 218.16935067994396,
            "lh": 0.0,
            "separation_factor": {"Ca_2+": 1.1111111111111112},
            "rate_coeff": {"Ca_2+": 0.00021159751754432275},
            "HTU": {"Ca_2+": 0.04782214714276131},
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

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 3993072.4698084667,
            "aggregate_fixed_operating_cost": 36893.31436153293,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 103.00465793454192,
            "aggregate_flow_NaCl": 22838957.969693653,
            "aggregate_flow_costs": {
                "electricity": 63205.71820179364,
                "NaCl": 2079295.2023431766,
            },
            "total_capital_cost": 3993072.4698084667,
            "total_operating_cost": 2084936.31694626,
            "LCOW": 0.17495285,
            "SEC": 0.0572305198,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 3993072.4698084667,
            "fixed_operating_cost": 36893.31436153293,
            "capital_cost_vessel": 101131.8811552528,
            "capital_cost_resin": 81985.14302562873,
            "capital_cost_regen_tank": 215778.261765144,
            "capital_cost_backwash_tank": 132704.75551115556,
            "operating_cost_hazardous": 0,
            "flow_mass_regen_soln": 22838957.969693653,
            "total_pumping_power": 103.00465793454192,
            "backwash_tank_vol": 174042.7639065449,
            "regeneration_tank_vol": 79251.61570744455,
            "direct_capital_cost": 1996536.2349042334,
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
            "aggregate_capital_cost": 5116213.693319044,
            "aggregate_fixed_operating_cost": 323306.0679291797,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 78.86531537032893,
            "aggregate_flow_NaOH": 279183.59700249793,
            "aggregate_flow_costs": {
                "electricity": 48393.334817541254,
                "NaOH": 555415.5212892868,
            },
            "total_capital_cost": 5116213.693319044,
            "total_operating_cost": 1020220.4492248963,
            "aggregate_direct_capital_cost": 2558106.8466595225,
            "maintenance_labor_chemical_operating_cost": 153486.4107995713,
            "total_fixed_operating_cost": 476792.478728751,
            "total_variable_operating_cost": 543427.9704961453,
            "total_annualized_cost": 1531841.8185568007,
            "LCOW": 0.10786919580206682,
            "SEC": 0.04381406413732131,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 5116213.693319045,
            "fixed_operating_cost": 323306.0679291797,
            "capital_cost_vessel": 75000.3204403362,
            "capital_cost_resin": 54924.687321091136,
            "capital_cost_regen_tank": 215778.261765144,
            "capital_cost_backwash_tank": 133603.4529501137,
            "operating_cost_hazardous": 276620.08370625216,
            "flow_mass_regen_soln": 279183.597002498,
            "total_pumping_power": 78.86531537032894,
            "backwash_tank_vol": 176401.06694050893,
            "regeneration_tank_vol": 79251.61570744455,
            "direct_capital_cost": 2558106.8466595225,
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
            "resin_surf_per_vol": 4444.444444444444,
            "c_norm": {"Cl_-": 0.25},
            "bed_vol_tot": 119.99999999999997,
            "bed_depth": 1.476,
            "bed_porosity": 0.5,
            "col_height": 3.1607902,
            "col_diam": 2.5435630784033805,
            "col_height_to_diam_ratio": 1.2426624001729332,
            "number_columns": 16,
            "t_breakthru": 4320000.0,
            "ebct": 240.0,
            "vel_bed": 0.00615,
            "service_flow_rate": 15,
            "N_Re": 4.15125,
            "N_Sc": {"Cl_-": 999.9999999999998},
            "N_Sh": {"Cl_-": 24.083093218519267},
            "N_Pe_particle": 0.09901383248136636,
            "N_Pe_bed": 216.5102470259211,
            "c_traps": {
                0: 0,
                1: 0.01,
                2: 0.06999999999999999,
                3: 0.13,
                4: 0.19,
                5: 0.25,
            },
            "tb_traps": {
                0: 0,
                1: 3344557.580567035,
                2: 3825939.149532492,
                3: 4034117.3864088706,
                4: 4188551.1118896296,
                5: 4320000.0,
            },
            "traps": {
                1: 0.0038710157182488833,
                2: 0.004457236749680154,
                3: 0.004818940668434703,
                4: 0.005719767610398483,
                5: 0.006694156338954019,
            },
            "c_norm_avg": {"Cl_-": 0.02556111708571624},
            "freundlich_n": 1.2,
            "mass_transfer_coeff": 0.15934630052514293,
            "bv": 18000,
            "bv_50": 20000,
            "pressure_drop": 7.151350934188799,
            "bed_vol": 7.499999999999998,
            "t_rinse": 1200.0,
            "t_waste": 1800.0,
            "bw_flow": 0.1129177958446251,
            "bed_expansion_frac": 0.46395000000000003,
            "rinse_flow": 0.4999999999999999,
            "t_cycle": 4321800.0,
            "bw_pump_power": 6.959523064801349,
            "rinse_pump_power": 30.816768130940375,
            "bed_expansion_h": 0.6847902,
            "main_pump_power": 30.81676813094038,
            "col_vol_per": 16.060925813008126,
            "col_vol_tot": 256.97481300813,
            "t_contact": 120.0,
            "vel_inter": 0.0123,
            "c_breakthru": {"Cl_-": 2.500000001615544e-07},
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
        ix.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 4485122.812315232,
            "aggregate_fixed_operating_cost": 6419597.45408913,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 30.81676813094038,
            "aggregate_flow_costs": {"electricity": 18909.78526050764},
            "total_capital_cost": 4485122.812315232,
            "total_operating_cost": 6571169.945193044,
            "capital_recovery_factor": 0.1,
            "HCl_cost": 0.4594594594594595,
            "NaOH_cost": 1.9666666666666666,
            "MeOH_cost": 3.395,
            "NaCl_cost": 0.09,
            "aggregate_direct_capital_cost": 2242561.406157616,
            "maintenance_labor_chemical_operating_cost": 134553.68436945695,
            "total_fixed_operating_cost": 6554151.138458587,
            "total_variable_operating_cost": 17018.806734456877,
            "total_annualized_cost": 7019682.226424567,
            "LCOW": 0.4943117934094987,
            "SEC": 0.017120426756094133,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 4485122.812315232,
            "fixed_operating_cost": 6419597.45408913,
            "capital_cost_vessel": 75000.32044033613,
            "capital_cost_resin": 54924.68732109111,
            "capital_cost_regen_tank": 0,
            "capital_cost_backwash_tank": 33836.27421335332,
            "operating_cost_hazardous": 0,
            "flow_mass_regen_soln": 0,
            "total_pumping_power": 30.81676813094038,
            "flow_vol_resin": 876.5999999999999,
            "single_use_resin_replacement_cost": 6419597.45408913,
            "direct_capital_cost": 2242561.406157616,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.ix.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
