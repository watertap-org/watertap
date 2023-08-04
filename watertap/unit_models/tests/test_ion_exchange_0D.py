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

        ix.service_flow_rate.fix(15)
        ix.langmuir[target_ion].fix(0.9)
        ix.resin_max_capacity.fix(3)
        ix.bed_depth.fix(1.7)
        ix.dimensionless_time.fix()
        ix.number_columns.fix(8)
        ix.resin_diam.fix()
        ix.resin_bulk_dens.fix()
        ix.bed_porosity.fix()
        ix.regen_dose.fix()

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
            "underdrain_h",
            "distributor_h",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "t_regen",
            "rinse_bv",
            "bw_rate",
            "t_bw",
            "service_to_regen_flow_ratio",
            "number_columns_redund",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_diam",
            "resin_bulk_dens",
            "resin_surf_per_vol",
            "regen_dose",
            "c_norm",
            "col_height_to_diam_ratio",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "number_columns",
            "t_breakthru",
            "t_contact",
            "ebct",
            "vel_bed",
            "vel_inter",
            "service_flow_rate",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "N_Pe_particle",
            "N_Pe_bed",
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "langmuir",
            "mass_removed",
            "num_transfer_units",
            "dimensionless_time",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 69
        assert number_total_constraints(m) == 42
        assert number_unused_variables(m) == 12

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

    @pytest.mark.component
    def test_solution(self, IX_lang):
        m = IX_lang
        ix = m.fs.ix
        target_ion = ix.config.target_ion

        results_dict = {
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.5547810762853227,
            "resin_unused_capacity": 1.4452189237146773,
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "langmuir": 0.9,
            "num_transfer_units": 35.54838294744622,
            "dimensionless_time": 1,
            "resin_surf_per_vol": 4285.714285714286,
            "col_height_to_diam_ratio": 1.0408526790314099,
            "bed_vol_tot": 120.00000000000003,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 3.3517855795370646,
            "number_columns": 8,
            "partition_ratio": 217.66935067994518,
            "fluid_mass_transfer_coeff": 3.456092786557271e-05,
            "t_breakthru": 52360.64416318684,
            "t_contact": 120.0,
            "mass_removed": 65300.80520398353,
            "vel_bed": 0.007083333333333333,
            "vel_inter": 0.014166666666666666,
            "service_flow_rate": 15,
            "N_Re": 4.958333333333333,
            "N_Sc": 1086.9565217391305,
            "N_Sh": 26.29635815858793,
            "N_Pe_particle": 0.10782790064157834,
            "N_Pe_bed": 261.86775870097597,
            "c_norm": 0.4919290557789296,
            "regen_dose": 300,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=1e-3) == value(var)

    @pytest.mark.component
    def test_costing(self, IX_lang):
        m = IX_lang
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(8894349.86900, rel=1e-3) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(2498819.7327, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(17788699.7380, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.30125629, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.0572452, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )


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

        # c0 = pyunits.convert(c0, to_units=pyunits.kg / pyunits.m**3)
        ix.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 1e-6,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        ix.freundlich_n.fix(1.2)
        ix.bv_50.fix(20000)
        ix.bv.fix(18000)
        ix.resin_bulk_dens.fix(0.72)
        ix.regen_dose.fix()
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
    def test_default_build(self, IX_fr):
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
            "underdrain_h",
            "distributor_h",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "t_regen",
            "rinse_bv",
            "bw_rate",
            "t_bw",
            "service_to_regen_flow_ratio",
            "number_columns_redund",
            "c_trap_min",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_diam",
            "resin_bulk_dens",
            "resin_surf_per_vol",
            "regen_dose",
            "c_norm",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "col_height_to_diam_ratio",
            "number_columns",
            "t_breakthru",
            "t_contact",
            "ebct",
            "vel_bed",
            "vel_inter",
            "service_flow_rate",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "N_Pe_particle",
            "N_Pe_bed",
            "c_traps",
            "tb_traps",
            "traps",
            "c_norm_avg",
            "c_breakthru",
            "freundlich_n",
            "mass_transfer_coeff",
            "bv",
            "bv_50",
            "bed_capacity_param",
            "kinetic_param",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 82
        assert number_total_constraints(m) == 52
        assert number_unused_variables(m) == 13

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
    def test_solution(self, IX_fr):
        m = IX_fr
        ix = m.fs.ix
        target_ion = ix.config.target_ion
        results_dict = {
            "resin_diam": 0.0006749999999999999,
            "resin_bulk_dens": 0.72,
            "resin_surf_per_vol": 4444.444444444445,
            "regen_dose": 300,
            "c_norm": {"Cl_-": 0.25},
            "bed_vol_tot": 120.00000000000001,
            "bed_depth": 1.476,
            "bed_porosity": 0.5,
            "col_height": 3.1607902,
            "col_diam": 2.5435630784033814,
            "col_height_to_diam_ratio": 1.242662400172933,
            "number_columns": 16,
            "t_breakthru": 4320000.0,
            "t_contact": 120.00000000000001,
            "ebct": 240.00000000000003,
            "vel_bed": 0.006149999999999999,
            "vel_inter": 0.012299999999999998,
            "service_flow_rate": 15,
            "N_Re": 4.151249999999999,
            "N_Sc": {"Cl_-": 999.9999999999998},
            "N_Sh": {"Cl_-": 24.083093218519274},
            "N_Pe_particle": 0.09901383248136636,
            "N_Pe_bed": 216.51024702592113,
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
                1: 3344557.580567041,
                2: 3825939.1495324974,
                3: 4034117.3864088757,
                4: 4188551.111889635,
                5: 4320000.000000004,
            },
            "traps": {
                1: 0.003871015718248886,
                2: 0.0044572367496801485,
                3: 0.004818940668434689,
                4: 0.005719767610398479,
                5: 0.00669415633895397,
            },
            "c_norm_avg": {"Cl_-": 0.02556111708571618},
            "c_breakthru": {"Cl_-": 2.500000002016729e-07},
            "freundlich_n": 1.2,
            "mass_transfer_coeff": 0.159346300525143,
            "bv": 18000,
            "bv_50": 20000,
            "bed_capacity_param": 311.9325370632754,
            "kinetic_param": 1.5934630052514297e-06,
        }

        for k, v in results_dict.items():
            var = getattr(ix, k)
            if isinstance(v, dict):
                for i, u in v.items():
                    assert pytest.approx(u, rel=1e-3) == value(var[i])
            else:
                assert pytest.approx(v, rel=1e-3) == value(var)

    @pytest.mark.component
    def test_costing(self, IX_fr):
        m = IX_fr
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        ix.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(9701947.4187, rel=1e-3) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(1448862.0602, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(19403894.837, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.238664, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.04382530, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )


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

        ix.freundlich_n.fix(1.2)
        ix.bv_50.fix(20000)
        ix.bv.fix(18000)
        ix.resin_bulk_dens.fix(0.72)
        ix.regen_dose.fix()
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
        assert m.fs.ix.config.regenerant is RegenerantChem.NaCl
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.anion
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.freundlich
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_default_build(self, IX_inert):
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
            "underdrain_h",
            "distributor_h",
            "Pe_p_A",
            "Pe_p_exp",
            "Sh_A",
            "Sh_exp_A",
            "Sh_exp_B",
            "Sh_exp_C",
            "bed_expansion_frac_A",
            "bed_expansion_frac_B",
            "bed_expansion_frac_C",
            "p_drop_A",
            "p_drop_B",
            "p_drop_C",
            "pump_efficiency",
            "t_regen",
            "rinse_bv",
            "bw_rate",
            "t_bw",
            "service_to_regen_flow_ratio",
            "number_columns_redund",
            "c_trap_min",
        ]

        for p in ix_params:
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_diam",
            "resin_bulk_dens",
            "resin_surf_per_vol",
            "regen_dose",
            "c_norm",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "col_height_to_diam_ratio",
            "number_columns",
            "t_breakthru",
            "t_contact",
            "ebct",
            "vel_bed",
            "vel_inter",
            "service_flow_rate",
            "N_Re",
            "N_Sc",
            "N_Sh",
            "N_Pe_particle",
            "N_Pe_bed",
            "c_traps",
            "tb_traps",
            "traps",
            "c_norm_avg",
            "c_breakthru",
            "freundlich_n",
            "mass_transfer_coeff",
            "bv",
            "bv_50",
            "bed_capacity_param",
            "kinetic_param",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 90
        assert number_total_constraints(m) == 56
        assert number_unused_variables(m) == 15

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
    def test_solution(self, IX_inert):
        m = IX_inert
        ix = m.fs.ix
        results_dict = {
            "resin_diam": 0.0006749999999999999,
            "resin_bulk_dens": 0.72,
            "resin_surf_per_vol": 4444.444444444445,
            "regen_dose": 300,
            "c_norm": {"Cl_-": 0.25},
            "bed_vol_tot": 120.0,
            "bed_depth": 1.4759999999999998,
            "bed_porosity": 0.5,
            "col_height": 3.1607901999999997,
            "col_diam": 2.543563078403381,
            "col_height_to_diam_ratio": 1.242662400172933,
            "number_columns": 16,
            "t_breakthru": 4320000.0,
            "t_contact": 120.0,
            "ebct": 240.0,
            "vel_bed": 0.006149999999999999,
            "vel_inter": 0.012299999999999998,
            "service_flow_rate": 15,
            "N_Re": 4.151249999999999,
            "N_Sc": {"Cl_-": 999.9999999999998},
            "N_Sh": {"Cl_-": 24.083093218519274},
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
                1: 3344557.5805670363,
                2: 3825939.149532493,
                3: 4034117.3864088715,
                4: 4188551.1118896315,
                5: 4320000.000000001,
            },
            "traps": {
                1: 0.0038710157182488838,
                2: 0.004457236749680154,
                3: 0.004818940668434693,
                4: 0.005719767610398498,
                5: 0.0066941563389539965,
            },
            "c_norm_avg": {"Cl_-": 0.025561117085716227},
            "c_breakthru": {"Cl_-": 2.5000000016147893e-07},
            "freundlich_n": 1.2,
            "mass_transfer_coeff": 0.159346300525143,
            "bv": 18000,
            "bv_50": 20000,
            "bed_capacity_param": 311.93253706327357,
            "kinetic_param": 1.59346300525143e-06,
        }

        for k, v in results_dict.items():
            var = getattr(ix, k)
            if isinstance(v, dict):
                for i, u in v.items():
                    assert pytest.approx(u, rel=1e-3) == value(var[i])
            else:
                assert pytest.approx(v, rel=1e-3) == value(var)

    @pytest.mark.component
    def test_costing(self, IX_inert):
        m = IX_inert
        ix = m.fs.ix

        m.fs.costing = WaterTAPCosting()
        ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ix.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            ix.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        ix.initialize()

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(9701947.4187, rel=1e-3) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(695243.5958, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(19403894.837, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.18559, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.04382530, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )
