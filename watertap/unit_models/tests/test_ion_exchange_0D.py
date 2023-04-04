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
    Objective,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    MCASStateBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
    IsothermType,
    DiffusionControlType,
)
from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof

from idaes.core import (
    EnergyBalanceType,
    MomentumBalanceType,
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
    set_scaling_factor,
)
from pyomo.util.check_units import assert_units_consistent
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestIonExchangeDefault:
    @pytest.fixture(scope="class")
    def IX_default(self):
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
    def test_config(self, IX_default):
        m = IX_default

        assert len(m.fs.ix.config) == 12

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert not m.fs.ix.config.hazardous_waste
        assert isinstance(m.fs.ix.regen_chem, RegenerantChem)
        assert m.fs.ix.regen_chem is RegenerantChem.NaCl
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.cation
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.langmuir
        assert isinstance(m.fs.ix.config.diffusion_control, DiffusionControlType)
        assert m.fs.ix.config.diffusion_control is DiffusionControlType.liquid
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_default_build(self, IX_default):
        m = IX_default
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
            "p_drop_psi_to_m",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
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
            print(p)
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "langmuir",
            "num_transfer_units",
            "dimensionless_time",
            "resin_surf_per_vol",
            "col_height_to_diam_ratio",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "number_columns",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "t_breakthru",
            "t_contact",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_bed",
            "vel_inter",
            "service_flow_rate",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "regen_dose",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 72
        assert number_total_constraints(m) == 45
        assert number_unused_variables(m) == 12

    @pytest.mark.unit
    def test_dof(self, IX_default):
        m = IX_default
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_default):
        m = IX_default
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
    def test_initialize(self, IX_default):
        m = IX_default
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, IX_default):
        m = IX_default
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, IX_default):
        m = IX_default
        ix = m.fs.ix
        target_ion = ix.config.target_ion
        assert (
            abs(
                value(
                    ix.mass_in[target_ion]
                    - ix.mass_out[target_ion]
                    - ix.mass_removed[target_ion]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_solution(self, IX_default):
        m = IX_default
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
            "mass_in": 130901.61040796709,
            "mass_removed": 130601.61040796709,
            "mass_out": 300.00000000000864,
            "vel_bed": 0.007083333333333333,
            "vel_inter": 0.014166666666666666,
            "service_flow_rate": 15,
            "Re": 4.958333333333333,
            "Sc": 1086.9565217391305,
            "Sh": 26.29635815858793,
            "Pe_p": 0.10782790064157834,
            "Pe_bed": 261.86775870097597,
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
    def test_costing(self, IX_default):
        m = IX_default
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

        assert pytest.approx(3296998.040376, rel=1e-5) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(2183260.52954, rel=1e-5) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(6593996.0807, rel=1e-5) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.2001943, rel=1e-5) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.0572452, rel=1e-5) == value(
            m.fs.costing.specific_energy_consumption
        )


class TestIonExchangeSolidDiffusion:
    @pytest.fixture(scope="class")
    def IX_solid_diff(self):
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
            "diffusion_control": "solid",
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

        ix.diff_resin_comp.fix(10e-11)
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
    def test_config(self, IX_solid_diff):
        m = IX_solid_diff

        assert len(m.fs.ix.config) == 12

        assert not m.fs.ix.config.dynamic
        assert not m.fs.ix.config.has_holdup
        assert m.fs.ix.config.property_package is m.fs.properties
        assert not m.fs.ix.config.hazardous_waste
        assert isinstance(m.fs.ix.regen_chem, RegenerantChem)
        assert m.fs.ix.regen_chem is RegenerantChem.NaCl
        assert isinstance(m.fs.ix.ion_exchange_type, IonExchangeType)
        assert m.fs.ix.ion_exchange_type is IonExchangeType.cation
        assert isinstance(m.fs.ix.config.isotherm, IsothermType)
        assert m.fs.ix.config.isotherm is IsothermType.langmuir
        assert isinstance(m.fs.ix.config.diffusion_control, DiffusionControlType)
        assert m.fs.ix.config.diffusion_control is DiffusionControlType.solid
        assert isinstance(m.fs.ix.config.energy_balance_type, EnergyBalanceType)
        assert m.fs.ix.config.energy_balance_type is EnergyBalanceType.none
        assert isinstance(m.fs.ix.config.momentum_balance_type, MomentumBalanceType)
        assert m.fs.ix.config.momentum_balance_type is MomentumBalanceType.pressureTotal

    @pytest.mark.unit
    def test_default_build(self, IX_solid_diff):
        m = IX_solid_diff
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
            "p_drop_psi_to_m",
            "holdup_A",
            "holdup_B",
            "holdup_exp",
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
            "phi_solid_coeff_A",
            "phi_solid_coeff_B",
            "phi_solid_coeff_C",
        ]

        for p in ix_params:
            print(p)
            assert hasattr(ix, p)
            param = getattr(ix, p)
            assert isinstance(param, Param)

        ix_vars = [
            "resin_max_capacity",
            "resin_eq_capacity",
            "resin_unused_capacity",
            "resin_diam",
            "resin_bulk_dens",
            "langmuir",
            "num_transfer_units",
            "dimensionless_time",
            "diff_resin_comp",
            "resin_surf_per_vol",
            "col_height_to_diam_ratio",
            "bed_vol_tot",
            "bed_depth",
            "bed_porosity",
            "col_height",
            "col_diam",
            "number_columns",
            "partition_ratio",
            "fluid_mass_transfer_coeff",
            "t_breakthru",
            "t_contact",
            "mass_in",
            "mass_removed",
            "mass_out",
            "vel_bed",
            "vel_inter",
            "service_flow_rate",
            "Re",
            "Sc",
            "Sh",
            "Pe_p",
            "Pe_bed",
            "c_norm",
            "regen_dose",
        ]

        for v in ix_vars:
            assert hasattr(ix, v)
            var = getattr(ix, v)
            assert isinstance(var, Var)

        # test statistics
        assert number_variables(m) == 73
        assert number_total_constraints(m) == 45
        assert number_unused_variables(m) == 12

    @pytest.mark.unit
    def test_dof(self, IX_solid_diff):
        m = IX_solid_diff
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, IX_solid_diff):
        m = IX_solid_diff
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
    def test_initialize(self, IX_solid_diff):
        m = IX_solid_diff
        initialization_tester(m, unit=m.fs.ix, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, IX_solid_diff):
        m = IX_solid_diff
        results = solver.solve(m, tee=True)
        assert_units_consistent(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, IX_solid_diff):
        m = IX_solid_diff
        ix = m.fs.ix
        target_ion = ix.config.target_ion
        assert (
            abs(
                value(
                    ix.mass_in[target_ion]
                    - ix.mass_out[target_ion]
                    - ix.mass_removed[target_ion]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_solution(self, IX_solid_diff):
        m = IX_solid_diff
        ix = m.fs.ix
        target_ion = ix.config.target_ion

        results_dict = {
            "resin_max_capacity": 3,
            "resin_eq_capacity": 1.6030726333501075,
            "resin_unused_capacity": 1.3969273666498925,
            "resin_diam": 0.0007,
            "resin_bulk_dens": 0.7,
            "langmuir": 0.9,
            "num_transfer_units": 659.549883435473,
            "dimensionless_time": 1,
            "diff_resin_comp": 1e-10,
            "resin_surf_per_vol": 4285.714285714286,
            "col_height_to_diam_ratio": 1.0408526790314099,
            "bed_vol_tot": 120.00000000000001,
            "bed_depth": 1.7,
            "bed_porosity": 0.5,
            "col_height": 3.488715,
            "col_diam": 3.3517855795370646,
            "number_columns": 8,
            "partition_ratio": 224.4301686690151,
            "fluid_mass_transfer_coeff": 3.456092786557271e-05,
            "t_breakthru": 53983.240480563625,
            "t_contact": 120.0,
            "mass_in": 134958.10120140907,
            "mass_removed": 134658.101201409,
            "mass_out": 300.00000000005065,
            "vel_bed": 0.007083333333333333,
            "vel_inter": 0.014166666666666666,
            "service_flow_rate": 15,
            "Re": 4.958333333333333,
            "Sc": 1086.9565217391305,
            "Sh": 26.29635815858793,
            "Pe_p": 0.10782790064157834,
            "Pe_bed": 261.86775870097597,
            "c_norm": 0.5080709442210705,
            "regen_dose": 300,
        }

        for v, val in results_dict.items():
            var = getattr(ix, v)
            if var.is_indexed():
                assert pytest.approx(val, rel=1e-3) == value(var[target_ion])
            else:
                assert pytest.approx(val, rel=1e-3) == value(var)

    @pytest.mark.component
    def test_costing(self, IX_solid_diff):
        m = IX_solid_diff
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

        assert pytest.approx(3296998.040376, rel=1e-5) == value(
            m.fs.costing.aggregate_capital_cost
        )
        assert pytest.approx(2129395.91276, rel=1e-5) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(6593996.0807, rel=1e-5) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(0.1964009, rel=1e-5) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.0572452, rel=1e-5) == value(
            m.fs.costing.specific_energy_consumption
        )
