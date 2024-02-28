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
    Var,
    Expression,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.pressure_exchanger import (
    PressureExchanger,
    PressureExchangeType,
)
import watertap.property_models.seawater_prop_pack as props
import watertap.property_models.multicomp_aq_sol_prop_pack as props_multi
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# # -----------------------------------------------------------------------------


@pytest.mark.unit
def test_config_no_mass_transfer():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = PressureExchanger(property_package=m.fs.properties)

    # check unit config arguments
    assert len(m.fs.unit.config) == 11

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.none
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.property_package is m.fs.properties

    # verify no mass transfer
    assert not m.fs.unit.config.has_mixing


@pytest.mark.unit
def test_config_mass_transfer():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = PressureExchanger(property_package=m.fs.properties, has_mixing=True)
    # check mass_transfer is added
    assert m.fs.unit.config.has_mixing


@pytest.mark.unit
def test_build(
    has_leakage=False, has_mixing=False, extra_variables=0, extra_constraint=0
):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = PressureExchanger(
        property_package=m.fs.properties, has_leakage=has_leakage, has_mixing=has_mixing
    )

    # test ports and state variables
    port_lst = [
        "low_pressure_inlet",
        "low_pressure_outlet",
        "high_pressure_inlet",
        "high_pressure_outlet",
    ]
    port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]

    for port_str in port_lst:
        assert hasattr(m.fs.unit, port_str)
        port = getattr(m.fs.unit, port_str)
        assert len(port.vars) == 3
        assert isinstance(port, Port)
        for var_str in port_vars_lst:
            assert hasattr(port, var_str)
            var = getattr(port, var_str)
            assert isinstance(var, Var)

    # test unit variables
    assert hasattr(m.fs.unit, "efficiency_pressure_exchanger")
    assert isinstance(m.fs.unit.efficiency_pressure_exchanger, Var)
    if not has_leakage:
        assert not hasattr(m.fs.unit, "leakage_vol")
    else:
        assert isinstance(m.fs.unit.leakage_vol, Var)

    if not has_mixing:
        assert not hasattr(m.fs.unit, "mixing_vol")
    else:
        assert isinstance(m.fs.unit.mixing_vol, Var)

    # test unit constraints

    unit_cons_lst = [
        "eq_pressure_transfer",
        "eq_equal_flow_vol",
        "eq_equal_low_pressure",
    ]

    if has_leakage or has_mixing:
        # unit_cons_lst.append("eq_mass_transfer_from_high_pressure_side")
        unit_cons_lst.append("eq_leakage")
        unit_cons_lst.append("eq_mixing")
        unit_cons_lst.append("eq_mass_transfer_term")

    for c in unit_cons_lst:
        assert hasattr(m.fs.unit, c)
        con = getattr(m.fs.unit, c)
        assert isinstance(con, Constraint)

    # test control volumes, only terms directly used by pressure exchanger

    cv_list = ["low_pressure_side", "high_pressure_side"]
    cv_var_lst = ["deltaP"]
    if has_mixing:
        cv_var_lst.append("mass_transfer_term")
    cv_exp_lst = ["work"]

    for cv_str in cv_list:
        assert hasattr(m.fs.unit, cv_str)
        cv = getattr(m.fs.unit, cv_str)
        for cv_var_str in cv_var_lst:
            cv_var = getattr(cv, cv_var_str)
            assert isinstance(cv_var, Var)
        for cv_exp_str in cv_exp_lst:
            cv_exp = getattr(cv, cv_exp_str)
            assert isinstance(cv_exp, Expression)

    # test state blocks, only terms directly used by pressure exchanger
    cv_stateblock_lst = ["properties_in", "properties_out"]
    stateblock_var_lst = ["pressure", "temperature"]
    stateblock_exp_lst = ["flow_vol"]
    for cv_str in cv_list:
        cv = getattr(m.fs.unit, cv_str)
        for cv_sb_str in cv_stateblock_lst:
            assert hasattr(cv, cv_sb_str)
            cv_sb = getattr(cv, cv_sb_str)
            for sb_var_str in stateblock_var_lst:
                assert hasattr(cv_sb[0], sb_var_str)
                sb_var = getattr(cv_sb[0], sb_var_str)
                assert isinstance(sb_var, Var)
            for sb_exp_str in stateblock_exp_lst:
                assert hasattr(cv_sb[0], sb_exp_str)
                sb_exp = getattr(cv_sb[0], sb_exp_str)
                assert isinstance(sb_exp, Expression)

    # test statistics
    assert number_variables(m.fs.unit) == 39 + extra_variables
    assert number_total_constraints(m.fs.unit) == 31 + extra_constraint
    assert number_unused_variables(m.fs.unit) == 0


@pytest.mark.unit
def test_build_without_mass_transfer():
    test_build(has_leakage=False, has_mixing=False)


class TestPressureExchanger_without_mass_transfer:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = PressureExchanger(property_package=m.fs.properties)

        # Specify inlet conditions
        temperature = 25 + 273.15
        flow_vol = 1e-3
        lowP_mass_frac_TDS = 0.035
        lowP_pressure = 101325
        highP_mass_frac_TDS = 0.07
        highP_pressure = 50e5

        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(flow_vol)
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(lowP_mass_frac_TDS)

        m.fs.unit.low_pressure_side.properties_in[0].pressure.fix(lowP_pressure)
        m.fs.unit.low_pressure_side.properties_in[0].temperature.fix(temperature)

        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(
            flow_vol
        )
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(highP_mass_frac_TDS)
        m.fs.unit.high_pressure_side.properties_in[0].pressure.fix(highP_pressure)
        m.fs.unit.high_pressure_side.properties_in[0].temperature.fix(temperature)

        # solve inlet conditions and only fix state variables (i.e. unfix flow_vol and mass_frac_phase)

        results = solver.solve(m.fs.unit.high_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.high_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        results = solver.solve(m.fs.unit.low_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        # Specify unit
        efficiency = 0.95
        m.fs.unit.efficiency_pressure_exchanger.fix(efficiency)
        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        initialization_tester(unit_frame)

    @pytest.mark.component
    def test_var_scaling(self, unit_frame):
        m = unit_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert pytest.approx(0.9877, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(3.582e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.755e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.pressure[0]
        )
        assert pytest.approx(0.9767, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(7.352e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.654e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.deltaP[0]
        )
        assert pytest.approx(4.654e3, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.work[0]
        )
        assert pytest.approx(-4.899e6, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.deltaP[0]
        )
        assert pytest.approx(-4.899e3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.work[0]
        )

        # testing solvent transfer
        assert pytest.approx(0.9767, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.9877, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        # testing solute transfer
        assert pytest.approx(7.352e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(3.582e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()


class TestPressureExchanger_with_high_pressure_difference:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = PressureExchanger(
            property_package=m.fs.properties,
            pressure_exchange_calculation=PressureExchangeType.high_pressure_difference,
        )

        # Specify inlet conditions
        temperature = 25 + 273.15
        flow_vol = 1e-3
        lowP_mass_frac_TDS = 0.035
        lowP_pressure = 101325
        highP_mass_frac_TDS = 0.07
        highP_pressure = 50e5

        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(flow_vol)
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(lowP_mass_frac_TDS)

        m.fs.unit.low_pressure_side.properties_in[0].pressure.fix(lowP_pressure)
        m.fs.unit.low_pressure_side.properties_in[0].temperature.fix(temperature)

        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(
            flow_vol
        )
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(highP_mass_frac_TDS)
        m.fs.unit.high_pressure_side.properties_in[0].pressure.fix(highP_pressure)
        m.fs.unit.high_pressure_side.properties_in[0].temperature.fix(temperature)

        # solve inlet conditions and only fix state variables (i.e. unfix flow_vol and mass_frac_phase)

        results = solver.solve(m.fs.unit.high_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.high_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        results = solver.solve(m.fs.unit.low_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        # Specify unit
        high_pressure_difference = 0.8 * pyunits.bar
        low_pressure_difference = 0 * pyunits.bar
        m.fs.unit.high_pressure_difference.fix(high_pressure_difference)
        m.fs.unit.low_pressure_difference.fix(low_pressure_difference)
        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        initialization_tester(unit_frame)

    @pytest.mark.component
    def test_var_scaling(self, unit_frame):
        m = unit_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert pytest.approx(0.9877, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(3.582e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.92e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.pressure[0]
        )
        assert pytest.approx(0.9767, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(7.352e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.819e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.deltaP[0]
        )
        assert pytest.approx(4.819e3, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.work[0]
        )
        assert pytest.approx(-4.899e6, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.deltaP[0]
        )
        assert pytest.approx(-4.899e3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.work[0]
        )

        # testing solvent transfer
        assert pytest.approx(0.9767, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.9877, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        # testing solute transfer
        assert pytest.approx(7.352e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(3.582e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()


class TestPressureExchanger_with_mass_transfer:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = PressureExchanger(
            property_package=m.fs.properties,
            has_leakage=True,
            has_mixing=True,
            pressure_exchange_calculation=PressureExchangeType.high_pressure_difference,
        )

        # Specify inlet conditions
        temperature = 25 + 273.15
        flow_vol = 1e-3
        lowP_mass_frac_TDS = 0.035
        lowP_pressure = 101325
        highP_mass_frac_TDS = 0.07
        highP_pressure = 50e5
        leakage_vol = 0.01
        mixing_vol = 0.035

        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(flow_vol)
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(lowP_mass_frac_TDS)

        m.fs.unit.low_pressure_side.properties_in[0].pressure.fix(lowP_pressure)
        m.fs.unit.low_pressure_side.properties_in[0].temperature.fix(temperature)

        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].fix(
            flow_vol
        )
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].fix(highP_mass_frac_TDS)
        m.fs.unit.high_pressure_side.properties_in[0].pressure.fix(highP_pressure)
        m.fs.unit.high_pressure_side.properties_in[0].temperature.fix(temperature)

        m.fs.unit.leakage_vol[0].fix(leakage_vol)
        m.fs.unit.mixing_vol[0].fix(mixing_vol)

        # solve inlet conditions and only fix state variables (i.e. unfix flow_vol and mass_frac_phase)

        results = solver.solve(m.fs.unit.high_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.high_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        results = solver.solve(m.fs.unit.low_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase["Liq"].unfix()
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp[
            "Liq", "TDS"
        ].unfix()

        # Specify unit
        high_pressure_difference = 0.8 * pyunits.bar
        low_pressure_difference = 0 * pyunits.bar
        m.fs.unit.high_pressure_difference.fix(high_pressure_difference)
        m.fs.unit.low_pressure_difference.fix(low_pressure_difference)
        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        initialization_tester(unit_frame)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert pytest.approx(0.9877, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(3.582e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.92e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.pressure[0]
        )
        assert pytest.approx(0.9868, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(7.352e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(4.819e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.deltaP[0]
        )
        assert pytest.approx(4.819e3, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.work[0]
        )
        assert pytest.approx(-4.899e6, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.deltaP[0]
        )
        assert pytest.approx(-4.948e3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.work[0]
        )

        # testing solvent transfer
        assert pytest.approx(3.559e-4, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.9871, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.9874, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        # testing solute transfer
        assert pytest.approx(-1.293e-3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "TDS"]
        )
        assert pytest.approx(7.222e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )
        assert pytest.approx(3.712e-2, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()


class TestPressureExchanger_with_ion_prop_pack:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props_multi.MCASParameterBlock(
            solute_list=["Na", "Ca", "Mg", "SO4", "Cl"],
            mw_data={
                "H2O": 0.018,
                "Na": 0.023,
                "Ca": 0.04,
                "Mg": 0.024,
                "SO4": 0.096,
                "Cl": 0.035,
            },
            ignore_neutral_charge=True,
        )
        m.fs.unit = PressureExchanger(
            property_package=m.fs.properties,
            has_leakage=True,
            has_mixing=True,
        )

        # Specify inlet conditions
        temperature = 25 + 273.15

        lowP_pressure = 101325
        highP_pressure = 50e5

        feed_flow_mass = 1 * pyunits.kg / pyunits.s
        lowP_mass_frac = {
            "Na": 11122e-6,
            "Ca": 382e-6,
            "Mg": 1394e-6,
            "SO4": 2136e-6,
            "Cl": 20300e-6,
        }

        highP_mass_frac = {
            "Na": 2 * 11122e-6,
            "Ca": 2 * 382e-6,
            "Mg": 2 * 1394e-6,
            "SO4": 2 * 2136e-6,
            "Cl": 2 * 20300e-6,
        }

        leakage_vol = 0.01
        mixing_vol = 0.035

        # Specify unit
        efficiency = 0.95
        m.fs.unit.efficiency_pressure_exchanger.fix(efficiency)

        m.fs.unit.leakage_vol[0].fix(leakage_vol)
        m.fs.unit.mixing_vol[0].fix(mixing_vol)

        LPS_in = m.fs.unit.low_pressure_side.properties_in[0]
        for ion, x in lowP_mass_frac.items():
            mol_comp_flow = (
                x * pyunits.kg / pyunits.kg * feed_flow_mass / LPS_in.mw_comp[ion]
            )

            LPS_in.flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)

        LPS_in.flow_mol_phase_comp["Liq", "H2O"].fix(
            (feed_flow_mass * (1 - sum(x for x in lowP_mass_frac.values())))
            / LPS_in.mw_comp["H2O"]
        )

        m.fs.unit.low_pressure_side.properties_in[0].pressure.fix(lowP_pressure)
        m.fs.unit.low_pressure_side.properties_in[0].temperature.fix(temperature)

        ### NOTE THE PRESSURE EXHCANGER EXPECTS EQUAL VOLUMETRIC FLOW RATES
        # so one of mass fraction (here we unfix H2O mass fraction)
        # s on high pressure or low pressure sides
        # should be left unfixed, as other wise you are fixing all the mass fractions,
        # and forcing equal flow, which leaves no degrees of freedom
        # to satisfy the equal flow constriants, which will lead to
        # a -1 DOF error and should not occur in practice where typicly at least one of
        # the inlet mass flow rates is colaculated.

        HPS_in = m.fs.unit.high_pressure_side.properties_in[0]
        for ion, x in highP_mass_frac.items():
            mol_comp_flow = (
                x * pyunits.kg / pyunits.kg * feed_flow_mass / HPS_in.mw_comp[ion]
            )

            HPS_in.flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)

        HPS_in.flow_mol_phase_comp["Liq", "H2O"] = (
            feed_flow_mass
            * (1 - sum(x for x in highP_mass_frac.values()))
            / HPS_in.mw_comp["H2O"]
        )

        m.fs.unit.high_pressure_side.properties_in[0].pressure.fix(highP_pressure)
        m.fs.unit.high_pressure_side.properties_in[0].temperature.fix(temperature)

        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        initialization_tester(unit_frame)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame

        assert pytest.approx(4.755e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_outlet.pressure[0]
        )

        assert pytest.approx(4.654e6, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.deltaP[0]
        )
        assert pytest.approx(4.654e3, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.work[0]
        )
        assert pytest.approx(-4.899e6, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.deltaP[0]
        )
        assert pytest.approx(-4.948e3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.work[0]
        )

        # testing solvent transfer
        assert pytest.approx(0.06733, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.5, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )
        assert pytest.approx(0.9634, rel=1e-3) == value(
            m.fs.unit.low_pressure_side.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )

        # testing solute transfer
        assert pytest.approx(-0.01659, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "Na"]
        )
        assert pytest.approx(-3.276e-4, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "Ca"]
        )
        assert pytest.approx(-1.992e-3, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "Mg"]
        )
        assert pytest.approx(-7.632e-4, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "SO4"]
        )
        assert pytest.approx(-1.989e-2, rel=1e-3) == value(
            m.fs.unit.high_pressure_side.mass_transfer_term[0, "Liq", "Cl"]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()
