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
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.uv_aop import Ultraviolet0D
import watertap.property_models.NDMA_prop_pack as props

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
from pyomo.util.check_units import assert_units_consistent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestUltraviolet:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NDMAParameterBlock()

        m.fs.unit = Ultraviolet0D(default={"property_package": m.fs.properties})

        # fully specify system
        feed_flow_mass = 1 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        exporure_time = 500 * pyunits.s
        inactivation_rate = 2.8 * pyunits.cm**2 / pyunits.J
        photolysis_rate_constant = 0.16 * pyunits.min**-1
        EEO = 0.25 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.3
        UVT = 0.9

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        # m.fs.unit.uv_dose.fix(uv_dose)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.exposure_time.fix(exporure_time)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
        m.fs.unit.photolysis_rate_constant["Liq", "NDMA"].fix(photolysis_rate_constant)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
        m.fs.unit.UVT.fix(UVT)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
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

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 35
        assert number_total_constraints(m) == 22
        assert number_unused_variables(m) == 0  # vars from property package parameters

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        [print(i[0]) for i in badly_scaled_var_lst]
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(0.999999999, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.001003, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_vol
        )
        assert pytest.approx(74e-9, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(0.999999999, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )
        assert pytest.approx(1.8248e-8, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(5000, rel=1e-3) == value(m.fs.unit.uv_dose)
        assert pytest.approx(0.04576, rel=1e-3) == value(m.fs.unit.UVA)
        assert pytest.approx(1.3333e-4, rel=1e-3) == value(
            m.fs.unit.reaction_rate_constant["Liq", "NDMA"]
        )
        assert pytest.approx(1829.525, rel=1e-3) == value(
            m.fs.unit.electricity_demand_phase_comp[0, "Liq", "NDMA"]
        )
