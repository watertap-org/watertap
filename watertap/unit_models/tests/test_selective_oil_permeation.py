#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
    assert_optimal_termination,
    units,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.selective_oil_permeation import SelectiveOilPermeation
import watertap.property_models.selective_oil_permeation_prop_pack as props
from watertap.core.util.initialization import assert_no_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.SopParameterBlock()
    m.fs.unit = SelectiveOilPermeation(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.is_isothermal
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.none
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)


class TestSelectiveOilPermeation:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.SopParameterBlock()
        m.fs.unit = SelectiveOilPermeation(property_package=m.fs.properties)

        # fully specify system
        m.fs.unit.feed_side.properties_in[0].temperature.fix(298.15)  # temp in K
        m.fs.unit.feed_side.properties_in[0].pressure.fix(2.52 * units.bar)
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            6.27e-2
        )  # H2O flow in kg/s
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "oil"].fix(
            4.94e-4
        )  # oil flow in kg/s
        m.fs.unit.feed_side.properties_out[0].pressure.fix(2.24 * units.bar)
        m.fs.unit.area.fix(1.4)  # membrane area in m^2
        m.fs.unit.pore_diameter.fix(4.7e-8)
        m.fs.unit.porosity.fix(0.4)
        m.fs.unit.membrane_thickness.fix(4e-5)
        m.fs.unit.permeability_constant.fix(150)
        m.fs.unit.module_diameter.fix(0.064)
        m.fs.unit.properties_permeate[0].pressure.fix(1 * units.bar)

        return m

    @pytest.mark.unit
    def test_build(self, unit_frame):
        m = unit_frame

        # test ports
        port_lst = ["inlet", "retentate", "permeate"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert (
                len(port.vars) == 3
            )  # number of state variables for SOP property package
            assert isinstance(port, Port)

        assert number_variables(m) == 35
        assert number_total_constraints(m) == 22
        assert number_unused_variables(m) == 0
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        m = unit_frame
        assert_no_degrees_of_freedom(m)

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e5, index=("Liq", "oil")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, unit_frame):
        m = unit_frame

        # set initial values (to ensure these variables are not badly scaled)
        m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp[
            "Liq", "oil"
        ].value = 1e-5
        m.fs.unit.properties_permeate[0].flow_mass_phase_comp["Liq", "oil"].value = 1e-6
        m.fs.unit.feed_side.properties_in[0].vol_frac_phase_comp[
            "Liq", "oil"
        ].value = 1e-3

        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        m.fs.unit.initialize(outlvl=idaeslog.DEBUG)
        initialization_tester(unit_frame)

    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, unit_frame):
        m = unit_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.component_list

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.properties_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-6
        )

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert pytest.approx(3.32237e-4, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].flow_mass_phase_comp["Liq", "oil"]
        )
        assert pytest.approx(1.38e5, rel=1e-3) == value(
            m.fs.unit.pressure_transmemb_avg[0]
        )
        assert pytest.approx(-2.8e4, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(3.0425e-07, rel=1e-3) == value(m.fs.unit.flux_vol_oil[0])
        assert pytest.approx(3.3246e-06, rel=1e-3) == value(
            m.fs.unit.flux_vol_oil_pure[0]
        )
        assert pytest.approx(0.091512, rel=1e-3) == value(
            m.fs.unit.effective_area_ratio[0]
        )
        assert pytest.approx(0.10799, rel=1e-3) == value(
            m.fs.unit.effective_area_ratio_num[0]
        )
        assert pytest.approx(1.1801, rel=1e-3) == value(
            m.fs.unit.effective_area_ratio_den[0]
        )
        assert pytest.approx(0.67254, rel=1e-3) == value(m.fs.unit.recovery_frac_oil[0])
        assert pytest.approx(0.019687, rel=1e-3) == value(
            m.fs.unit.liquid_velocity_in[0]
        )
        assert pytest.approx(3.32237e-4, rel=1e-3) == value(
            m.fs.unit.mass_transfer_oil[0]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()
