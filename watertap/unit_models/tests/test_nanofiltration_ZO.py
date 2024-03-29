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
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.nanofiltration_ZO import NanofiltrationZO

import watertap.property_models.multicomp_aq_sol_prop_pack as props

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
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.MCASParameterBlock(solute_list=["Na_+", "Ca_2+", "Mg_2+", "SO4_2-", "Cl_-"])
    m.fs.unit = NanofiltrationZO(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.none
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.MCASParameterBlock(solute_list=["Na_+", "Ca_2+", "Mg_2+", "SO4_2-", "Cl_-"])
    m.fs.unit = NanofiltrationZO(
        property_package=m.fs.properties, has_pressure_change=True
    )

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)


class TestNanofiltration:
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties =  props.MCASParameterBlock(solute_list=["Na_+", "Ca_2+", "Mg_2+", "SO4_2-", "Cl_-"],
                                                    material_flow_basis=props.MaterialFlowBasis.mass)

        m.fs.unit = NanofiltrationZO(property_package=m.fs.properties)

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac = {
            "Na_+": 11122e-6,
            "Ca_2+": 382e-6,
            "Mg_2+": 1394e-6,
            "SO4_2-": 2136e-6,
            "Cl_-": 20316.88e-6,
        }
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            feed_flow_mass * (1 - sum(x for x in feed_mass_frac.values()))
        )
        for j in feed_mass_frac:
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j].fix(
                feed_flow_mass * feed_mass_frac[j]
            )
        m.fs.unit.feed_side.properties_in[0].pressure.fix(101325)
        m.fs.unit.feed_side.properties_in[0].temperature.fix(298.15)

        m.fs.unit.flux_vol_solvent.fix(1.67e-6)
        m.fs.unit.area.fix(500)
        m.fs.unit.properties_permeate[0].pressure.fix(101325)

        m.fs.unit.rejection_phase_comp[0, "Liq", "Na_+"].fix(0.01)
        m.fs.unit.rejection_phase_comp[0, "Liq", "Ca_2+"].fix(0.79)
        m.fs.unit.rejection_phase_comp[0, "Liq", "Mg_2+"].fix(0.94)
        m.fs.unit.rejection_phase_comp[0, "Liq", "SO4_2-"].fix(0.87)
        m.fs.unit.rejection_phase_comp[0, "Liq", "Cl_-"].fix(0.15)  
       
        m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(defined_state=True, adjust_by_ion='Cl_-')

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
            )  # number of state variables for NaCl property package
            assert isinstance(port, Port)

        # test statistics
        assert number_variables(m) == 117
        assert number_total_constraints(m) == 86
        assert number_unused_variables(m) == 16
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        m = unit_frame
        assert_no_degrees_of_freedom(m)

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e4, index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "Mg_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "SO4_2-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, unit_frame):
        initialization_tester(unit_frame)

    @pytest.mark.component
    def test_var_scaling(self, unit_frame):
        m = unit_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        solver.options = {"nlp_scaling_method": "user-scaling"}
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
            <= 1e-5
        )

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert pytest.approx(0.8350, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(2.002, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].conc_mol_phase_comp["Liq", "Ca_2+"]
        )
        assert pytest.approx(487.116, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].conc_mol_phase_comp["Liq", "Cl_-"]
        )
        assert pytest.approx(3.441, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].conc_mol_phase_comp["Liq", "Mg_2+"]
        )
        assert pytest.approx(478.9, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].conc_mol_phase_comp["Liq", "Na_+"]
        )
        assert pytest.approx(2.891, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].conc_mol_phase_comp["Liq", "SO4_2-"]
        )
        assert pytest.approx(0.86, rel=1e-3) == value(
            m.fs.unit.recovery_vol_phase[0, "Liq"]
        )

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()
