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
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.nanofiltration_0D import NanoFiltration0D
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.solvers import get_solver
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

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestNanoFiltration:
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = NanoFiltration0D(
            property_package=m.fs.properties, has_pressure_change=True
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 6e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 1e5
        membrane_area = 50 * feed_flow_mass
        A = 3.77e-11
        B = 4.724e-5
        sigma = 0.28
        pressure_atmospheric = 101325

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.sigma.fix(sigma)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        return m

    @pytest.mark.unit
    def test_config(self, NF_frame):
        m = NF_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.has_pressure_change
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ["inlet", "retentate", "permeate"]
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
            "A_comp",
            "B_comp",
            "sigma",
            "dens_solvent",
            "flux_mass_phase_comp_in",
            "flux_mass_phase_comp_out",
            "avg_conc_mass_phase_comp_in",
            "avg_conc_mass_phase_comp_out",
            "area",
            "deltaP",
            "mass_transfer_phase_comp",
            "flux_mass_phase_comp_avg",
            "eq_mass_transfer_term",
            "eq_permeate_production",
            "eq_flux_in",
            "eq_flux_out",
            "eq_avg_conc_in",
            "eq_avg_conc_out",
            "eq_connect_mass_transfer",
            "eq_connect_enthalpy_transfer",
            "eq_permeate_isothermal",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "feed_side"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "pressure_osm_phase",
            "osm_coeff",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "enth_mass_phase",
            "eq_pressure_osm_phase",
            "eq_osm_coeff",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
            "eq_enth_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)
        # permeate
        assert hasattr(m.fs.unit, "properties_permeate")
        blk = getattr(m.fs.unit, "properties_permeate")
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # test statistics
        assert number_variables(m) == 73
        assert number_total_constraints(m) == 45
        assert number_unused_variables(m) == 7  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, NF_frame):
        initialization_tester(NF_frame)

    @pytest.mark.component
    def test_var_scaling(self, NF_frame):
        m = NF_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

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

        assert (
            abs(
                value(
                    flow_mass_inlet
                    * b.feed_side.properties_in[0].enth_mass_phase["Liq"]
                    - flow_mass_retentate
                    * b.feed_side.properties_out[0].enth_mass_phase["Liq"]
                    - flow_mass_permeate
                    * b.properties_permeate[0].enth_mass_phase["Liq"]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame
        assert pytest.approx(1.079e-2, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(3.435e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.5396, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(1.717e-2, rel=1e-3) == value(
            m.fs.unit.properties_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
