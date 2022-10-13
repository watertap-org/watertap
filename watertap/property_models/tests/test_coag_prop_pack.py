###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import pytest
from watertap.property_models.coagulation_prop_pack import (
    CoagulationParameterBlock,
    CoagulationStateBlock,
)
from watertap.property_models.tests.property_test_harness import PropertyAttributeError
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    units as pyunits,
    Suffix,
    Constraint,
    SolverFactory,
    SolverStatus,
    TerminationCondition,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.solvers import get_solver

__author__ = "Austin Ladshaw"

solver = get_solver()

# -----------------------------------------------------------------------------
# Start test class
class TestCoagulationPropPack:
    @pytest.fixture(scope="class")
    def coag_obj(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.properties = CoagulationParameterBlock()

        return model

    @pytest.mark.unit
    def test_build_properties(self, coag_obj):
        model = coag_obj

        # Check to make sure we have the correct components and phases
        assert isinstance(model.fs.properties.component_list, Set)
        for j in model.fs.properties.component_list:
            assert j in ["H2O", "TDS", "TSS", "Sludge"]
        assert isinstance(model.fs.properties.phase_list, Set)
        for p in model.fs.properties.phase_list:
            assert p in ["Liq"]

        # Check for existance of the state block object
        assert model.fs.properties._state_block_class is CoagulationStateBlock

        # Assert that we have the expected set of parameters
        assert isinstance(model.fs.properties.cp, Param)
        assert isinstance(model.fs.properties.ref_dens_liq, Param)
        assert isinstance(model.fs.properties.dens_slope, Param)
        assert isinstance(model.fs.properties.dens_param_A, Param)
        assert isinstance(model.fs.properties.dens_param_B, Param)
        assert isinstance(model.fs.properties.dens_param_C, Param)
        assert isinstance(model.fs.properties.ref_pressure_correction, Param)
        assert isinstance(model.fs.properties.ref_pressure_slope, Param)
        assert isinstance(model.fs.properties.mu_A, Param)
        assert isinstance(model.fs.properties.mu_B, Param)
        assert isinstance(model.fs.properties.mu_C, Param)

    @pytest.mark.unit
    def test_metadata(self, coag_obj):
        model = coag_obj

        # Create the state block and pull in the default metadata
        model.fs.stream = model.fs.properties.build_state_block([0])
        metadata = model.fs.properties.get_metadata().properties

        # check that properties are not built if not demanded
        for v_name in metadata:
            if metadata[v_name]["method"] is not None:
                if model.fs.stream[0].is_property_constructed(v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was found "
                        "on the stateblock without being demanded".format(v_name=v_name)
                    )

        # check that properties are built if demanded
        for v_name in metadata:
            if metadata[v_name]["method"] is not None:
                if not hasattr(model.fs.stream[0], v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was not built "
                        "when demanded".format(v_name=v_name)
                    )

        # check the other stateblock functions
        res = model.fs.stream[0].define_state_vars()
        assert res["flow_mass_phase_comp"] == model.fs.stream[0].flow_mass_phase_comp
        assert res["pressure"] == model.fs.stream[0].pressure
        assert res["temperature"] == model.fs.stream[0].temperature
        assert (
            model.fs.stream[0].get_material_flow_terms("Liq", "H2O")
            == model.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert (
            model.fs.stream[0].default_material_balance_type()
            == MaterialBalanceType.componentPhase
        )
        assert (
            model.fs.stream[0].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )
        assert model.fs.stream[0].get_material_flow_basis() == MaterialFlowBasis.mass
        assert (
            model.fs.stream[0].get_enthalpy_flow_terms("Liq")
            == model.fs.stream[0].enth_flow["Liq"]
        )

        # NOTE: At this point, all property methods have been touched, so
        #       all will now be built

    @pytest.mark.unit
    def test_stats(self, coag_obj):
        model = coag_obj

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert degrees_of_freedom(model) == 6
        assert_units_consistent(model)

        # Fix some variables and verify it closes the system of equations
        model.fs.stream[0].temperature.fix(298)
        model.fs.stream[0].pressure.fix(101325)
        model.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
        model.fs.stream[0].flow_mass_phase_comp["Liq", "TSS"].fix(0.01)
        model.fs.stream[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.01)
        model.fs.stream[0].flow_mass_phase_comp["Liq", "Sludge"].fix(0.001)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, coag_obj):
        model = coag_obj

        # Set some scaling factors and look for 'bad' scaling
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TSS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "Sludge")
        )
        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(iscale.unscaled_constraints_generator(model))
        assert len(unscaled_constraint_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, coag_obj):
        model = coag_obj

        # call the 'calculate_state', which will call initialize and return results
        #       pass var_args as the fixed states from 'test_stats'
        args = {
            ("temperature", None): 298,
            ("pressure", None): 101325,
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.01,
            ("flow_mass_phase_comp", ("Liq", "TSS")): 0.01,
            ("flow_mass_phase_comp", ("Liq", "Sludge")): 0.001,
        }
        results = model.fs.stream.calculate_state(var_args=args)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, coag_obj):
        model = coag_obj

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, coag_obj):
        model = coag_obj

        assert value(
            model.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"]
        ) == pytest.approx(1, rel=1e-4)
        assert value(
            model.fs.stream[0].flow_mass_phase_comp["Liq", "TSS"]
        ) == pytest.approx(0.01, rel=1e-4)
        assert value(
            model.fs.stream[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(0.01, rel=1e-4)
        assert value(
            model.fs.stream[0].flow_mass_phase_comp["Liq", "Sludge"]
        ) == pytest.approx(0.001, rel=1e-4)

        assert value(
            model.fs.stream[0].mass_frac_phase_comp["Liq", "H2O"]
        ) == pytest.approx(0.979431, rel=1e-4)
        assert value(
            model.fs.stream[0].mass_frac_phase_comp["Liq", "TSS"]
        ) == pytest.approx(0.009794, rel=1e-4)
        assert value(
            model.fs.stream[0].mass_frac_phase_comp["Liq", "TDS"]
        ) == pytest.approx(0.009794, rel=1e-4)
        assert value(
            model.fs.stream[0].mass_frac_phase_comp["Liq", "Sludge"]
        ) == pytest.approx(0.0009794, rel=1e-4)

        assert value(model.fs.stream[0].dens_mass_phase["Liq"]) == pytest.approx(
            1013.95727, rel=1e-4
        )
        assert value(model.fs.stream[0].visc_d_phase["Liq"]) == pytest.approx(
            0.0008944, rel=1e-4
        )

        assert value(
            model.fs.stream[0].conc_mass_phase_comp["Liq", "H2O"]
        ) == pytest.approx(993.1021, rel=1e-4)
        assert value(
            model.fs.stream[0].conc_mass_phase_comp["Liq", "TSS"]
        ) == pytest.approx(9.931021, rel=1e-4)
        assert value(
            model.fs.stream[0].conc_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(9.931021, rel=1e-4)
        assert value(
            model.fs.stream[0].conc_mass_phase_comp["Liq", "Sludge"]
        ) == pytest.approx(0.9931021, rel=1e-4)


class TestCoagulationPropPackFailures:
    @pytest.fixture(scope="class")
    def coag_obj_fail(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.properties = CoagulationParameterBlock()

        return model

    @pytest.mark.unit
    def test_default_scaling(self, coag_obj_fail):
        model = coag_obj_fail

        model.fs.stream = model.fs.properties.build_state_block([0])

        # call scaling without setting defaults
        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(iscale.unscaled_constraints_generator(model))
        assert len(unscaled_constraint_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.unit
    def test_state_fail(self, coag_obj_fail):
        model = coag_obj_fail

        args = {
            ("temperature", None): 298,
            ("pressure", None): 101325,
            ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
            ("flow_mass_phase_comp", ("Liq", "TDS")): 0.01,
            ("flow_mass_phase_comp", ("Liq", "TSS")): 0.01,
            ("flow_mass_phase_comp", ("Liq", "Sludge")): 0.001,
        }
        with pytest.raises(ValueError):
            results = model.fs.stream.calculate_state(var_args=args)
