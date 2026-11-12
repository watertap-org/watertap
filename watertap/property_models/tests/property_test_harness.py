#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
from math import fabs
import copy

from pyomo.environ import (
    ConcreteModel,
    Block,
    Set,
    Var,
    Constraint,
    value,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, ControlVolume0DBlock, MaterialBalanceType
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
from idaes.core.util.testing import assert_solution_equivalent
from watertap.core.solvers import get_solver


# -----------------------------------------------------------------------------
class PropertyAttributeError(AttributeError):
    """
    WaterTAP exception for generic attribute errors arising from property package testing.
    """


class PropertyValueError(ValueError):
    """
    WaterTAP exception for generic value errors arising from property package testing.
    """


class PropertyRuntimeError(RuntimeError):
    """
    WaterTAP exception for generic runtime errors arising from property package testing.
    """


def _legacy_testing_translator(blk, default_solution, check_property_constructed=False):
    """
    Translates tests written using WaterTAP's test harness and the new standard
    testing function assert_solution_equivalent.

    Args:
        blk: Block with variables whose values are being tested
        default_solution: Dictionary of cached values being tested
        check_property_constructed: Whether to check if a property is
            constructed before adding it to be tested.

    Returns:
        None
    """
    expected_results = {}
    for (v_name, ind), val in default_solution.items():
        if check_property_constructed and not blk.is_property_constructed(v_name):
            continue
        if v_name not in expected_results:
            expected_results[v_name] = {}
        # TODO revise behavior for near-zero but not exactly zero values
        if val == 0:
            expected_results[v_name][ind] = (val, None, 1e-8)
        else:
            expected_results[v_name][ind] = (val, 1e-3, None)

    assert_solution_equivalent(blk, expected_results)


class PropertyTestHarness:

    def configure_class(self, m):
        self.expected_on_demand = set()
        self.configure()

        # attaching objects to model to carry through in pytest frame
        assert not hasattr(m, "_test_objs")
        m._test_objs = Block()
        m._test_objs.stateblock_statistics = self.stateblock_statistics
        m._test_objs.default_solution = self.default_solution
        m._test_objs.expected_on_demand = self.expected_on_demand

    def configure(self):
        """
        Placeholder method to allow user to setup test harness.

        The configure function must set the attributes:

        prop_pack: property package parameter block

        param_args: dictionary for property parameter arguments

        scaling_args: dictionary for standard scaling arguments
            keys = (string name of variable, tuple index), values = scaling factor

        stateblock_statistics: dictionary of model statistics
            {'number_variables': VALUE,
             'number_total_constraints': VALUE,
             'number_unused_variables': VALUE,
             'default_degrees_of_freedom': VALUE}

        default_solution: dictionary of property values at default initialized state variables
            keys = (string name of variable, tuple index), values = value
        """

    @pytest.fixture(scope="class")
    @classmethod
    def frame_stateblock(cls):
        m = ConcreteModel()
        self = cls()
        self.configure_class(m)

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], **self.param_args)

        for (v_str, ind), sf in self.scaling_args.items():
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        return m

    @pytest.mark.unit
    def test_components_phases(self, frame_stateblock):
        m = frame_stateblock
        assert isinstance(m.fs.properties.component_list, Set)
        assert isinstance(m.fs.properties.phase_list, Set)
        assert isinstance(m.fs.properties._phase_component_set, Set)

    @pytest.mark.unit
    def test_parameters(self, frame_stateblock):
        m = frame_stateblock
        # test that the parameter variables are fixed
        for v in m.fs.properties.component_data_objects(Var):
            if not v.is_fixed():
                raise PropertyRuntimeError(
                    "Variable {v_name} is unfixed, but all variables on the "
                    "property parameter block must be fixed by default".format(
                        v_name=v.name
                    )
                )

    @pytest.mark.unit
    def test_state_variables(self, frame_stateblock):
        m = frame_stateblock
        # state variables
        state_vars_dict = m.fs.stream[0].define_state_vars()
        # check minimum number of state variables (i.e. flow, temperature, pressure)
        if len(state_vars_dict) < 3:
            raise PropertyValueError(
                "Property package has {num_state_vars} state variables, which is less "
                "than the minimum of 3 (i.e. flow, temperature, and pressure)".format(
                    num_state_vars=len(state_vars_dict)
                )
            )

        # check that state variables are in metadata
        metadata = m.fs.properties.get_metadata().properties
        metadata_strs = [v.name for v in metadata.list_supported_properties()]
        for sv_name in state_vars_dict:
            if sv_name not in metadata_strs:
                raise PropertyAttributeError(
                    "State variable {sv_name} is not included in the "
                    "metadata".format(sv_name=sv_name)
                )

    @pytest.mark.unit
    def test_permanent_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties
        for v in metadata.list_supported_properties():
            if metadata[v.name].method is None:
                if not hasattr(m.fs.stream[0], v.name):
                    raise PropertyAttributeError(
                        "Property {v_name} is included in the metadata but is not found "
                        "on the stateblock".format(v_name=v.name)
                    )

    @pytest.mark.unit
    def test_on_demand_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties

        # check that properties are not built if not demanded
        for v in metadata.list_supported_properties():
            if metadata[v.name].method is not None:
                if v.name in m._test_objs.expected_on_demand:
                    continue
                if m.fs.stream[0].is_property_constructed(v.name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was found "
                        "on the stateblock without being demanded".format(v_name=v.name)
                    )

        # check that properties are built if demanded
        for v in metadata.list_supported_properties():
            if metadata[v.name].method is not None:
                if not hasattr(m.fs.stream[0], v.name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was not built "
                        "when demanded".format(v_name=v.name)
                    )

    @pytest.mark.unit
    def test_stateblock_statistics(self, frame_stateblock):
        m = frame_stateblock
        blk = m.fs.stream[0]
        stats = m._test_objs.stateblock_statistics
        if number_variables(blk) != stats["number_variables"]:
            raise PropertyValueError(
                "The number of variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_variables(blk), num_test=stats["number_variables"]
                )
            )
        if number_total_constraints(blk) != stats["number_total_constraints"]:
            raise PropertyValueError(
                "The number of constraints were {num}, but {num_test} was "
                "expected ".format(
                    num=number_total_constraints(blk),
                    num_test=stats["number_total_constraints"],
                )
            )
        if number_unused_variables(blk) != stats["number_unused_variables"]:
            raise PropertyValueError(
                "The number of unused variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_unused_variables(blk),
                    num_test=stats["number_unused_variables"],
                )
            )
        if degrees_of_freedom(blk) != stats["default_degrees_of_freedom"]:
            raise PropertyValueError(
                "The number of degrees of freedom were {num}, but {num_test} was "
                "expected ".format(
                    num=degrees_of_freedom(blk),
                    num_test=stats["default_degrees_of_freedom"],
                )
            )

    @pytest.mark.unit
    def test_units_consistent(self, frame_stateblock):
        m = frame_stateblock
        # TODO replace with structural issues
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_scaling(self, frame_stateblock):
        m = frame_stateblock
        calculate_scaling_factors(m.fs.stream[0])

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.stream[0]))
        if len(unscaled_var_list) != 0:
            unscaled_var_name_list = [v.name for v in unscaled_var_list]
            raise PropertyAttributeError(
                "The following variable(s) are unscaled: {lst}".format(
                    lst=unscaled_var_name_list
                )
            )

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m.fs.stream[0]))
        if len(unscaled_constraint_list) != 0:
            unscaled_constraint_name_list = [c.name for c in unscaled_constraint_list]
            raise PropertyAttributeError(
                "The following constraint(s) are unscaled: {lst}".format(
                    lst=unscaled_constraint_name_list
                )
            )

    @staticmethod
    def check_constraint_status(blk, tol=1e-8):
        # TODO replace with assert_no_numerical_issues()
        c_violated_lst = []

        # implementation modified from Pyomo 5.7.3, pyomo.utils.infeasible.log_infeasible_constraints
        for c in blk.component_data_objects(
            ctype=Constraint, active=True, descend_into=True
        ):
            c_body_value = value(c.body, exception=False)
            c_lb_value = value(c.lower, exception=False)
            c_ub_value = value(c.upper, exception=False)

            c_undefined = False
            equality_violated = False
            lb_violated = False
            ub_violated = False

            if c_body_value is None:
                c_undefined = True
            else:
                if c.equality:
                    if fabs(c_lb_value - c_body_value) >= tol:
                        equality_violated = True
                else:
                    if c.has_lb() and c_lb_value - c_body_value >= tol:
                        lb_violated = True
                    if c.has_ub() and c_body_value - c_ub_value >= tol:
                        ub_violated = True

            if not any((c_undefined, equality_violated, lb_violated, ub_violated)):
                continue
            else:
                c_violated_lst.append(c.name)

        if len(c_violated_lst) > 0:
            raise PropertyRuntimeError(
                "Default initialization did not converge, the following constraint(s) "
                "are undefined or violate the equality or inequality: {violated_lst}".format(
                    violated_lst=c_violated_lst
                )
            )

    @pytest.mark.component
    def test_default_initialization(self, frame_stateblock):
        m = frame_stateblock

        # initialize
        m.fs.stream.initialize()

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.stream[0])

        # check results
        _legacy_testing_translator(m.fs.stream[0], m._test_objs.default_solution)

    @pytest.mark.component
    def test_badly_scaled(self, frame_stateblock):
        m = frame_stateblock
        badly_scaled_var_list = list(
            badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8)
        )
        if len(badly_scaled_var_list) != 0:
            lst = []
            for var, val in badly_scaled_var_list:
                lst.append((var.name, val))
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst)
            )

    @pytest.mark.component
    def test_default_solve(self, frame_stateblock):
        m = frame_stateblock

        m.fs.stream.fix_initialization_states()

        # solve model
        opt = get_solver()
        results = opt.solve(m.fs.stream[0])
        assert_optimal_termination(results)

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.stream[0])

        # check results
        _legacy_testing_translator(m.fs.stream[0], m._test_objs.default_solution)

    @pytest.mark.unit
    def test_list_and_print_properties(self, frame_stateblock):
        m = frame_stateblock
        m.fs.properties.list_properties()
        m.fs.properties.print_properties()

    @pytest.fixture(scope="class")
    @classmethod
    def frame_control_volume(cls):
        m = ConcreteModel()
        self = cls()
        self.configure_class(m)

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = self.prop_pack()

        param_args = copy.copy(self.param_args)
        # Defined state is set by the control volume
        # so we can ignore it here.
        param_args.pop("defined_state", None)

        m.fs.cv = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=m.fs.properties,
            property_package_args=self.param_args,
        )
        m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
        # Since we have no VLE, we need phase-component balances
        # to have a square problem
        m.fs.cv.add_material_balances(balance_type=MaterialBalanceType.componentPhase)
        m.fs.cv.add_energy_balances()
        m.fs.cv.add_momentum_balances()

        for (v_str, ind), sf in self.scaling_args.items():
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        return m

    @pytest.mark.component
    def test_property_control_volume(self, frame_control_volume):
        m = frame_control_volume

        # scale model
        calculate_scaling_factors(m)

        # initialize control volume
        flags = m.fs.cv.initialize(hold_state=True)
        # Control volume is not solved by the initialization
        solver_obj = get_solver()
        results = solver_obj.solve(m.fs.cv)
        assert_optimal_termination(results)

        # Unfix inlet state vars
        m.fs.cv.release_state(flags)

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.cv)

        for blk in [m.fs.cv.properties_in[0], m.fs.cv.properties_out[0]]:
            try:
                _legacy_testing_translator(
                    blk,
                    m._test_objs.default_solution,
                    check_property_constructed=True,
                )
            except AssertionError as err:
                raise AssertionError(
                    f"Values in {blk.name} do not match expected values."
                ) from err


class PropertyRegressionTest:
    def configure_class(self):
        self.solver = None  # string for solver, if None use IDAES default
        self.optarg = None  # dictionary for solver options, if None use IDAES default
        self.configure()

    def configure(self):
        """
        Placeholder method to allow user to setup the property regression test.

        The configure function must set the attributes:

        prop_pack: property package parameter block

        param_args: dictionary for property parameter arguments

        solver: string name for solver, if not provided or None will use IDAES default

        optarg: dictionary of solver options, if not provided or None will use IDAES default

        scaling_args: dictionary of scaling arguments
            keys = (string name of variable, tuple index), values = scaling factor

        state_args: dictionary of specified state variables
            keys = (string name of variable, tuple index), values = value

        regression_solution: dictionary of property values for the specified state variables
            keys = (string name of variable, tuple index), values = value
        """

    @pytest.mark.component
    def test_property_regression(self):
        self.configure_class()

        # create model
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], **self.param_args)

        # set default scaling
        for (v_str, ind), sf in self.scaling_args.items():
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        # set state variables
        for (v_str, ind), val in self.state_args.items():
            var = getattr(m.fs.stream[0], v_str)
            var[ind].fix(val)

        # touch all properties
        metadata = m.fs.properties.get_metadata().properties
        for p in metadata.list_supported_properties():
            getattr(m.fs.stream[0], p.name)

        # scale model
        calculate_scaling_factors(m)

        # solve model
        opt = get_solver(self.solver, self.optarg)
        results = opt.solve(m)
        assert_optimal_termination(results)

        # check results
        _legacy_testing_translator(
            m.fs.stream[0], default_solution=self.regression_solution
        )

        # check if any variables are badly scaled
        lst = []
        for var, val in badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8):
            lst.append((var.name, val))
            print(var.name, var.value)
        if lst:
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst)
            )


class PropertyCalculateStateTest:
    def configure_class(self):
        self.solver = None  # string for solver, if None use IDAES default
        self.optarg = None  # dictionary for solver options, if None use IDAES default
        self.configure()

    def configure(self):
        """
        Placeholder method to allow user to setup the property regression test.

        The configure function must set the attributes:

        prop_pack: property package parameter block

        param_args: dictionary for property parameter arguments

        solver: string name for solver, if not provided or None will use IDAES default

        optarg: dictionary of solver options, if not provided or None will use IDAES default

        scaling_args: dictionary of scaling arguments
            keys = (string name of variable, tuple index), values = scaling factor

        var_args: dictionary of specified variables
            keys = (string name of variable, tuple index), values = value

        state_solution: dictionary of the values for the state variables
            keys = (string name of variable, tuple index), values = value
        """

    @pytest.mark.component
    def test_calculate_state(self):
        self.configure_class()

        # create model
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], **self.param_args)
        # set default scaling
        for (v_str, ind), sf in self.scaling_args.items():
            if not isinstance(ind, str):
                try:
                    l = len(ind)
                except TypeError:
                    # Not a tuple
                    l = None
                if l is not None and l == 1:
                    # Tuple of length 1, set_default_scaling
                    # expects a string/number, not a tuple
                    # with one element
                    ind = ind[0]
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        # touch all properties in scaling
        for (v_name, ind), val in self.scaling_args.items():
            getattr(m.fs.stream[0], v_name)

        # touch all properties in var_args
        for (v_name, ind), val in self.var_args.items():
            getattr(m.fs.stream[0], v_name)

        # scale model
        calculate_scaling_factors(m)

        # calculate state
        results = m.fs.stream.calculate_state(
            var_args=self.var_args, solver=self.solver, optarg=self.optarg
        )
        assert_optimal_termination(results)

        # check results
        _legacy_testing_translator(m.fs.stream[0], default_solution=self.state_solution)

        # check if any variables are badly scaled
        lst = []
        for var, val in badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8):
            lst.append((var.name, val))
            print(var.name, var.value)
        if lst:
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst)
            )
