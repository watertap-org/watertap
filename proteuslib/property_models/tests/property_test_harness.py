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
from math import fabs

from pyomo.environ import (ConcreteModel,
                           Block,
                           Set,
                           Var,
                           Constraint,
                           value,
                           assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import (FlowsheetBlock,
                        ControlVolume0DBlock)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
class PropertyAttributeError(AttributeError):
    """
    ProteusLib exception for generic attribute errors arising from property package testing.
    """
    pass


class PropertyValueError(ValueError):
    """
    ProteusLib exception for generic value errors arising from property package testing.
    """
    pass


class PropertyRuntimeError(RuntimeError):
    """
    ProteusLib exception for generic runtime errors arising from property package testing.
    """
    pass


class PropertyTestHarness():
    def configure_class(self, m):
        self.configure()

        # attaching objects to model to carry through in pytest frame
        assert not hasattr(m, '_test_objs')
        m._test_objs = Block()
        m._test_objs.stateblock_statistics = self.stateblock_statistics
        m._test_objs.default_solution = self.default_solution

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
        pass

    @pytest.fixture(scope='class')
    def frame_stateblock(self):
        m = ConcreteModel()
        self.configure_class(m)

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], default=self.param_args)

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
                        v_name=v.name))

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
                    num_state_vars=len(state_vars_dict)))

        # check that state variables are in metadata
        metadata = m.fs.properties.get_metadata().properties
        for sv_name in state_vars_dict:
            if sv_name not in metadata:
                raise PropertyAttributeError(
                    "State variable {sv_name} is not included in the "
                    "metadata".format(sv_name=sv_name))

    @pytest.mark.unit
    def test_permanent_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties
        for v_name in metadata:
            if metadata[v_name]['method'] is None:
                if not hasattr(m.fs.stream[0], v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is included in the metadata but is not found "
                        "on the stateblock".format(v_name=v_name))

    @pytest.mark.unit
    def test_on_demand_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties

        # check that properties are not built if not demanded
        for v_name in metadata:
            if metadata[v_name]['method'] is not None:
                if m.fs.stream[0].is_property_constructed(v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was found "
                        "on the stateblock without being demanded".format(v_name=v_name))

        # check that properties are built if demanded
        for v_name in metadata:
            if metadata[v_name]['method'] is not None:
                if not hasattr(m.fs.stream[0], v_name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was not built "
                        "when demanded".format(v_name=v_name))

    @pytest.mark.unit
    def test_stateblock_statistics(self, frame_stateblock):
        m = frame_stateblock
        blk = m.fs.stream[0]
        stats = m._test_objs.stateblock_statistics
        if number_variables(blk) != stats['number_variables']:
            raise PropertyValueError(
                "The number of variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_variables(blk),
                    num_test=stats['number_variables']))
        if number_total_constraints(blk) != stats['number_total_constraints']:
            raise PropertyValueError(
                "The number of constraints were {num}, but {num_test} was "
                "expected ".format(
                    num=number_total_constraints(blk),
                    num_test=stats['number_total_constraints']))
        if number_unused_variables(blk) != stats['number_unused_variables']:
            raise PropertyValueError(
                "The number of unused variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_unused_variables(blk),
                    num_test=stats['number_unused_variables']))
        if degrees_of_freedom(blk) != stats['default_degrees_of_freedom']:
            raise PropertyValueError(
                "The number of degrees of freedom were {num}, but {num_test} was "
                "expected ".format(
                    num=degrees_of_freedom(blk),
                    num_test=stats['default_degrees_of_freedom']))

    @pytest.mark.unit
    def test_units_consistent(self, frame_stateblock):
        m = frame_stateblock
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
                    lst=unscaled_var_name_list))

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m.fs.stream[0]))
        if len(unscaled_constraint_list) != 0:
            unscaled_constraint_name_list = [c.name for c in unscaled_constraint_list]
            raise PropertyAttributeError(
                "The following constraint(s) are unscaled: {lst}".format(
                    lst=unscaled_constraint_name_list))

    @staticmethod
    def check_constraint_status(blk, tol=1e-8):
        c_violated_lst = []

        # implementation modified from Pyomo 5.7.3, pyomo.utils.infeasible.log_infeasible_constraints
        for c in blk.component_data_objects(
                ctype=Constraint, active=True, descend_into=True):
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
                    violated_lst=c_violated_lst))

    @pytest.mark.component
    def test_default_initialization(self, frame_stateblock):
        m = frame_stateblock

        # initialize
        m.fs.stream.initialize()

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.stream[0])

        # check results
        for (v_name, ind), val in m._test_objs.default_solution.items():
            var = getattr(m.fs.stream[0], v_name)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise PropertyValueError(
                    "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
                    "but it has a value of {val_t}. \nUpdate default_solution dict in the "
                    "configure function that sets up the PropertyTestHarness".format(
                        v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))

    @pytest.mark.component
    def test_badly_scaled(self, frame_stateblock):
        m = frame_stateblock
        badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
        if len(badly_scaled_var_list) != 0:
            lst = []
            for (var, val) in badly_scaled_var_list:
                lst.append((var.name, val))
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst))

    @pytest.mark.component
    def test_default_solve(self, frame_stateblock):
        m = frame_stateblock

        # fix state variables
        fix_state_vars(m.fs.stream)

        # solve model
        opt = get_solver()
        results = opt.solve(m.fs.stream[0])
        assert_optimal_termination(results)

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.stream[0])

        # check results
        for (v_name, ind), val in m._test_objs.default_solution.items():
            var = getattr(m.fs.stream[0], v_name)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise PropertyValueError(
                    "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
                    "but it has a value of {val_t}. \nUpdate default_solution dict in the "
                    "configure function that sets up the PropertyTestHarness".format(
                        v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))

    @pytest.fixture(scope='class')
    def frame_control_volume(self):
        m = ConcreteModel()
        self.configure_class(m)

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = self.prop_pack()

        m.fs.cv = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": m.fs.properties,
            "property_package_args": self.param_args})
        m.fs.cv.add_state_blocks(
            has_phase_equilibrium=False)
        m.fs.cv.add_material_balances()
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
        m.fs.cv.initialize()

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        self.check_constraint_status(m.fs.cv)

        # check results, properties_in and properties_out should be the same as default initialize
        for (v_name, ind), val in m._test_objs.default_solution.items():
            for sb in [m.fs.cv.properties_in[0], m.fs.cv.properties_out[0]]:
                if sb.is_property_constructed(v_name):
                    var = getattr(sb, v_name)  # get property if it was created
                else:
                    continue  # skip property if it was not created
                if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                    raise PropertyValueError(
                        "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
                        "but it has a value of {val_t}. \nUpdate default_solution dict in the "
                        "configure function that sets up the PropertyTestHarness".format(
                            v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))


class PropertyRegressionTest():
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
        pass

    @pytest.mark.component
    def test_property_regression(self):
        self.configure_class()

        # create model
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], default=self.param_args)

        # set default scaling
        for (v_str, ind), sf in self.scaling_args.items():
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        # set state variables
        for (v_str, ind), val in self.state_args.items():
            var = getattr(m.fs.stream[0], v_str)
            var[ind].fix(val)

        # touch all properties
        metadata = m.fs.properties.get_metadata().properties
        for v_str in metadata.keys():
            getattr(m.fs.stream[0], v_str)

        # scale model
        calculate_scaling_factors(m)

        # solve model
        opt = get_solver(self.solver, self.optarg)
        results = opt.solve(m)
        assert_optimal_termination(results)

        # check results
        for (v_str, ind), val in self.regression_solution.items():
            var = getattr(m.fs.stream[0], v_str)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise PropertyValueError(
                    "Variable {v_str} with index {ind} is expected to have a value of {val} +/- 0.1%, but it "
                    "has a value of {val_t}. \nUpdate regression_solution in the configure function "
                    "that sets up the PropertyRegressionTest".format(
                        v_str=v_str, ind=ind, val=val, val_t=value(var[ind])))

        # check if any variables are badly scaled
        badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
        if len(badly_scaled_var_list) != 0:
            lst = []
            for (var, val) in badly_scaled_var_list:
                lst.append((var.name, val))
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst))


class PropertyCalculateStateTest():
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
        pass

    @pytest.mark.component
    def test_calculate_state(self):
        self.configure_class()

        # create model
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = self.prop_pack()
        m.fs.stream = m.fs.properties.build_state_block([0], default=self.param_args)

        # set default scaling
        for (v_str, ind), sf in self.scaling_args.items():
            m.fs.properties.set_default_scaling(v_str, sf, index=ind)

        # touch all properties in var_args
        for (v_name, ind), val in self.var_args.items():
            getattr(m.fs.stream[0], v_name)

        # scale model
        calculate_scaling_factors(m)

        # calculate state
        results = m.fs.stream.calculate_state(var_args=self.var_args, solver=self.solver, optarg=self.optarg)
        assert_optimal_termination(results)

        # check results
        for (v_str, ind), val in self.state_solution.items():
            var = getattr(m.fs.stream[0], v_str)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise PropertyValueError(
                    "Variable {v_str} with index {ind} is expected to have a value of {val} +/- 0.1%, but it "
                    "has a value of {val_t}. \nUpdate state_solution in the configure function "
                    "that sets up the PropertyCalculateStateTest".format(
                        v_str=v_str, ind=ind, val=val, val_t=value(var[ind])))

        # check if any variables are badly scaled
        badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
        if len(badly_scaled_var_list) != 0:
            lst = []
            for (var, val) in badly_scaled_var_list:
                lst.append((var.name, val))
                print(var.name, var.value)
            raise PropertyValueError(
                "The following variable(s) are poorly scaled: {lst}".format(lst=lst))