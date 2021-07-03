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

from pyomo.environ import (ConcreteModel,
                           Block,
                           Set,
                           Param,
                           Var,
                           Expression,
                           Constraint,
                           value,
                           units as pyunits,
                           SolverStatus)
from pyomo.util.check_units import assert_units_consistent, check_units_equivalent
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
from idaes.core.util.exceptions import (PropertyPackageError,
                                        PropertyNotSupportedError,
                                        ConfigurationError)
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
solver.options["nlp_scaling_method"] = "user-scaling"


# -----------------------------------------------------------------------------
class PropertyTestHarness(object):
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
        for v in m.fs.properties.component_objects(Var):
            assert v.is_fixed()

    @pytest.mark.unit
    def test_state_variables(self, frame_stateblock):
        m = frame_stateblock
        # state variables
        state_vars_dict = m.fs.stream[0].define_state_vars()
        assert len(state_vars_dict) >= 3  # minimum number of state variables: flow, temperature, pressure

        # check that state variables are in metadata
        metadata = m.fs.properties.get_metadata().properties
        for sv_name in state_vars_dict:
            assert sv_name in metadata

    @pytest.mark.unit
    def test_permanent_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties
        for v_name in metadata:
            if metadata[v_name]['method'] is None:
                assert hasattr(m.fs.stream[0], v_name)

    @pytest.mark.unit
    def test_on_demand_properties(self, frame_stateblock):
        m = frame_stateblock
        metadata = m.fs.properties.get_metadata().properties

        # check that properties are not built if not demanded
        for v_name in metadata:
            if metadata[v_name]['method'] is not None:
                assert not m.fs.stream[0].is_property_constructed(v_name)

        # check that properties are built if demanded
        for v_name in metadata:
            if metadata[v_name]['method'] is not None:
                assert hasattr(m.fs.stream[0], v_name)

    @pytest.mark.unit
    def test_stateblock_statistics(self, frame_stateblock):
        m = frame_stateblock
        blk = m.fs.stream[0]
        assert number_variables(blk) == \
               m._test_objs.stateblock_statistics['number_variables']
        assert number_total_constraints(blk) == \
               m._test_objs.stateblock_statistics['number_total_constraints']
        assert number_unused_variables(blk) == \
               m._test_objs.stateblock_statistics['number_unused_variables']
        assert degrees_of_freedom(blk) == \
               m._test_objs.stateblock_statistics['default_degrees_of_freedom']

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
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m.fs.stream[0]))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_default_initialization(self, frame_stateblock):
        m = frame_stateblock

        # initialize
        m.fs.stream.initialize()

        # check results
        for (v_name, ind), val in m._test_objs.default_solution.items():
            var = getattr(m.fs.stream[0], v_name)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise Exception(
                    "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
                    "but it has a value of {val_t}. \nPlease update default_solution dict in the "
                    "configure function that sets up the PropertyTestHarness".format(
                        v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))

    @pytest.mark.component
    def test_badly_scaled(self, frame_stateblock):
        m = frame_stateblock
        badly_scaled_var_list = list(badly_scaled_var_generator(m))
        assert len(badly_scaled_var_list) == 0

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

        # check results, properties_in and properties_out should be the same as default initialize
        for (v_name, ind), val in m._test_objs.default_solution.items():
            for sb in [m.fs.cv.properties_in[0], m.fs.cv.properties_out[0]]:
                if sb.is_property_constructed(v_name):
                    var = getattr(sb, v_name)  # get property if it was created
                else:
                    continue  # skip property if it was not created
                if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                    raise Exception(
                        "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
                        "but it has a value of {val_t}. \nPlease update default_solution dict in the "
                        "configure function that sets up the PropertyTestHarness".format(
                            v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))



class PropertyRegressionTest(object):
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
        assert results.solver.status == SolverStatus.ok

        # check results
        for (v_str, ind), val in self.regression_solution.items():
            var = getattr(m.fs.stream[0], v_str)
            if not pytest.approx(val, rel=1e-3) == value(var[ind]):
                raise Exception(
                    "Variable {v_str} with index {ind} is expected to have a value of {val} +/- 0.1%, but it "
                    "has a value of {val_t}. \nPlease update regression_solution in the configure function "
                    "that sets up the PropertyRegressionTest".format(
                        v_str=v_str, ind=ind, val=val, val_t=value(var[ind])))
