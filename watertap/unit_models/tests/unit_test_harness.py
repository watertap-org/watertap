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
class UnitAttributeError(AttributeError):
    """
    ProteusLib exception for generic attribute errors arising from unit model testing.
    """
    pass


class UnitValueError(ValueError):
    """
    ProteusLib exception for generic value errors arising from unit model testing.
    """
    pass


class UnitRuntimeError(RuntimeError):
    """
    ProteusLib exception for generic runtime errors arising from unit model testing.
    """
    pass


class UnitTestHarness():
    def configure_class(self):
        self.solver = None  # string for solver, if None use IDAES default
        self.optarg = None  # dictionary for solver options, if None use IDAES default

        self.configure()
        blk = self.unit_model

        # attaching objects to model to carry through in pytest frame
        assert not hasattr(blk, '_test_objs')
        blk._test_objs = Block()
        blk._test_objs.solver = self.solver
        blk._test_objs.optarg = self.optarg
        blk._test_objs.stateblock_statistics = self.unit_statistics
        blk._test_objs.unit_solution = self.unit_solution

    def configure(self):
        """
        Placeholder method to allow user to setup test harness.

        The configure function must set the attributes:

        unit_model: pyomo unit model block (e.g. m.fs.unit), the block should have:
                    1) all variables and constraints scaled
                    2) zero degrees of freedom, i.e. fully specified

        unit_statistics: dictionary of model statistics
            {'number_config_args': VALUE,
             'number_variables': VALUE,
             'number_total_constraints': VALUE,
             'number_unused_variables': VALUE}

        unit_solution: dictionary of property values for the specified state variables
            keys = (string name of variable, tuple index), values = value
        """
        pass

    @pytest.fixture(scope='class')
    def frame_unit(self):
        self.configure_class()
        return self.unit_model

    @pytest.mark.unit
    def test_unit_statistics(self, frame_unit):
        blk = frame_unit
        stats = blk._test_objs.stateblock_statistics

        if number_variables(blk) != stats['number_variables']:
            raise UnitValueError(
                "The number of variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_variables(blk),
                    num_test=stats['number_variables']))
        if number_total_constraints(blk) != stats['number_total_constraints']:
            raise UnitValueError(
                "The number of constraints were {num}, but {num_test} was "
                "expected ".format(
                    num=number_total_constraints(blk),
                    num_test=stats['number_total_constraints']))
        if number_unused_variables(blk) != stats['number_unused_variables']:
            raise UnitValueError(
                "The number of unused variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_unused_variables(blk),
                    num_test=stats['number_unused_variables']))

    @pytest.mark.unit
    def test_units_consistent(self, frame_unit):
        assert_units_consistent(frame_unit)

    @pytest.mark.unit  # TODO: consider getting rid of this
    def test_scaling(self, frame_unit):
        blk = frame_unit

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(blk))
        if len(unscaled_var_list) != 0:
            unscaled_var_name_list = [v.name for v in unscaled_var_list]
            raise UnitAttributeError(
                "The following variable(s) are unscaled: {lst}".format(
                    lst=unscaled_var_name_list))

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(blk))
        if len(unscaled_constraint_list) != 0:
            unscaled_constraint_name_list = [c.name for c in unscaled_constraint_list]
            raise UnitAttributeError(
                "The following constraint(s) are unscaled: {lst}".format(
                    lst=unscaled_constraint_name_list))

    @pytest.mark.unit
    def test_dof(self, frame_unit):
        if degrees_of_freedom(frame_unit) != 0:  # TODO: change to check_dof from WaterTAP utilities
            raise UnitAttributeError(
                "The unit has {dof} degrees of freedom when 0 is required."
                "".format(dof=degrees_of_freedom(frame_unit)))

    @staticmethod
    def check_constraint_status(blk, tol=1e-7):
        c_violated_lst = []

        # implementation modified from Pyomo 5.7.3, pyomo.util.infeasible.log_infeasible_constraints
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
            raise UnitRuntimeError(
                "Default initialization did not converge, the following constraint(s) "
                "are undefined or violate the equality or inequality: {violated_lst}".format(
                    violated_lst=c_violated_lst))

    @pytest.mark.component
    def test_initialization(self, frame_unit):
        blk = frame_unit

        # initialize
        blk.initialize(solver=blk._test_objs.solver, optarg=blk._test_objs.optarg)

        # check convergence
        # TODO: update this when IDAES API is updated to return solver status for initialize()
        blk.permeate_side.eq_mass_frac_permeate_io.display()
        self.check_constraint_status(blk)


    @pytest.mark.component
    def test_solve(self, frame_unit, capsys):
        blk = frame_unit

        # solve unit
        opt = get_solver(solver=blk._test_objs.solver, options=blk._test_objs.optarg)
        results = opt.solve(blk)

        # check solve
        assert_optimal_termination(results)

        # check results
        blk.report()
        captured = capsys.readouterr()
        report_out = \
"""
====================================================================================
Unit : unit                                                                Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                                            : Value       : Fixed : Bounds
                                     Membrane Area :      50.000 :  True : (0.1, 1000.0)
                                   Membrane Length :      20.000 :  True : (0.1, 500.0)
                                    Membrane Width :      2.5000 : False : (0.1, 500.0)
                   NaCl Concentration @Inlet,Bulk  :      35.751 : False : (0.001, 2000.0)
     NaCl Concentration @Inlet,Membrane-Interface  :      43.381 : False : (0.001, 2000.0)
                  NaCl Concentration @Outlet,Bulk  :      46.102 : False : (0.001, 2000.0)
    NaCl Concentration @Outlet,Membrane-Interface  :      50.832 : False : (0.001, 2000.0)
                      NaCl Permeate Concentration  :     0.36827 : False : (0.001, 2000.0)
        Osmotic Pressure @Inlet,Membrane-Interface :  3.4820e+06 : False : (500.0, 50000000.0)
                     Osmotic Pressure @Outlet,Bulk :  3.7086e+06 : False : (500.0, 50000000.0)
      Osmotic Pressure @Outlet,Membrane-Interface  :  4.1055e+06 : False : (500.0, 50000000.0)
                                   Pressure Change : -1.7701e+05 : False : (None, None)
                            Reynolds Number @Inlet :      473.82 : False : (10, 5000.0)
                           Reynolds Number @Outlet :      362.00 : False : (10, 5000.0)
                        Solvent Mass Recovery Rate :     0.22864 : False : (0.01, 0.999999)
                                   Velocity @Inlet :     0.23035 : False : (0.01, 5)
                                  Velocity @Outlet :     0.17821 : False : (0.01, 5)
                        Volumetric Flowrate @Inlet :  0.00097899 : False : (1e-08, None)
                       Volumetric Flowrate @Outlet :  0.00075741 : False : (1e-08, None)
                          Volumetric Recovery Rate :     0.22653 : False : (0.01, 0.999999)

------------------------------------------------------------------------------------
    Stream Table
                                          Feed Inlet  Feed Outlet  Permeate Outlet
    flow_mass_phase_comp ('Liq', 'H2O')      0.96500     0.74436        0.22064   
    flow_mass_phase_comp ('Liq', 'NaCl')    0.035000    0.034918     8.1671e-05   
    temperature                               298.15      298.15         298.15   
    pressure                              5.0000e+06  4.8230e+06     1.0132e+05   
====================================================================================
"""
        # print(len(captured.out), len(report_out))
        # print(repr(captured))
        # print('***captured.out***')
        # print(captured.out)
        # print('***report_out***')
        # print(report_out)
        print('***repr captured***')
        print(repr(captured))
        print('***repr captured.err***')
        print(repr(captured.err))
        print('***repr captured.out***')
        print(repr(captured.out))
        print('***repr report_out***')
        print(repr(report_out))

        # differences = difflib.ndiff(captured.out, report_out)
        #
        # for difference in differences:
        #     print(difference)
        import difflib
        cases = [(captured.out, report_out)]
        for a, b in cases:
            print('{} => {}'.format(a, b))
            for i, s in enumerate(difflib.ndiff(a, b)):
                if s[0] == ' ':
                    continue
                elif s[0] == '-':
                    print(u'Delete "{}" from position {}'.format(s[-1], i))
                elif s[0] == '+':
                    print(u'Add "{}" to position {}'.format(s[-1], i))
                print(a[i-50:i])
            print()

        assert captured.out == report_out

#         assert captured.out == \
# """
# ====================================================================================
# Unit : unit                                                                Time: 0.0
# ------------------------------------------------------------------------------------
#     Unit Performance
#
#     Variables:
#
#     Key                                            : Value       : Fixed : Bounds
#                                      Membrane Area :      50.000 :  True : (0.1, 1000.0)
#                                    Membrane Length :      20.000 :  True : (0.1, 500.0)
#                                     Membrane Width :      2.5000 : False : (0.1, 500.0)
#                    NaCl Concentration @Inlet,Bulk  :      35.751 : False : (0.001, 2000.0)
#      NaCl Concentration @Inlet,Membrane-Interface  :      43.381 : False : (0.001, 2000.0)
#                   NaCl Concentration @Outlet,Bulk  :      46.102 : False : (0.001, 2000.0)
#     NaCl Concentration @Outlet,Membrane-Interface  :      50.832 : False : (0.001, 2000.0)
#                       NaCl Permeate Concentration  :     0.36827 : False : (0.001, 2000.0)
#         Osmotic Pressure @Inlet,Membrane-Interface :  3.4820e+06 : False : (500.0, 50000000.0)
#                      Osmotic Pressure @Outlet,Bulk :  3.7086e+06 : False : (500.0, 50000000.0)
#       Osmotic Pressure @Outlet,Membrane-Interface  :  4.1055e+06 : False : (500.0, 50000000.0)
#                                    Pressure Change : -1.7701e+05 : False : (None, None)
#                             Reynolds Number @Inlet :      473.82 : False : (10, 5000.0)
#                            Reynolds Number @Outlet :      362.00 : False : (10, 5000.0)
#                         Solvent Mass Recovery Rate :     0.22864 : False : (0.01, 0.999999)
#                                    Velocity @Inlet :     0.23035 : False : (0.01, 5)
#                                   Velocity @Outlet :     0.17821 : False : (0.01, 5)
#                         Volumetric Flowrate @Inlet :  0.00097899 : False : (1e-08, None)
#                        Volumetric Flowrate @Outlet :  0.00075741 : False : (1e-08, None)
#                           Volumetric Recovery Rate :     0.22653 : False : (0.01, 0.999999)
#
# ------------------------------------------------------------------------------------
#     Stream Table
#                                           Feed Inlet  Feed Outlet  Permeate Outlet
#     flow_mass_phase_comp ('Liq', 'H2O')      0.96500     0.74436        0.22064
#     flow_mass_phase_comp ('Liq', 'NaCl')    0.035000    0.034918     8.1671e-05
#     temperature                               298.15      298.15         298.15
#     pressure                              5.0000e+06  4.8230e+06     1.0132e+05
# ====================================================================================
# """

#
#         # check results
#         for (v_name, ind), val in blk._test_objs.unit_solution.items():
#             var = getattr(blk, v_name)
#             if not pytest.approx(val, rel=1e-3) == value(var[ind]):
#                 raise UnitValueError(
#                     "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
#                     "but it has a value of {val_t}. \nUpdate unit_solution dict in the "
#                     "configure function that sets up the UnitTestHarness".format(
#                         v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))
#
#     @pytest.mark.component
#     def test_badly_scaled(self, frame_stateblock):
#         m = frame_stateblock
#         badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
#         if len(badly_scaled_var_list) != 0:
#             lst = []
#             for (var, val) in badly_scaled_var_list:
#                 lst.append((var.name, val))
#             raise PropertyValueError(
#                 "The following variable(s) are poorly scaled: {lst}".format(lst=lst))
#
#     @pytest.mark.component
#     def test_default_solve(self, frame_stateblock):
#         m = frame_stateblock
#
#         # fix state variables
#         fix_state_vars(m.fs.stream)
#
#         # solve model
#         opt = get_solver()
#         results = opt.solve(m.fs.stream[0])
#         assert results.solver.status == SolverStatus.ok
#         assert results.solver.termination_condition == TerminationCondition.optimal
#
#         # check convergence
#         # TODO: update this when IDAES API is updated to return solver status for initialize()
#         self.check_constraint_status(m.fs.stream[0])
#
#         # check results
#         for (v_name, ind), val in m._test_objs.default_solution.items():
#             var = getattr(m.fs.stream[0], v_name)
#             if not pytest.approx(val, rel=1e-3) == value(var[ind]):
#                 raise PropertyValueError(
#                     "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
#                     "but it has a value of {val_t}. \nUpdate default_solution dict in the "
#                     "configure function that sets up the PropertyTestHarness".format(
#                         v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))
#
#     @pytest.fixture(scope='class')
#     def frame_control_volume(self):
#         m = ConcreteModel()
#         self.configure_class(m)
#
#         m.fs = FlowsheetBlock(default={"dynamic": False})
#         m.fs.properties = self.prop_pack()
#
#         m.fs.cv = ControlVolume0DBlock(default={
#             "dynamic": False,
#             "has_holdup": False,
#             "property_package": m.fs.properties,
#             "property_package_args": self.param_args})
#         m.fs.cv.add_state_blocks(
#             has_phase_equilibrium=False)
#         m.fs.cv.add_material_balances()
#         m.fs.cv.add_energy_balances()
#         m.fs.cv.add_momentum_balances()
#
#         for (v_str, ind), sf in self.scaling_args.items():
#             m.fs.properties.set_default_scaling(v_str, sf, index=ind)
#
#         return m
#
#     @pytest.mark.component
#     def test_property_control_volume(self, frame_control_volume):
#         m = frame_control_volume
#
#         # scale model
#         calculate_scaling_factors(m)
#
#         # initialize control volume
#         m.fs.cv.initialize()
#
#         # check convergence
#         # TODO: update this when IDAES API is updated to return solver status for initialize()
#         self.check_constraint_status(m.fs.cv)
#
#         # check results, properties_in and properties_out should be the same as default initialize
#         for (v_name, ind), val in m._test_objs.default_solution.items():
#             for sb in [m.fs.cv.properties_in[0], m.fs.cv.properties_out[0]]:
#                 if sb.is_property_constructed(v_name):
#                     var = getattr(sb, v_name)  # get property if it was created
#                 else:
#                     continue  # skip property if it was not created
#                 if not pytest.approx(val, rel=1e-3) == value(var[ind]):
#                     raise PropertyValueError(
#                         "Variable {v_name} with index {ind} is expected to have a value of {val} +/- 0.1%, "
#                         "but it has a value of {val_t}. \nUpdate default_solution dict in the "
#                         "configure function that sets up the PropertyTestHarness".format(
#                             v_name=v_name, ind=ind, val=val, val_t=value(var[ind])))
#
#
# class PropertyRegressionTest():
#     def configure_class(self):
#         self.solver = None  # string for solver, if None use IDAES default
#         self.optarg = None  # dictionary for solver options, if None use IDAES default
#         self.configure()
#
#     def configure(self):
#         """
#         Placeholder method to allow user to setup the property regression test.
#
#         The configure function must set the attributes:
#
#         prop_pack: property package parameter block
#
#         param_args: dictionary for property parameter arguments
#
#         solver: string name for solver, if not provided or None will use IDAES default
#
#         optarg: dictionary of solver options, if not provided or None will use IDAES default
#
#         scaling_args: dictionary of scaling arguments
#             keys = (string name of variable, tuple index), values = scaling factor
#
#         state_args: dictionary of specified state variables
#             keys = (string name of variable, tuple index), values = value
#
#         regression_solution: dictionary of property values for the specified state variables
#             keys = (string name of variable, tuple index), values = value
#         """
#         pass
#
#     @pytest.mark.component
#     def test_property_regression(self):
#         self.configure_class()
#
#         # create model
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(default={"dynamic": False})
#         m.fs.properties = self.prop_pack()
#         m.fs.stream = m.fs.properties.build_state_block([0], default=self.param_args)
#
#         # set default scaling
#         for (v_str, ind), sf in self.scaling_args.items():
#             m.fs.properties.set_default_scaling(v_str, sf, index=ind)
#
#         # set state variables
#         for (v_str, ind), val in self.state_args.items():
#             var = getattr(m.fs.stream[0], v_str)
#             var[ind].fix(val)
#
#         # touch all properties
#         metadata = m.fs.properties.get_metadata().properties
#         for v_str in metadata.keys():
#             getattr(m.fs.stream[0], v_str)
#
#         # scale model
#         calculate_scaling_factors(m)
#
#         # solve model
#         opt = get_solver(self.solver, self.optarg)
#         results = opt.solve(m)
#         assert results.solver.status == SolverStatus.ok
#         assert results.solver.termination_condition == TerminationCondition.optimal
#
#         # check results
#         for (v_str, ind), val in self.regression_solution.items():
#             var = getattr(m.fs.stream[0], v_str)
#             if not pytest.approx(val, rel=1e-3) == value(var[ind]):
#                 raise PropertyValueError(
#                     "Variable {v_str} with index {ind} is expected to have a value of {val} +/- 0.1%, but it "
#                     "has a value of {val_t}. \nUpdate regression_solution in the configure function "
#                     "that sets up the PropertyRegressionTest".format(
#                         v_str=v_str, ind=ind, val=val, val_t=value(var[ind])))
#
#         # check if any variables are badly scaled
#         badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
#         if len(badly_scaled_var_list) != 0:
#             lst = []
#             for (var, val) in badly_scaled_var_list:
#                 lst.append((var.name, val))
#             raise PropertyValueError(
#                 "The following variable(s) are badly scaled: {lst}".format(lst=lst))
