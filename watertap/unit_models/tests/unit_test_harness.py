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

from pyomo.environ import Block, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
class UnitAttributeError(AttributeError):
    """
    WaterTAP exception for generic attribute errors arising from unit model testing.
    """


class UnitValueError(ValueError):
    """
    WaterTAP exception for generic value errors arising from unit model testing.
    """


class UnitRuntimeError(RuntimeError):
    """
    WaterTAP exception for generic runtime errors arising from unit model testing.
    """


class UnitTestHarness:
    def configure_class(self):
        self.solver = None  # string for solver, if None use WaterTAP default
        self.optarg = (
            None  # dictionary for solver options, if None use WaterTAP default
        )

        self.configure()
        blk = self.unit_model_block

        # attaching objects to model to carry through in pytest frame
        assert not hasattr(blk, "_test_objs")
        blk._test_objs = Block()
        blk._test_objs.solver = self.solver
        blk._test_objs.optarg = self.optarg
        blk._test_objs.stateblock_statistics = self.unit_statistics
        blk._test_objs.unit_solutions = self.unit_solutions

    def configure(self):
        """
        Placeholder method to allow user to setup test harness.

        The configure function must set the attributes:

        unit_model: pyomo unit model block (e.g. m.fs.unit), the block should
            have zero degrees of freedom, i.e. fully specified

        unit_statistics: dictionary of model statistics
            {'number_config_args': VALUE,
             'number_variables': VALUE,
             'number_total_constraints': VALUE,
             'number_unused_variables': VALUE}

        unit_solutions: dictionary of property values for the specified state variables
            keys = (string name of variable, tuple index), values = value
        """
        pass

    @pytest.fixture(scope="class")
    def frame_unit(self):
        self.configure_class()
        return self.unit_model_block

    @pytest.mark.unit
    def test_unit_statistics(self, frame_unit):
        blk = frame_unit
        stats = blk._test_objs.stateblock_statistics

        if number_variables(blk) != stats["number_variables"]:
            raise UnitValueError(
                "The number of variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_variables(blk), num_test=stats["number_variables"]
                )
            )
        if number_total_constraints(blk) != stats["number_total_constraints"]:
            raise UnitValueError(
                "The number of constraints were {num}, but {num_test} was "
                "expected ".format(
                    num=number_total_constraints(blk),
                    num_test=stats["number_total_constraints"],
                )
            )
        if number_unused_variables(blk) != stats["number_unused_variables"]:
            raise UnitValueError(
                "The number of unused variables were {num}, but {num_test} was "
                "expected ".format(
                    num=number_unused_variables(blk),
                    num_test=stats["number_unused_variables"],
                )
            )

    @pytest.mark.unit
    def test_units_consistent(self, frame_unit):
        assert_units_consistent(frame_unit)

    @pytest.mark.unit
    def test_dof(self, frame_unit):
        if degrees_of_freedom(frame_unit) != 0:
            raise UnitAttributeError(
                "The unit has {dof} degrees of freedom when 0 is required."
                "".format(dof=degrees_of_freedom(frame_unit))
            )

    @pytest.mark.component
    def test_initialization(self, frame_unit):
        blk = frame_unit

        # initialize
        try:
            blk.initialize(solver=blk._test_objs.solver, optarg=blk._test_objs.optarg)
        except InitializationError:
            raise InitializationError(
                "The unit has failed to initialize successfully. Please check the output logs for more information."
            )

    @pytest.mark.component
    def test_unit_solutions(self, frame_unit):
        self.configure_class()
        blk = frame_unit
        solutions = blk._test_objs.unit_solutions

        # solve unit
        if blk._test_objs.solver is None:
            opt = get_solver()
        else:
            opt = get_solver(
                solver=blk._test_objs.solver, options=blk._test_objs.optarg
            )
        results = opt.solve(blk)

        # check solve
        try:
            assert len(list(iscale.badly_scaled_var_generator(blk, zero=1e-8))) == 0
        except AssertionError:
            badly_scaled_var_list = iscale.badly_scaled_var_generator(blk, zero=1e-8)
            for x in badly_scaled_var_list:
                print(
                    f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}"
                )
            raise AssertionError(
                "Check the output logs for variables that need to be re-scaled"
            )
        assert_optimal_termination(results)

        # check results
        for key, val in solutions.items():
            pyo_obj = blk.find_component(key)
            pyo_obj_value = pyo_obj()
            try:
                if val == 0:
                    assert pytest.approx(pyo_obj_value, abs=1e-08) == val
                else:
                    assert pytest.approx(pyo_obj_value, rel=1e-03) == val
            except:
                raise AssertionError(
                    f"{pyo_obj}: Expected {val}, but got {pyo_obj_value}"
                )
