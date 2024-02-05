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

from pyomo.environ import Block, assert_optimal_termination, ComponentMap, value
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
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
        self.unit_solutions = ComponentMap()

        # arguments for badly scaled variables
        self.default_large = 1e4
        self.default_small = 1e-3
        self.default_zero = 1e-10

        # arguments for solver tolerance
        self.default_absolute_tolerance = 1e-12
        self.default_relative_tolerance = 1e-6

        self.configure()
        blk = self.unit_model_block

        # attaching objects to model to carry through in pytest frame
        # TODO: Consider removing these objects and directly calling self
        assert not hasattr(blk, "_test_objs")
        blk._test_objs = Block()
        blk._test_objs.solver = self.solver
        blk._test_objs.optarg = self.optarg
        blk._test_objs.unit_solutions = self.unit_solutions

    def configure(self):
        """
        Placeholder method to allow user to setup test harness.

        The configure function must set the attributes:

        unit_model: pyomo unit model block (e.g. m.fs.unit), the block should
            have zero degrees of freedom, i.e. fully specified

        unit_solutions: dictionary of property values for the specified state variables
        """

    @pytest.fixture(scope="class")
    def frame_unit(self):
        self.configure_class()
        return self.unit_model_block

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
        badly_scaled_vars = list(
            iscale.badly_scaled_var_generator(
                blk,
                large=self.default_large,
                small=self.default_small,
                zero=self.default_zero,
            )
        )
        if badly_scaled_vars:
            lines = [
                f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}"
                for x in badly_scaled_vars
            ]
            msg = "One or more badly scaled variables found:\n"
            msg += "\n".join(lines)
            raise AssertionError(msg)

        assert_optimal_termination(results)

        # check results

        for var, val in solutions.items():
            comp_obj = None
            try:
                val = float(val)
            except:
                # expect the same API as pytest.approx
                comp_obj = val
                val = comp_obj.expected
            if comp_obj is None:
                comp_obj = pytest.approx(
                    val,
                    abs=self.default_absolute_tolerance,
                    rel=self.default_relative_tolerance,
                )
            if not comp_obj == value(var):
                raise AssertionError(f"{var}: Expected {val}, got {value(var)} instead")
