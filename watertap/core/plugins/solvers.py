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

# The function _constraint_autoscale_large_jac is modified from the function
# idaes.core.util.scaling.constraint_autoscale_large_jac
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import logging
import abc

import pyomo.environ as pyo
from pyomo.common.collections import Bunch
from pyomo.common.dependencies import attempt_import
from pyomo.common.modeling import unique_component_name
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.interfaces.external_grey_box import ExternalGreyBoxBlock

from idaes.core.util.scaling import (
    get_scaling_factor,
    set_scaling_factor,
    unset_scaling_factor,
)
from idaes.logger import getLogger

IPython, IPython_available = attempt_import("IPython")

_log = getLogger("watertap.core")

_pyomo_nl_writer_log = logging.getLogger("pyomo.repn.plugins.nl_writer")


def _pyomo_nl_writer_logger_filter(record):
    msg = record.getMessage()
    if "scaling_factor" in msg and "model contains export suffix" in msg:
        return False
    return True


class _WaterTAPSolverWrapper(abc.ABC):

    name = None
    base_solver = None

    @abc.abstractmethod
    def _set_options(self, solver):
        raise NotImplementedError

    @abc.abstractmethod
    def _get_pyomo_nlp(self, model):
        raise NotImplementedError

    def __init__(self, **kwds):
        kwds["name"] = self.name
        self.options = Bunch()
        for opt_key, opt_val in kwds.get("options", {}).items():
            setattr(self.options, opt_key, opt_val)

    def __getattr__(self, attr):
        # if not available here, ask the base_solver
        try:
            return getattr(pyo.SolverFactory(self.base_solver), attr)
        except AttributeError:
            raise

    def solve(self, blk, *args, **kwds):

        solver = pyo.SolverFactory(self.base_solver)
        self._tee = kwds.get("tee", False)

        self._original_options = self.options

        self.options = Bunch()
        self.options.update(self._original_options)
        self.options.update(kwds.pop("options", {}))

        # Set the default watertap options
        if "tol" not in self.options:
            self.options["tol"] = 1e-08
        if "constr_viol_tol" not in self.options:
            self.options["constr_viol_tol"] = 1e-08
        if "acceptable_constr_viol_tol" not in self.options:
            self.options["acceptable_constr_viol_tol"] = 1e-08
        if "bound_relax_factor" not in self.options:
            self.options["bound_relax_factor"] = 0.0
        if "honor_original_bounds" not in self.options:
            self.options["honor_original_bounds"] = "no"

        if not self._is_user_scaling():
            self._set_options(solver)
            try:
                return solver.solve(blk, *args, **kwds)
            finally:
                self._options_cleanup()

        if self._tee:
            print(
                f"{self.name}: {self.base_solver} with user variable scaling and IDAES jacobian constraint scaling"
            )
        # past here we need the AmplInterface
        if not AmplInterface.available():
            raise RuntimeError("Pynumero not available.")

        _pyomo_nl_writer_log.addFilter(_pyomo_nl_writer_logger_filter)
        nlp = self._scale_constraints(blk)

        # set different default for `alpha_for_y` if this is an LP
        # see: https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_alpha_for_y
        if nlp is not None:
            if nlp.nnz_hessian_lag() == 0:
                if "alpha_for_y" not in self.options:
                    self.options["alpha_for_y"] = "bound-mult"

        # Now set the options to be used by Ipopt
        # as we've popped off the above in _get_option
        self._set_options(solver)

        try:
            return solver.solve(blk, *args, **kwds)
        finally:
            self._cleanup()

    def _scale_constraints(self, blk):
        # These options are typically available with gradient-scaling, and they
        # have corresponding options in the IDAES constraint_autoscale_large_jac
        # function. Here we use their Ipopt names and default values, see
        # https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_NLP_Scaling
        max_grad = self._get_option("nlp_scaling_max_gradient", 100)
        min_scale = self._get_option("nlp_scaling_min_value", 1e-08)

        # These options are custom for the IDAES constraint_autoscale_large_jac
        # function. We expose them as solver options as this has become part
        # of the solve process.
        ignore_variable_scaling = self._get_option("ignore_variable_scaling", False)
        ignore_constraint_scaling = self._get_option("ignore_constraint_scaling", False)

        self._cache_scaling_factors(blk)

        # NOTE: This function sets the scaling factors on the
        #       constraints. Hence we cache the constraint scaling
        #       factors and reset them to their original values
        #       so that repeated calls to solve change the scaling
        #       each time based on the initial values, just like in Ipopt.
        self._add_dummy_objective(blk)
        try:
            nlp = self._get_pyomo_nlp(blk)
            _constraint_autoscale_large_jac(
                nlp,
                ignore_constraint_scaling=ignore_constraint_scaling,
                ignore_variable_scaling=ignore_variable_scaling,
                max_grad=max_grad,
                min_scale=min_scale,
            )
        except Exception as err:
            nlp = None
            if str(err) == "Error in AMPL evaluation":
                print(
                    "ipopt-watertap: Issue in AMPL function evaluation; Jacobian constraint scaling not applied."
                )
                halt_on_ampl_error = self.options.get("halt_on_ampl_error", "yes")
                if halt_on_ampl_error == "no":
                    print(
                        "ipopt-watertap: halt_on_ampl_error=no, so continuing with optimization."
                    )
                else:
                    self._cleanup()
                    raise RuntimeError(
                        "Error in AMPL evaluation.\n"
                        "Re-run ipopt with:\n"
                        '1. solver options = {"halt_on_ampl_error" : "yes", "nlp_scaling_method" : "gradient-based"\n'
                        "2. set keyword argument symbolic_solver_labels=True in the pyomo solve function call to see the affected function.\n"
                    )
            else:
                print("Error in constraint_autoscale_large_jac")
                self._cleanup()
                raise

        return nlp

    def _options_cleanup(self):
        self.options = self._original_options
        del self._original_options

    def _cleanup(self):
        self._options_cleanup()
        self._reset_scaling_factors()
        self._remove_dummy_objective()
        _pyomo_nl_writer_log.removeFilter(_pyomo_nl_writer_logger_filter)

    def _cache_scaling_factors(self, blk):
        self._scaling_cache = [
            (c, get_scaling_factor(c))
            for c in blk.component_data_objects(
                pyo.Constraint, active=True, descend_into=True
            )
        ]

    def _reset_scaling_factors(self):
        for c, s in self._scaling_cache:
            if s is None:
                unset_scaling_factor(c)
            else:
                set_scaling_factor(c, s)
        del self._scaling_cache

    def _add_dummy_objective(self, blk):
        # Pynumero requires an objective, but I don't, so let's see if we have one
        n_obj = 0
        for c in blk.component_data_objects(pyo.Objective, active=True):
            n_obj += 1
        # Add an objective if there isn't one
        if n_obj == 0:
            self._dummy_objective = pyo.Objective(expr=0)
            name = unique_component_name(blk, "objective")
            blk.add_component(name, self._dummy_objective)
        else:
            self._dummy_objective = None

    def _remove_dummy_objective(self):
        if self._dummy_objective is not None:
            # delete dummy objective
            blk = self._dummy_objective.parent_block()
            blk.del_component(self._dummy_objective)
        del self._dummy_objective

    def _get_option(self, option_name, default_value):
        # NOTE: The options are already copies of the original,
        #       so it is safe to pop them so they don't get sent to Ipopt.
        option_value = self.options.pop(option_name, None)
        if option_value is None:
            option_value = default_value
        else:
            if self._tee:
                print(f"{self.name}: {option_name}={option_value}")
        return option_value

    def _is_user_scaling(self):
        if "nlp_scaling_method" not in self.options:
            self.options["nlp_scaling_method"] = "user-scaling"
        if self.options["nlp_scaling_method"] != "user-scaling":
            if self._tee:
                print(
                    "The ipopt-watertap solver is designed to be run with user-scaling. "
                    f"Ipopt with nlp_scaling_method={self.options['nlp_scaling_method']} will be used instead"
                )
            return False
        return True


@pyo.SolverFactory.register(
    "ipopt-watertap",
    doc="The Ipopt NLP solver, with user-based variable and automatic Jacobian constraint scaling",
)
class IpoptWaterTAP(_WaterTAPSolverWrapper):

    name = "ipopt-watertap"
    base_solver = "ipopt"

    def _set_options(self, solver):
        for k, v in self.options.items():
            solver.options[k] = v

    def _get_pyomo_nlp(self, blk):
        return PyomoNLP(blk)


@pyo.SolverFactory.register(
    "cyipopt-watertap",
    doc="The Ipopt NLP solver, with user-based variable and automatic Jacobian constraint scaling",
)
class CyIpoptWaterTAP(_WaterTAPSolverWrapper):

    name = "cyipopt-watertap"
    base_solver = "cyipopt"

    def _set_options(self, solver):
        for k, v in self.options.items():
            solver.config.options[k] = v

    def _get_pyomo_nlp(self, blk):
        greyboxes = []
        try:
            for greybox in blk.component_objects(
                ExternalGreyBoxBlock, descend_into=True
            ):
                greybox.parent_block().reclassify_component_type(greybox, pyo.Block)
                greyboxes.append(greybox)

            nlp = PyomoNLP(blk)

        finally:
            for greybox in greyboxes:
                greybox.parent_block().reclassify_component_type(
                    greybox, ExternalGreyBoxBlock
                )

        return nlp


def _constraint_autoscale_large_jac(
    nlp,
    ignore_constraint_scaling=False,
    ignore_variable_scaling=False,
    max_grad=100,
    min_scale=1e-6,
):
    """Automatically scale constraints based on the Jacobian.  This function
    imitates Ipopt's default constraint scaling.  This scales constraints down
    to avoid extremely large values in the Jacobian.  This function also returns
    the unscaled and scaled Jacobian matrixes and the Pynumero NLP which can be
    used to identify the constraints and variables corresponding to the rows and
    comlumns.

    Args:
        nlp: model to scale
        ignore_constraint_scaling: ignore existing constraint scaling
        ignore_variable_scaling: ignore existing variable scaling
        max_grad: maximum value in Jacobian after scaling, subject to minimum
            scaling factor restriction.
        min_scale: minimum scaling factor allowed, keeps constraints from being
            scaled too much.

    Returns:
        None
    """
    jac = nlp.evaluate_jacobian().tocsr()
    # Get lists of variables and constraints to translate Jacobian indexes
    # save them on the NLP for later, since generating them seems to take a while
    nlp.clist = clist = nlp.get_pyomo_constraints()
    nlp.vlist = vlist = nlp.get_pyomo_variables()
    # Create a scaled Jacobian to account for variable scaling, for now ignore
    # constraint scaling
    jac_scaled = jac.copy()
    for i, c in enumerate(clist):
        for j in jac_scaled[i].indices:
            v = vlist[j]
            if ignore_variable_scaling:
                sv = 1
            else:
                sv = get_scaling_factor(v, default=1)
            jac_scaled[i, j] = jac_scaled[i, j] / sv
    # calculate constraint scale factors
    for i, c in enumerate(clist):
        sc = get_scaling_factor(c, default=1)
        if ignore_constraint_scaling or get_scaling_factor(c) is None:
            sc = 1
            row = jac_scaled[i]
            for d in row.indices:
                row[0, d] = abs(row[0, d])
            mg = row.max()
            if mg > max_grad:
                sc = max(min_scale, max_grad / mg)
            set_scaling_factor(c, sc)
        for j in jac_scaled[i].indices:
            # update the scaled jacobian
            jac_scaled[i, j] = jac_scaled[i, j] * sc
    return None


class _BaseDebugSolverWrapper:

    # defined by the derived class,
    # created on the fly
    _base_solver = None
    _debug_solver_name = None

    def __init__(self, **kwds):

        kwds["name"] = self._debug_solver_name
        self.options = Bunch()
        if kwds.get("options") is not None:
            for key in kwds["options"]:
                setattr(self.options, key, kwds["options"][key])

        self._value_cache = pyo.ComponentMap()

    def restore_initial_values(self, blk):
        for var in blk.component_data_objects(pyo.Var, descend_into=True):
            var.set_value(self._value_cache[var], skip_validation=True)

    def _cache_initial_values(self, blk):
        for v in blk.component_data_objects(pyo.Var, descend_into=True):
            self._value_cache[v] = v.value

    def solve(self, blk, *args, **kwds):

        if not IPython_available:
            raise ImportError(f"The DebugSolverWrapper requires ipython.")

        solver = pyo.SolverFactory(self._base_solver)

        for k, v in self.options.items():
            solver.options[k] = v

        self._cache_initial_values(blk)

        try:
            results = solver.solve(blk, *args, **kwds)
        except:
            results = None
        if results is not None and pyo.check_optimal_termination(results):
            return results

        # prevent circular imports
        from watertap.core.util import model_debug_mode

        # deactivate the model debug mode so we don't
        # nest this environment within itself
        model_debug_mode.deactivate()

        self.restore_initial_values(blk)
        debug = self

        # else there was a problem
        print(f"Solver debugging mode: the block {blk.name} failed to solve.")
        print(f"{blk.name} is called `blk` in this context.")
        print(f"The solver {solver.name} is available in the variable `solver`.")
        print(f"The Initial values have be restored into the block.")
        print(
            f"You can restore them anytime by calling `debug.restore_initial_values(blk)`."
        )
        print(
            f"The model has been loaded into an IDAES DiagnosticsToolbox instance called `dt`."
        )
        from idaes.core.util.model_diagnostics import DiagnosticsToolbox

        dt = DiagnosticsToolbox(blk)
        # dt.report_structural_issues()
        IPython.embed(colors="neutral")

        # activate the model debug mode
        # to keep the state the same
        model_debug_mode.activate()

        return results


def create_debug_solver_wrapper(solver_name):

    assert pyo.SolverFactory(solver_name).available()

    debug_solver_name = f"debug-solver-wrapper-{solver_name}"

    @pyo.SolverFactory.register(
        debug_solver_name,
        doc=f"Debug solver wrapper for {solver_name}",
    )
    class DebugSolverWrapper(_BaseDebugSolverWrapper):
        _base_solver = solver_name
        _debug_solver_name = debug_solver_name

    return debug_solver_name
