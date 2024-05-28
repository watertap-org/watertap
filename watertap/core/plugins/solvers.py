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

import logging

import pyomo.environ as pyo
from pyomo.common.collections import Bunch
from pyomo.solvers.plugins.solvers.IPOPT import IPOPT

import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    get_scaling_factor,
    set_scaling_factor,
    unset_scaling_factor,
)
from idaes.logger import getLogger

_log = getLogger("watertap.core")

_pyomo_nl_writer_log = logging.getLogger("pyomo.repn.plugins.nl_writer")


def _pyomo_nl_writer_logger_filter(record):
    msg = record.getMessage()
    if "scaling_factor" in msg and "model contains export suffix" in msg:
        return False
    return True


@pyo.SolverFactory.register(
    "ipopt-watertap",
    doc="The Ipopt NLP solver, with user-based variable and automatic Jacobian constraint scaling",
)
class IpoptWaterTAP:

    name = "ipopt-watertap"
    _base_solver = IPOPT

    def __init__(self, **kwds):
        kwds["name"] = self.name
        self.options = Bunch()
        for opt_key, opt_val in kwds.get("options", {}).items():
            setattr(self.options, opt_key, opt_val)

    def __getattr__(self, attr):
        # if not available here, ask the _base_solver
        try:
            return getattr(self._base_solver(), attr)
        except AttributeError:
            raise

    def solve(self, blk, *args, **kwds):

        solver = self._base_solver()
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
            for k, v in self.options.items():
                solver.options[k] = v
            try:
                return solver.solve(blk, *args, **kwds)
            finally:
                self._options_cleanup()

        if self._tee:
            print(
                "ipopt-watertap: Ipopt with user variable scaling and IDAES jacobian constraint scaling"
            )

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
        for k, v in self.options.items():
            solver.options[k] = v

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
        try:
            _, _, nlp = iscale.constraint_autoscale_large_jac(
                blk,
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
                        "1. solver options = {\"halt_on_ampl_error\" : \"yes\", \"nlp_scaling_method\" : \"gradient-based\"\n"
                        "2. set keyword argument symbolic_solver_labels=True in the pyomo solve function call to see the affected function.\n"                    )
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

    def _get_option(self, option_name, default_value):
        # NOTE: The options are already copies of the original,
        #       so it is safe to pop them so they don't get sent to Ipopt.
        option_value = self.options.pop(option_name, None)
        if option_value is None:
            option_value = default_value
        else:
            if self._tee:
                print(f"ipopt-watertap: {option_name}={option_value}")
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
