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

import pyomo.environ as pyo
from pyomo.core.base.block import _BlockData
from pyomo.core.kernel.block import IBlock
from pyomo.solvers.plugins.solvers.IPOPT import IPOPT

import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    get_scaling_factor,
    set_scaling_factor,
    unset_scaling_factor,
)
from idaes.logger import getLogger

_log = getLogger("watertap.core")


@pyo.SolverFactory.register(
    "ipopt-watertap",
    doc="The Ipopt NLP solver, with user-based variable and automatic Jacobian constraint scaling",
)
class IpoptWaterTAP(IPOPT):
    def __init__(self, **kwds):
        kwds["name"] = "ipopt-watertap"
        self._cleanup_needed = False
        super().__init__(**kwds)

    def _presolve(self, *args, **kwds):
        if len(args) > 1 or len(args) == 0:
            raise TypeError(
                f"IpoptWaterTAP.solve takes 1 positional argument but {len(args)} were given"
            )
        if not isinstance(args[0], (_BlockData, IBlock)):
            raise TypeError(
                "IpoptWaterTAP.solve takes 1 positional argument: a Pyomo ConcreteModel or Block"
            )

        self._tee = kwds.get("tee", False)

        # Set the default watertap options
        if "tol" not in self.options:
            self.options["tol"] = 1e-08
        if "constr_viol_tol" not in self.options:
            self.options["constr_viol_tol"] = 1e-08

        if not self._is_user_scaling():
            self._cleanup_needed = False
            return super()._presolve(*args, **kwds)

        if self._tee:
            print(
                "ipopt-watertap: Ipopt with user variable scaling and IDAES jacobian constraint scaling"
            )

        bound_relax_factor = self._get_option("bound_relax_factor", 1e-10)
        if bound_relax_factor < 0.0:
            raise ValueError(
                f"Option bound_relax_factor must be non-negative; bound_relax_factor={bound_relax_factor}"
            )

        # we are doing this ourselves, don't want Ipopt to also do it
        # also effectively turns "honor_original_bounds" off
        self.options["bound_relax_factor"] = 0.0

        # raise an error if "honor_original_bounds" is set to "yes" (for now)
        if self.options.get("honor_original_bounds", "no") == "yes":
            raise ValueError(
                f"""Option honor_original_bounds must be set to "no" -- ipopt-watertap does not presently implement this option"""
            )

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

        self._model = args[0]
        self._cache_scaling_factors()
        self._cache_and_set_relaxed_bounds(bound_relax_factor)
        self._cleanup_needed = True

        # NOTE: This function sets the scaling factors on the
        #       constraints. Hence we cache the constraint scaling
        #       factors and reset them to their original values
        #       so that repeated calls to solve change the scaling
        #       each time based on the initial values, just like in Ipopt.
        try:
            iscale.constraint_autoscale_large_jac(
                self._model,
                ignore_constraint_scaling=ignore_constraint_scaling,
                ignore_variable_scaling=ignore_variable_scaling,
                max_grad=max_grad,
                min_scale=min_scale,
            )
        except Exception as err:
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
                        "Run ipopt with halt_on_ampl_error=yes and symbolic_solver_labels=True to see the affected function."
                    )
            else:
                print("Error in constraint_autoscale_large_jac")
                self._cleanup()
                raise

        try:
            # this creates the NL file, among other things
            return super()._presolve(*args, **kwds)
        except:
            self._cleanup()
            raise

    def _cleanup(self):
        if self._cleanup_needed:
            self._reset_scaling_factors()
            self._reset_bounds()
            # remove our reference to the model
            del self._model

    def _postsolve(self):
        self._cleanup()
        return super()._postsolve()

    def _cache_scaling_factors(self):
        self._scaling_cache = [
            (c, get_scaling_factor(c))
            for c in self._model.component_data_objects(
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

    def _cache_and_set_relaxed_bounds(self, bound_relax_factor):
        self._bound_cache = pyo.ComponentMap()
        val = pyo.value
        for v in self._model.component_data_objects(
            pyo.Var, active=True, descend_into=True
        ):
            # we could hit a variable more
            # than once because of References
            if v in self._bound_cache:
                continue
            if v.lb is None and v.ub is None:
                continue
            self._bound_cache[v] = (v.lb, v.ub)
            sf = get_scaling_factor(v, default=1)
            if v.lb is not None:
                v.lb = val(
                    (v.lb * sf - bound_relax_factor * max(1, abs(val(v.lb * sf)))) / sf
                )
            if v.ub is not None:
                v.ub = val(
                    (v.ub * sf + bound_relax_factor * max(1, abs(val(v.ub * sf)))) / sf
                )

    def _reset_bounds(self):
        for v, (lb, ub) in self._bound_cache.items():
            v.lb = lb
            v.ub = ub
        del self._bound_cache

    def _get_option(self, option_name, default_value):
        # NOTE: options get reset to their original value at the end of the
        #       OptSolver.solve. The options in _presolve (where this is called)
        #       are already copies of the original, so it is safe to pop them so
        #       they don't get sent to Ipopt.
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


## reconfigure IDAES to use the ipopt-watertap solver
import idaes

_default_solver_config_value = idaes.cfg.get("default_solver")
_idaes_default_solver = _default_solver_config_value._default

_default_solver_config_value.set_default_value("ipopt-watertap")
if not _default_solver_config_value._userSet:
    _default_solver_config_value.reset()
