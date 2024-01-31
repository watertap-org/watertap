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

import logging

import pyomo.environ as pyo
from pyomo.common.collections import Bunch
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.common.modeling import unique_component_name
from pyomo.core.base.block import _BlockData
from pyomo.core.kernel.block import IBlock
from pyomo.solvers.plugins.solvers.IPOPT import IPOPT

from pyomo.core.plugins.transform.add_slack_vars import AddSlackVariables
from pyomo.core.plugins.transform.hierarchy import IsomorphicTransformation

from pyomo.contrib.incidence_analysis import get_incident_variables, IncidenceMethod

from pyomo.opt import WriterFactory

import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    get_scaling_factor,
    set_scaling_factor,
    unset_scaling_factor,
)
from idaes.logger import getLogger

_log = getLogger("watertap.core")
_default_nl_writer = WriterFactory.get_class("nl")

_pyomo_nl_writer_log = logging.getLogger("pyomo.repn.plugins.nl_writer")


def _pyomo_nl_writer_logger_filter(record):
    msg = record.getMessage()
    if "scaling_factor" in msg and "model contains export suffix" in msg:
        return False
    return True


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

        # until proven otherwise
        self._cleanup_needed = False

        self._tee = kwds.get("tee", False)

        # Set the default watertap options
        if "tol" not in self.options:
            self.options["tol"] = 1e-08
        if "constr_viol_tol" not in self.options:
            self.options["constr_viol_tol"] = 1e-08
        if "bound_relax_factor" not in self.options:
            self.options["bound_relax_factor"] = 0.0
        if "honor_original_bounds" not in self.options:
            self.options["honor_original_bounds"] = "no"

        if not self._is_user_scaling():
            super()._presolve(*args, **kwds)
            self._cleanup()
            return

        if self._tee:
            print(
                "ipopt-watertap: Ipopt with user variable scaling and IDAES jacobian constraint scaling"
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
        self._cleanup_needed = True
        _pyomo_nl_writer_log.addFilter(_pyomo_nl_writer_logger_filter)

        # NOTE: This function sets the scaling factors on the
        #       constraints. Hence we cache the constraint scaling
        #       factors and reset them to their original values
        #       so that repeated calls to solve change the scaling
        #       each time based on the initial values, just like in Ipopt.
        try:
            _, _, nlp = iscale.constraint_autoscale_large_jac(
                self._model,
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
                        "Run ipopt with halt_on_ampl_error=yes and symbolic_solver_labels=True to see the affected function."
                    )
            else:
                print("Error in constraint_autoscale_large_jac")
                self._cleanup()
                raise

        # set different default for `alpha_for_y` if this is an LP
        # see: https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_alpha_for_y
        if nlp is not None:
            if nlp.nnz_hessian_lag() == 0:
                if "alpha_for_y" not in self.options:
                    self.options["alpha_for_y"] = "bound-mult"

        try:
            # this creates the NL file, among other things
            return super()._presolve(*args, **kwds)
        except:
            self._cleanup()
            raise

    def _cleanup(self):
        if self._cleanup_needed:
            self._reset_scaling_factors()
            # remove our reference to the model
            del self._model
            _pyomo_nl_writer_log.removeFilter(_pyomo_nl_writer_logger_filter)

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


class _VariableBoundsAsConstraints(IsomorphicTransformation):
    """Replace all variables bounds and domain information with constraints.

       Leaves fixed Vars untouched (for now)
    """

    def _apply_to(self, instance, **kwds):

        boundconstrblockname = unique_component_name(instance, "_variable_bounds")
        instance.add_component(boundconstrblockname, pyo.Block())
        boundconstrblock = instance.component(boundconstrblockname)

        for v in instance.component_data_objects(pyo.Var, descend_into=True):
            if v.fixed:
                continue
            lb, ub = v.bounds
            if lb is None and ub is None:
                continue
            var_name = v.getname(fully_qualified=True)
            if lb is not None:
                con_name = "lb_for_"+var_name
                con = pyo.Constraint(expr=(lb, v, None))
                boundconstrblock.add_component(con_name, con)
            if ub is not None:
                con_name = "ub_for_"+var_name
                con = pyo.Constraint(expr=(None, v, ub))
                boundconstrblock.add_component(con_name, con)

            # now we deactivate the variable bounds / domain
            v.domain = pyo.Reals
            v.setlb(None)
            v.setub(None)


@pyo.SolverFactory.register(
    "ipopt-watertap",
    doc="The Ipopt NLP solver, with user-based variable and automatic Jacobian constraint scaling",
)
class WaterTAPSolver:

    def __init__(self, **kwds):

        self.options = Bunch()
        if kwds.get("options") is not None:
            for key in kwds["options"]:
                setattr(self.options, key, kwds["options"][key])

    def solve(self, model, *args, **kwds):

        base_solver = IpoptWaterTAP()

        if "constr_viol_tol" not in self.options:
            self.options["constr_viol_tol"] = 1e-08

        for k,v in self.options.items():
            base_solver.options[k] = v

        # first, cache the values we get
        _value_cache = ComponentMap()
        for v in model.component_data_objects(pyo.Var, descend_into=True):
            _value_cache[v] = v.value

        try:
            results = base_solver.solve(model, *args, **kwds)
            if pyo.check_optimal_termination(results):
                return results
        except:
            pass


        # Something went wrong. We'll try to find a feasible point
        # *or* prove that the model is infeasible
        if model.parent_block() is None:
            common_name = ""
        else:
            common_name = model.name

        modified_model = model.clone()

        _modified_model_var_to_original_model_var = ComponentMap()
        _modified_model_value_cache = ComponentMap()

        for v in model.component_data_objects(pyo.Var, descend_into=True):
            modified_model_var = modified_model.find_component(v.name[len(common_name)+1:])

            _modified_model_var_to_original_model_var[modified_model_var] = v
            _modified_model_value_cache[modified_model_var] = _value_cache[v]
            modified_model_var.set_value(_value_cache[v], skip_validation=True)

        # TODO: For WT / IDAES models, we should probably be more
        #       selective in *what* we elasticize. E.g., it probably
        #       does not make sense to elasticize property calculations
        #       and maybe certain other equality constraints calculating
        #       values. Maybe we shouldn't elasticize *any* equality 
        #       constraints.
        #       For example, elasticizing the calculation of mass fraction
        #       makes absolutely no sense and will just be noise for the
        #       modeler to sift through. We could try to sort the constraints
        #       such that we look for those with linear coefficients `1` on
        #       some term and leave those be.
        # move the variable bounds to the constraints
        _VariableBoundsAsConstraints().apply_to(modified_model)

        _var_to_constraints_map = ComponentMap()
        _constraint_to_vars_map = ComponentMap()
        for c in modified_model.component_data_objects(pyo.Constraint, descend_into=True):
            _constraint_to_vars_map[c] = []
            for v in get_incident_variables(c.body, method=IncidenceMethod.standard_repn, linear_only=False, include_fixed=False):
                _constraint_to_vars_map[c].append(v)
                if v not in _var_to_constraints_map:
                    _var_to_constraints_map[v] = []
                _var_to_constraints_map[v].append(c)

        AddSlackVariables().apply_to(modified_model)
        slack_block = modified_model._core_add_slack_variables
        # don't care about dual feasibility so much
        #slack_block._slack_objective.expr /= 1e3

        for v in slack_block.component_data_objects(pyo.Var):
            v.fix(0)
        # start with variable bounds -- these are the easist to interpret
        for c in modified_model._variable_bounds.component_data_objects(pyo.Constraint, descend_into=True):
            plus = slack_block.component(f"_slack_plus_{c.name}")
            minus = slack_block.component(f"_slack_minus_{c.name}")
            assert not (plus is None and minus is None)
            if plus is not None:
                plus.unfix()
            if minus is not None:
                minus.unfix()

        # collect the initial set of infeasibilities
        # we keep the constraint_viol_tol tight
        results = base_solver.solve(modified_model, *args, **kwds)

        if pyo.check_optimal_termination(results):
            msg = f"Model {model.name} may be infeasible. A feasible solution was found with the following variable bounds relaxed:"
            for v in slack_block.component_data_objects(pyo.Var):
                if v.value > self.options["constr_viol_tol"]:
                    c = _get_constraint(modified_model, v)
                    assert "_variable_bounds" in c.name
                    name = c.local_name
                    if "lb" in name:
                        msg += f"\n\tlb of var {name[7:]} by {v.value}"
                    elif "ub" in name:
                        msg += f"\n\tub of var {name[7:]} by {v.value}"
                    else:
                        raise RuntimeError("unrecongized var name")
            raise Exception(msg)

        for v in slack_block.component_data_objects(pyo.Var):
            v.unfix()
        results = base_solver.solve(modified_model, *args, **kwds)

        # this first solve should be feasible, by definition
        pyo.assert_optimal_termination(results)

        # TODO: Elasticizing too much at once seems to cause Ipopt trouble.
        #       After an initial sweep, we should just fix one elastic variable
        #       and put everything else on a stack of "constraints to elasticize".
        #       We elastisize one constraint at a time and fix one constraint at a time.
        #       After fixing an elastic variable, we elasticize a single constraint it
        #       appears in and put the remaining constraints on the stack. If the resulting problem
        #       is feasible, we keep going "down the tree". If the resulting problem is
        #       infeasible or cannot be solved, we elasticize a single constraint from
        #       the top of the stack.
        #       The algorithm stops when the stack is empty and the subproblem is infeasible.
        #       Along the way, any time the current problem is infeasible we can check to
        #       see if the current set of constraints in the filter is as a collection of
        #       infeasible constraints -- to terminate early.
        # Phase 1 -- build the initial set of constraints, or prove feasibility
        stack = []
        elastic_filter = ComponentSet()
        proved_feasibility = False
        while pyo.check_optimal_termination(results) or stack:
            if pyo.check_optimal_termination(results):
                max_viol = 0.0
                max_pair = None
                for v in slack_block.component_data_objects(pyo.Var):
                    if v.value > self.options["constr_viol_tol"]:
                        constr = _get_constraint(modified_model, v)
                        stack.append((v, constr))
                        if v.value > max_viol:
                            max_viol = v.value
                            max_pair = (v, constr)
            else: # infeasible, pop
                max_pair = stack.pop()
            if max_pair is None:
                proved_feasibility = True
                break
            else:
                max_pair[0].fix(0)
                elastic_filter.add(max_pair[1])
                print(f"Adding constraint {max_pair[1].name} to the filter")
            # for var, val in _modified_model_value_cache.items():
            #     var.set_value(val, skip_validation=True)
            # fix everything but the max, then only unfix those
            # in the graph related to the enforced constraint
            for v in slack_block.component_data_objects(pyo.Var):
                if v.value > self.options["constr_viol_tol"]:
                    continue
                else:
                    for vs, _ in stack:
                        if vs is v:
                            break
                    else: # no break
                        v.fix(0)
            for c in [max_pair[1]]:
                for v in _constraint_to_vars_map[c]:
                    for c_prime in _var_to_constraints_map[v]:
                        if c_prime not in elastic_filter:
                            plus = slack_block.component(f"_slack_plus_{c_prime.name}")
                            minus = slack_block.component(f"_slack_minus_{c_prime.name}")
                            assert not (plus is None and minus is None)
                            if plus is not None:
                                plus.unfix()
                            if minus is not None:
                                minus.unfix()
                            print(f"relaxing constraint {c_prime.name}")
            try:
                results = base_solver.solve(modified_model, *args, **kwds)
            except:
                break

        if proved_feasibility:
            # load the feasible solution into the original model
            for modified_model_var, v in _modified_model_var_to_original_model_var.items():
                v.set_value(modified_model_var.value, skip_validation=True)
            base_solver.options["bound_push"] = 0.0
            results = base_solver.solve(model, *args, **kwds)
            return results

        # Phase 2 -- deletion filter
        print(f"PHASE 2, {len(elastic_filter)} Constraints to consider.")
        # temporarily switch to nl_v1 writer
        WriterFactory.register("nl")(WriterFactory.get_class("nl_v1"))

        # remove slacks by fixing them to 0
        for v in slack_block.component_data_objects(pyo.Var):
            v.fix(0)
        for o in modified_model.component_data_objects(pyo.Objective, descend_into=True):
            o.deactivate()

        # mark all constraints not in the filter as inactive
        for c in modified_model.component_data_objects(pyo.Constraint):
            if c in elastic_filter:
                continue
            else:
                c.deactivate()

        # TODO: From here on, change Ipopt's options
        #       so it will terminate on a feasible solution
        #       (e.g., do not consider dual or complimentary
        #              or overall `tol`)
        base_solver.options["tol"] = 1e-08
        base_solver.options["dual_inf_tol"] = 1e20
        base_solver.options["compl_inf_tol"] = 1e20
        try:
            results = base_solver.solve(modified_model, *args, **kwds)
        except:
            results = None
        assert results is None or not pyo.check_optimal_termination(results)

        deletion_filter = []
        guards = []
        for constr in elastic_filter:
            constr.deactivate()
            for var, val in _modified_model_value_cache.items():
                var.set_value(val, skip_validation=True)

            try:
                results = base_solver.solve(modified_model, *args, **kwds)
            except:
                math_failure = True
            else:
                math_failure = False

            if math_failure:
                constr.activate()
                guards.append(constr)
            elif pyo.check_optimal_termination(results):
                constr.activate()
                deletion_filter.append(constr)
            else: # still infeasible without this constraint
                pass

        print("Constraints / bounds in set:")
        _print_results(deletion_filter)
        print("Constraints / bounds in guards:")
        _print_results(guards)
        WriterFactory.register("nl")(_default_nl_writer)

        raise Exception("Found model to be infeasible, see conflicting constraints above")


def _print_results(constr_list):
    for c in constr_list:
        c_name = c.name
        if "_variable_bounds" in c_name:
            name = c.local_name
            if "lb" in name:
                print(f"\tlb of var {name[7:]}")
            elif "ub" in name:
                print(f"\tub of var {name[7:]}")
            else:
                raise RuntimeError("unrecongized var name")
        else:
            print(f"\tconstraint: {c_name}")

def _get_constraint(modified_model, v):
    if "_slack_plus_" in v.name:
        constr = modified_model.find_component(v.local_name[len("_slack_plus_"):])
        if constr is None:
            raise RuntimeError("Bad constraint name {v.local_name[len('_slack_plus_'):]}")
        return constr
    elif "_slack_minus_" in v.name:
        constr = modified_model.find_component(v.local_name[len("_slack_minus_"):])
        if constr is None:
            raise RuntimeError("Bad constraint name {v.local_name[len('_slack_minus_'):]}")
        return constr
    else:
        raise RuntimeError("Bad var name {v.name}")


## reconfigure IDAES to use the ipopt-watertap solver
import idaes

_default_solver_config_value = idaes.cfg.get("default_solver")
_idaes_default_solver = _default_solver_config_value._default

_default_solver_config_value.set_default_value("ipopt-watertap")
if not _default_solver_config_value._userSet:
    _default_solver_config_value.reset()
