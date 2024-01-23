"""
Minimal Intractable System (MIS) finder

See: https://www.sce.carleton.ca/faculty/chinneck/docs/CPAIOR07InfeasibilityTutorial.pdf

"""
import pyomo.environ as pyo

from pyomo.core.plugins.transform.add_slack_vars import AddSlackVariables

from pyomo.core.plugins.transform.hierarchy import IsomorphicTransformation

from pyomo.common.modeling import unique_component_name
from pyomo.common.collections import ComponentMap, ComponentSet

from idaes.core.solvers import get_solver


class VariableBoundsAsConstraints(IsomorphicTransformation):
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


def compute_MIS(model, solver=None, tee=True, epsilon=1e-8):

    # hold the original harmless
    m = model.clone()

    if solver is None:
        solver = get_solver()
    elif isinstance(solver, str):
        solver = pyo.SolverFactory(solver)
    else:
        # assume we have a solver
        assert solver.available()

    # first, cache the values we get
    _value_cache = ComponentMap()
    for v in m.component_data_objects(pyo.Var, descend_into=True):
        _value_cache[v] = v.value

    # move the variable bounds to the constraints
    VariableBoundsAsConstraints().apply_to(m)
    # now add slacks to every constraint (including variable bounds)
    AddSlackVariables().apply_to(m)
    slack_block = m._core_add_slack_variables

    # collect the initial set of infeasibilities
    results = solver.solve(m, tee=tee)

    # this first solve should be feasible, by definition
    pyo.assert_optimal_termination(results)

    # TODO: For an interior point method, we should not
    #       need this loop -- anything potentially
    #       contributing to the infeasibility will
    #       be caught above
    # Phase 1 -- build the initial set of constraints
    elastic_filter = ComponentSet()
    while pyo.check_optimal_termination(results):
        fixed_var = False
        for v in slack_block.component_data_objects(pyo.Var):
            if v.value > epsilon:
                v.fix(0)
                fixed_var = True
                if "_slack_plus_" in v.name:
                    elastic_filter.add(m.find_component(v.local_name[len("_slack_plus_"):]))
                elif "_slack_minus_" in v.name:
                    elastic_filter.add(m.find_component(v.local_name[len("_slack_minus_"):]))
                else:
                    raise RuntimeError("Bad var name {v.name}")
        if not fixed_var:
            raise RuntimeError("Could not detect infeasibility")
        for var, val in _value_cache.items():
            var.value = val
        results = solver.solve(m, tee=tee)

    # Phase 2 -- deletion filter

    # remove slacks by fixing them to 0
    for v in slack_block.component_data_objects(pyo.Var):
        v.fix(0)

    # mark all constraints not in the filter as inactive
    for c in m.component_data_objects(pyo.Constraint):
        if c in elastic_filter:
            continue
        else:
            c.deactivate()

    deletion_filter = []
    guards = []
    for constr in elastic_filter:
        constr.deactivate()
        for var, val in _value_cache.items():
            var.value = val

        try:
            results = solver.solve(m, tee=tee)
        except:
            math_failure = True
        else:
            math_failure = False

        if math_failure:
            constr.activate()
            gaurds.append(constr)
        elif pyo.check_optimal_termination(results):
            constr.activate()
            deletion_filter.append(constr)
        else: # still infeasible without this constraint
            pass

    print("Constraints / bounds in set:")
    _print_results(deletion_filter)
    print("Constraints / bounds in guards:")
    _print_results(guards)
    return deletion_filter, guards

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
