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

__all__ = [
    "get_solver",
]

# from watertap_solvers import get_solver
from idaes.core.solvers import get_solver as idaes_get_solver


def get_solver(
    solver="ipopt_v2",
    solver_options: dict = None,
    writer_config: dict = None,
    options: dict = None,
):
    if options is not None:
        if solver_options is not None:
            raise RuntimeError()
        else:
            solver_options = options

    if solver_options is None:
        solver_options = {}

    if writer_config is None:
        writer_config = {}

    solver = "ipopt_v2"

    solver_options["nlp_scaling_method"] = "gradient-based"

    if "tol" not in solver_options:
        solver_options["tol"] = 1e-08
    if "constr_viol_tol" not in solver_options:
        solver_options["constr_viol_tol"] = 1e-08
    if "acceptable_constr_viol_tol" not in solver_options:
        solver_options["acceptable_constr_viol_tol"] = 1e-08
    if "bound_relax_factor" not in solver_options:
        solver_options["bound_relax_factor"] = 0.0
    if "honor_original_bounds" not in solver_options:
        solver_options["honor_original_bounds"] = "no"
    if "linear_solver" not in solver_options:
        solver_options["linear_solver"] = "ma27"
    if "max_iter" not in solver_options:
        solver_options["max_iter"] = 3000

    if "scale_model" not in writer_config:
        writer_config["scale_model"] = True
    if "linear_presolve" not in writer_config:
        writer_config["linear_presolve"] = False
    if "skip_trivial_constraints" not in writer_config:
        writer_config["skip_trivial_constraints"] = False

    return idaes_get_solver(
        solver=solver, solver_options=solver_options, writer_config=writer_config
    )
