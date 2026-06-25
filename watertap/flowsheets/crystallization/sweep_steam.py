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
from parameter_sweep import LinearSample, parameter_sweep
from watertap.core.solvers import get_solver
from watertap.flowsheets.crystallization.crystallizer_TVC import (
    build,
    set_operating_conditions,
    initialize_system,
    optimize_set_up,
    solve,
)


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {}
    opt_function = solve
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["Steam Cost"] = m.fs.costing.steam.cost
    return outputs, optimize_kwargs, opt_function


def run_steam_cost_sweep(nx=200, output_filename="steam_cost_sweep.csv"):
    solver = get_solver()
    m = build()
    set_operating_conditions(m)
    initialize_system(m, solver=solver)
    optimize_set_up(m)

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}

    sweep_params["steam_cost"] = LinearSample(m.fs.costing.steam.cost, 0.00, 0.008, nx)

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_steam_cost_sweep()
