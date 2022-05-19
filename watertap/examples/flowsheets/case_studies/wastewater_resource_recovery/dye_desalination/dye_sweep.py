import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination as dye_desalination


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = dye_desalination.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW_comp
    # outputs["Total_Cost"] = m.fs.costing.total_annualized_cost
    # outputs["LCODS"] = m.fs.costing.LCODS

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True):

    m = dye_desalination.main()

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)
    m.fs.costing.dye_disposal_cost.unfix()

    sweep_params = {}

    if case_num == 1:
        sweep_params["disposal_cost"] = LinearSample(
            m.fs.costing.dye_disposal_cost, 5, 10, nx
        )
    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    output_filename = "sensitivity_" + str(case_num) + ".csv"
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


def main(case_num=1, nx=11, interpolate_nan_outputs=False):
    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    comm, rank, num_procs = _init_mpi()  # TODO - resolve MPI communicator init issue

    global_results, sweep_params, m = run_analysis(
        case_num, nx, interpolate_nan_outputs
    )
    print(global_results)

    return global_results, m


if __name__ == "__main__":
    results, model = main()
