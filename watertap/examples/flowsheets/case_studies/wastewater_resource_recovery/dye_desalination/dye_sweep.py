import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination as dye_desalination


def set_up_sensitivity(m):
    """
    Initialize optimization sweep for specific cost parameters and baselines
    :param m:
    :return:
    """
    outputs = {}
    optimize_kwargs = {"check_termination": False}  # None
    opt_function = dye_desalination.solve

    # add costing parameters

    # add baselines

    # create outputs

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True):
    """
    :param case_num:
    :param nx:
    :param interpolate_nan_outputs:
    :return:
    """
    m = dye_desalination.main()

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num != 0:
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

    return global_results, sweep_params


def main(case_num=1, nx=1, interpolate_nan_outputs=True):
    global_results, sweep_params = run_analysis(case_num, nx, interpolate_nan_outputs)
    print(global_results)
    return


if __name__ == "__main__":
    main()
