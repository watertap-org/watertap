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

    return


def main(case_num=1, nx=1, interpolate_nan_outputs=True):
    run_analysis(case_num, nx, interpolate_nan_outputs)
    return


if __name__ == "__main__":
    main()
