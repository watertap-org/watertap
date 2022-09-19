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

import os, sys

from watertap.tools.parameter_sweep import LinearSample, parameter_sweep

import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.peracetic_acid_disinfection.peracetic_acid_disinfection as peracetic_acid_disinfection


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = peracetic_acid_disinfection.solve

    # create outputs
    outputs["LCOT"] = m.fs.costing.LCOT

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True, save_path=None):
    m = peracetic_acid_disinfection.main()[0]
    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)
    sweep_params = {}

    if case_num == 1:
        sweep_params["disinfection_solution_cost"] = LinearSample(
            m.fs.costing.disinfection_solution_cost, 2, 10, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=save_path,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    print(global_results)

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(*sys.argv[1:])
