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

import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.groundwater_treatment.groundwater_treatment as groundwater_treatment


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = groundwater_treatment.solve

    # create outputs
    outputs["LCOT"] = m.fs.costing.LCOT

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True):
    m = groundwater_treatment.main()[0]
    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)
    sweep_params = {}

    if case_num == 1:
        sweep_params["filtration_media_disposal_cost"] = LinearSample(
            m.fs.costing.filtration_media_disposal_cost, 10, 1000, nx
        )
    elif case_num == 2:
        sweep_params["filtration_media_cost"] = LinearSample(
            m.fs.costing.filtration_media_cost, 100, 5000, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    output_filename = "sensitivity_" + str(case_num) + ".csv"
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        output_filename,
    )

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_path,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    print(global_results)

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(*sys.argv[1:])
