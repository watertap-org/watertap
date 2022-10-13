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
import os
import sys
from watertap.tools.parameter_sweep import (
    LinearSample,
    parameter_sweep,
)
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.electrochemical_nutrient_removal.electrochemical_nutrient_removal as electrochemical_nutrient_removal


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = electrochemical_nutrient_removal.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["LCOS"] = m.fs.costing.LCOS

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True, results_path=None):

    m = electrochemical_nutrient_removal.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # sensitivity analysis
        sweep_params["MgCl2_cost"] = LinearSample(
            m.fs.costing.magnesium_chloride_cost, 135, 538, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=results_path,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(*sys.argv[1:])
