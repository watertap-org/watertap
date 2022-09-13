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

import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1595_photothermal_membrane_candoP.amo_1595 as amo_1595


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = amo_1595.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["LCOW_with_revenue"] = m.fs.costing.LCOW_with_revenue
    outputs["LCOT"] = m.fs.costing.LCOT
    outputs["LCOT_with_revenue"] = m.fs.costing.LCOT_with_revenue

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=2, interpolate_nan_outputs=True, save_outputs=None):
    m = amo_1595.main()[0]
    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)
    sweep_params = {}

    if case_num == 1:
        sweep_params["bcp_cost"] = LinearSample(m.fs.costing.bcp_cost, 0.1, 1, nx)
    elif case_num == 2:
        sweep_params["nitrous_oxide_cost"] = LinearSample(
            m.fs.costing.nitrous_oxide_cost, 0, 5, nx
        )
    elif case_num == 3:
        sweep_params["water_cost"] = LinearSample(m.fs.costing.water_cost, 0, 2, nx)
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    # run sweep
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=save_outputs,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    print(global_results)

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(*sys.argv[1:])
