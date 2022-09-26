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
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.biomembrane_filtration.biomembrane_filtration as biomembrane_filtration


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = biomembrane_filtration.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True, save_outputs=None):

    m = biomembrane_filtration.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        m.fs.mabr.removal_frac_mass_comp[0, "ammonium_as_nitrogen"].unfix()
        sweep_params["ammonium_removal"] = LinearSample(
            m.fs.mabr.removal_frac_mass_comp[0, "ammonium_as_nitrogen"], 0.5, 0.9, nx
        )
    elif case_num == 2:
        m.fs.dmbr.removal_frac_mass_comp[0, "nitrate"].unfix()
        sweep_params["nitrate_removal"] = LinearSample(
            m.fs.dmbr.removal_frac_mass_comp[0, "nitrate"], 0.7, 0.95, nx
        )
    elif case_num == 3:
        m.fs.dmbr.removal_frac_mass_comp[0, "tss"].unfix()
        sweep_params["tss_removal"] = LinearSample(
            m.fs.dmbr.removal_frac_mass_comp[0, "tss"], 0.6, 0.95, nx
        )
    elif case_num == 4:
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.05, 0.12, nx
        )
    else:
        raise ValueError("case_num = %d not recognized." % case_num)

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=save_outputs,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis(*sys.argv[1:])
