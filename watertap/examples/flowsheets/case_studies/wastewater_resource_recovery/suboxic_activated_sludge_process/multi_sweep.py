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
from watertap.tools.parameter_sweep import (
    LinearSample,
    parameter_sweep,
)
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.suboxic_activated_sludge_process.suboxic_activated_sludge_process as suboxic_activated_sludge_process


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"fail_flag": False}
    opt_function = suboxic_activated_sludge_process.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=11, interpolate_nan_outputs=True):

    m = suboxic_activated_sludge_process.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        m.fs.suboxicASM.removal_frac_mass_comp[0, "bod"].unfix()
        sweep_params["bod_removal"] = LinearSample(
            m.fs.suboxicASM.removal_frac_mass_comp[0, "bod"], 0.9363, 0.9904, nx
        )
    elif case_num == 2:
        m.fs.suboxicASM.removal_frac_mass_comp[0, "tss"].unfix()
        sweep_params["tss_removal"] = LinearSample(
            m.fs.suboxicASM.removal_frac_mass_comp[0, "tss"], 0.9405, 0.9940, nx
        )
    elif case_num == 3:
        m.fs.suboxicASM.removal_frac_mass_comp[0, "tkn"].unfix()
        sweep_params["tkn_removal"] = LinearSample(
            m.fs.suboxicASM.removal_frac_mass_comp[0, "tkn"], 0.8846, 0.9615, nx
        )
    elif case_num == 4:
        m.fs.costing.electricity_cost.unfix()
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.05, 0.12, nx
        )
    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    output_filename = "sensitivity_" + str(case_num) + ".csv"

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
        interpolate_nan_outputs=interpolate_nan_outputs,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    results, sweep_params, m = run_analysis()
