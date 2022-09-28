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
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.swine_wwt.swine_wwt as swine_wwt


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": False}
    opt_function = swine_wwt.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.levelized_costs.LCOW
    outputs["LCOT"] = m.fs.costing.levelized_costs.LCOT
    outputs["LCOH2"] = m.fs.costing.levelized_costs.LCOH2
    outputs["LCOCOD"] = m.fs.costing.levelized_costs.LCOCOD
    outputs["LCOP"] = m.fs.costing.levelized_costs.LCOP
    outputs["LCOVFA"] = m.fs.costing.levelized_costs.LCOVFA
    outputs["LCON"] = m.fs.costing.levelized_costs.LCON

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num=1, nx=5, interpolate_nan_outputs=True, save_outputs=None):
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    m = swine_wwt.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # AnMBR-MEC Unit Capex and Opex
        sweep_params["AnMBR-MEC Unit Capex"] = LinearSample(
            m.fs.costing.anaerobic_mbr_mec.unit_capex, 0.08, 8, nx
        )
        sweep_params["AnMBR-MEC Unit Opex"] = LinearSample(
            m.fs.costing.anaerobic_mbr_mec.unit_capex, 0.006, 6, nx
        )
    elif case_num == 2:
        # AnMBR-MEC Water Recovery
        sweep_params["AnMBR-MEC H2O Recovery"] = LinearSample(
            m.fs.mbr_mec.recovery_frac_mass_H2O, 0.1, 1, nx
        )
        # AnMBR-MEC Conversion to ffCOD
        sweep_params["AnMBR-MEC Conversion to ffCOD"] = LinearSample(
            m.fs.mbr_mec.reaction_conversion[0, "cod_to_nonbiodegradable_cod"],
            0.3,
            0.7,
            nx,
        )
    elif case_num == 3:
        # VFA Recovery Unit Capex and heating costs
        sweep_params["VFA Rec Unit Capex"] = LinearSample(
            m.fs.costing.vfa_recovery.unit_capex, 0.3, 30, nx
        )
        sweep_params["Cost of Heat"] = LinearSample(m.fs.costing.heat_cost, 0, 0.15, nx)
    else:
        raise ValueError(f"{case_num} is not yet implemented")

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
    if len(sys.argv) == 1:
        print(
            "Usage: python multi_sweep.py case_number number_of_samples interpolate_nan_outputs"
        )
        print(
            f"Results will be written to {os.path.dirname(os.path.abspath(__file__))}"
        )
    else:
        results, sweep_params, m = run_analysis(*sys.argv[1:])
