#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as RO
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    ERDtype,
)


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"fail_flag": False}
    opt_function = RO.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW

    return outputs, optimize_kwargs, opt_function


def run_analysis(
    case_num=1, nx=11, interpolate_nan_outputs=True, withERD=True, output_filename=None
):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)
    withERD = bool(withERD)

    # select flowsheet configuration
    if withERD:
        m = RO.main(erd_type=ERDtype.pump_as_turbine)[0]
    else:
        m = RO.main(erd_type=ERDtype.pressure_exchanger)[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    # choose parameter sweep from case structure
    sweep_params = {}

    if case_num == 1:
        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].unfix()
        m.fs.feed.properties[0].flow_vol_phase["Liq"].unfix()

        sweep_params["salt_concentration"] = LinearSample(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"], 0.001, 0.01, nx
        )
        sweep_params["flowrate"] = LinearSample(
            m.fs.feed.properties[0].flow_vol_phase["Liq"],
            0.03,
            0.11,
            nx,
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

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
