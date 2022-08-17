import sys
import os
import time

from pyomo.environ import Constraint
from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
from watertap.tools.sweep_visualizer import line_plot, contour_plot
import watertap.examples.flowsheets.uvaop.sim_uv_aop as uvaop


def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {"check_termination": True}
    opt_function = uvaop.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    # outputs["Total_Cost"] = m.fs.costing.total_annualized_cost

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True):
    m = uvaop.main()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}

    if case_num == 1:
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].unfix()
        sweep_params["EEO"] = LinearSample(m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"], 0.1, 5, nx)

    elif case_num == 2:
        m.fs.unit.uv_intensity.unfix()
        sweep_params["uv_intensity"] = LinearSample(m.fs.unit.uv_intensity, 0.5, 2, nx)

    else:
        raise ValueError("case_num = %d not recognized." % (case_num))

    output_filename = merge_path("sensitivity_" + str(case_num) + ".csv")
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


def merge_path(filename):
    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), filename)
    return source_path


def visualize_results(
    case_num,
    plot_type,
    xlabel,
    ylabel,
    xunit=None,
    yunit=None,
    zlabel=None,
    zunit=None,
    levels=None,
    cmap=None,
    isolines=None,
):
    data_file = merge_path("sensitivity_" + str(case_num) + ".csv")
    if plot_type == "line":
        fig, ax = line_plot(data_file, xlabel, ylabel, xunit, yunit)

    elif plot_type == "contour":
        fig, ax = contour_plot(
            data_file, xlabel, ylabel, zlabel, xunit, yunit, zunit, isolines=isolines
        )
    else:
        raise ValueError("Plot type not yet implemented")
    return fig


def main(case_num=1, nx=10, interpolate_nan_outputs=False):
    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    # comm, rank, num_procs = _init_mpi()

    global_results, sweep_params, m = run_analysis(
        case_num, nx, interpolate_nan_outputs
    )
    print(global_results)

    # visualize results
    data_file = merge_path("sensitivity_" + str(case_num) + ".csv")
    line_plot(data_file, "# EEO", "LCOW")
    # contour_plot(
    #     data_file, "# dye_cost", "waste_disposal", "LCOW", isolines=[1, 2], cmap="GnBu"
    # )

    return global_results, m


if __name__ == "__main__":
    results, model = main()