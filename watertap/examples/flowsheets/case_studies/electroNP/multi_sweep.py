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
from watertap.tools.parameter_sweep import (
    LinearSample,
    parameter_sweep,
)
import watertap.examples.flowsheets.case_studies.electroNP.electroNP_flowsheet as electroNP_flowsheet
from pyomo.environ import units as pyunits
import os
from watertap.tools.parameter_sweep.sweep_visualizer import line_plot, contour_plot
import matplotlib.pyplot as plt


def set_up_sensitivity(m):
    outputs = {}
    # optimize_kwargs = {"fail_flag": False}
    optimize_kwargs = {}
    opt_function = electroNP_flowsheet.solve

    # create outputs
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["ElectroNP Capital Cost"] = m.fs.electroNP.costing.capital_cost
    outputs["Electricity"] = pyunits.convert(
        m.fs.costing.aggregate_flow_costs["electricity"] / m.fs.AD.inlet.flow_vol[0],
        to_units=pyunits.USD_2018 / pyunits.m**3,
    )

    return outputs, optimize_kwargs, opt_function


def run_analysis(case_num, nx, interpolate_nan_outputs=True, output_filename=None):

    if output_filename is None:
        output_filename = "sensitivity_" + str(case_num) + ".csv"

    m = electroNP_flowsheet.build_flowsheet()[0]

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)

    sweep_params = {}
    if case_num == 1:
        # sensitivity analysis
        m.fs.costing.electroNP.sizing_cost.unfix()
        sweep_params["sizing_cost"] = LinearSample(
            m.fs.costing.electroNP.sizing_cost, 1, 30, nx
        )
    elif case_num == 2:
        m.fs.costing.electroNP.sizing_cost.unfix()
        sweep_params["sizing_cost"] = LinearSample(
            m.fs.costing.electroNP.sizing_cost, 0, 30, nx
        )
        m.fs.costing.electroNP.HRT.unfix()
        sweep_params["HRT"] = LinearSample(m.fs.costing.electroNP.HRT, 0.1, 5, nx)
    elif case_num == 3:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0.04, 2.5, nx
        )
    elif case_num == 4:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0.02, 2.2, nx
        )
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.02, 0.3, nx
        )
    elif case_num == 5:
        sweep_params["electricity_cost"] = LinearSample(
            m.fs.costing.electricity_cost, 0.06, 0.1, nx
        )
        sweep_params["MgCl2_cost"] = LinearSample(
            m.fs.costing.magnesium_chloride_cost, 0.06, 0.1, nx
        )
    elif case_num == 6:
        sweep_params["energy_consumption"] = LinearSample(
            m.fs.electroNP.energy_electric_flow_mass, 0, 2.5, nx
        )
        sweep_params["MgCl2_dosage"] = LinearSample(
            m.fs.electroNP.magnesium_chloride_dosage, 0.1, 0.5, nx
        )
    else:
        raise ValueError(f"{case_num} is not yet implemented")

    # output_filename = "sensitivity_" + str(case_num) + ".csv"

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
    data_file = merge_path("interpolated_sensitivity_" + str(case_num) + ".csv")
    if plot_type == "line":
        fig, ax = line_plot(data_file, xlabel, ylabel, xunit, yunit)

    elif plot_type == "contour":
        fig, ax = contour_plot(
            data_file,
            xlabel,
            ylabel,
            zlabel,
            xunit,
            yunit,
            zunit,
            levels=levels,
            cmap=cmap,
            isolines=isolines,
        )
    else:
        raise ValueError("Plot type not yet implemented")
    return fig, ax


def main(case_num, nx=11, interpolate_nan_outputs=True):
    # when from the command line
    case_num = int(case_num)
    nx = int(nx)
    interpolate_nan_outputs = bool(interpolate_nan_outputs)

    # comm, rank, num_procs = _init_mpi()

    global_results, sweep_params, m = run_analysis(
        case_num, nx, interpolate_nan_outputs
    )
    # print(global_results)

    # visualize results

    # case 1
    # fig, ax = visualize_results(
    #     case_num,
    #     plot_type="line",
    #     xlabel="# sizing_cost",
    #     ylabel="ElectroNP Capital Cost",
    # )
    # # ax.plot(1.25, 127.792, 'ro')
    # ax.set_xlabel("Sizing Cost ($/m3)")
    # ax.set_ylabel("ElectroNP Capital Cost ($/(m3/hr)")

    # case 2
    # fig, ax = visualize_results(
    #     case_num,
    #     plot_type="contour",
    #     xlabel="# sizing_cost",
    #     ylabel="HRT",
    #     zlabel="ElectroNP Capital Cost",
    #     isolines=[2000, 5000],
    #     cmap="GnBu",
    # )
    # ax.plot(1.25, 1.3333, "ko")
    # ax.set_xlabel("Sizing Cost ($/m3)")
    # ax.set_ylabel("HRT (h)")
    # ax.set_title("ElectroNP Capital Cost ($/(m3/hr)")

    # case 3
    # fig, ax = visualize_results(
    #     case_num,
    #     plot_type="line",
    #     xlabel="# energy_consumption",
    #     ylabel="LCOW",
    # )
    # # ax.plot(0.044, 6.0104, 'ro')
    # # plt.axvline(0.044, ls='--')
    # ax.set_xlabel("Electricity Intensity (kWh/kg P)")
    # ax.set_ylabel("LCOW ($/m3)")

    # case 4
    fig, ax = visualize_results(
        case_num,
        plot_type="contour",
        xlabel="# energy_consumption",
        ylabel="electricity_cost",
        zlabel="LCOW",
        isolines=[2000, 5000],
        cmap="GnBu",
    )
    ax.plot(0.044, 0.08, "ko")
    ax.set_xlim([0, 2])
    ax.set_ylim([0.02, 0.2])
    ax.set_xlabel("Electricity Intensity (kWh/kg P)")
    ax.set_ylabel("Electricity Cost ($/kWh)")
    ax.set_title("LCOW ($/m3)")

    # case 5
    # fig, ax = visualize_results(
    #     case_num,
    #     plot_type="contour",
    #     xlabel="# electricity_cost",
    #     ylabel="MgCl2_cost",
    #     zlabel="LCOW",
    #     cmap="GnBu",
    # )
    # ax.plot(0.08, 0.0786, 'ko')
    # ax.set_ylim([0.07, 0.1])
    # ax.set_xlim([0.06, 0.1])
    # ax.set_xlabel("Electricity Cost ($/kWh)")
    # ax.set_ylabel("MgCl2 Cost ($/kg P)")
    # ax.set_title("LCOW ($/m3)")

    # case 6
    # fig, ax = visualize_results(
    #     case_num,
    #     plot_type="contour",
    #     xlabel="# energy_consumption",
    #     ylabel="MgCl2_dosage",
    #     zlabel="LCOW",
    #     cmap="GnBu",
    # )
    # ax.plot(0.044, 0.388, 'ko')
    # ax.set_ylim([0.2, 0.5])
    # ax.set_xlabel("Electricity Intensity (kWh/kg P)")
    # ax.set_ylabel("MgCl2 Dosage (kg MgCl2/kg P)")
    # ax.set_title("LCOW ($/m3)")

    return global_results, m


if __name__ == "__main__":
    results, model = main(case_num=4, nx=17)
    # results, sweep_params, m = run_analysis()
