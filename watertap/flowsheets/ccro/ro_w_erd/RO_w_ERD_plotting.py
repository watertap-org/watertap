import pprint
from collections import defaultdict
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator

if __name__ == "__main__":

    sweep_file = "/Users/ksitterl/Documents/Python/watertap/watertap/watertap/flowsheets/ccro/ro_w_erd/output/ro_w_erd_analysisType_two_stage_1_pump_salinity_recovery_sweep.h5"
    dm = PsDataManager(sweep_file)

    dm.register_data_key("fs.system_recovery", "System Recovery (%)", "%")
    dm.register_data_key("fs.costing.LCOW", "LCOW")
    dm.register_data_key("fs.costing.SEC", "SEC")
    dm.register_data_key(
        "fs.stage[1].pump.control_volume.properties_out[0.0].pressure",
        "Pump 1 Pressure",
        "bar",
    )
    # dm.register_data_key(
    #     "fs.stage[2].pump.control_volume.properties_out[0.0].pressure",
    #     "Pump 2 Pressure",
    #     "bar",
    # )
    dm.load_data()
    dm.display()

    xvar = "System Recovery (%)"
    yvar = "Salinity (g/L)"

    fig_init = dict(width=4, height=4, nrows=1, ncols=1, save_data=True)

    input_maps = defaultdict(list)
    xticks = list()
    yticks = list()

    for zvar in ["LCOW", "SEC", "Pump 1 Pressure", "Pump 2 Pressure"]:
        try:
            
            for d, key in dm.keys():
                if key == xvar:

                    input_maps[zvar].append(dm[d, zvar].data)
                    xticks = dm[d, xvar].data
                    yticks.append(d[-1] * 1000)
        except KeyError:
            print(f"Key {zvar} not found in data manager")

    # print(input_maps)

    ###############################################
    zvar = "LCOW"
    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    input_map = np.array(input_maps[zvar])
    print(input_map)

    xdata = np.array(xticks)
    ydata = np.array(yticks[::-1])

    fig.plot_map(xdata=xdata, ydata=ydata, zdata=input_map, ax_idx=0, fix_nans=False)

    fig.set_axis_ticklabels(
        xlabel=xvar,
        ylabel=yvar,
        ax_idx=0,
        xticklabels=[int(x) for x in xticks],
        yticklabels=[int(y) for y in yticks],
        ylims=[min(yticks), max(yticks)],
        xlims=[min(xticks), max(xticks)],
    )

    fig.ax[0].set_title(zvar)

    fig_save = sweep_file.replace(".h5", f".{zvar}.png")
    fig.save_fig(name=fig_save)

    ###############################################

    zvar = "SEC"
    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    input_map = np.array(input_maps[zvar])

    xdata = np.array(xticks)
    ydata = np.array(yticks[::-1])

    fig.plot_map(xdata=xdata, ydata=ydata, zdata=input_map, ax_idx=0, fix_nans=True)

    fig.set_axis_ticklabels(
        xlabel=xvar,
        ylabel=yvar,
        ax_idx=0,
        xticklabels=[int(x) for x in xticks],
        yticklabels=[int(y) for y in yticks],
        ylims=[min(yticks), max(yticks)],
        xlims=[min(xticks), max(xticks)],
    )

    fig.ax[0].set_title(zvar)

    fig_save = sweep_file.replace(".h5", f".{zvar}.png")
    fig.save_fig(name=fig_save)

    ###############################################
    zvar = "Pump 1 Pressure"
    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    input_map = np.array(input_maps[zvar])

    xdata = np.array(xticks)
    ydata = np.array(yticks)

    fig.plot_map(xdata=xdata, ydata=ydata, zdata=input_map, ax_idx=0, fix_nans=True)

    fig.set_axis_ticklabels(
        xlabel=xvar,
        ylabel=yvar,
        ax_idx=0,
        xticklabels=[int(x) for x in xticks],
        yticklabels=[int(y) for y in yticks],
        ylims=[min(yticks), max(yticks)],
        xlims=[min(xticks), max(xticks)],
    )

    fig.ax[0].set_title(zvar)

    fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
    fig.save_fig(name=fig_save)

    ###############################################
    # zvar = "Pump 2 Pressure"
    # fig = FigureGenerator(save_data=True)
    # fig.init_figure(**fig_init)

    # input_map = np.array(input_maps[zvar])

    # xdata = np.array(xticks)
    # ydata = np.array(yticks)

    # fig.plot_map(xdata=xdata, ydata=ydata, zdata=input_map, ax_idx=0, fix_nans=True)

    # fig.set_axis_ticklabels(
    #     xlabel=xvar,
    #     ylabel=yvar,
    #     ax_idx=0,
    #     xticklabels=[int(x) for x in xticks],
    #     yticklabels=[int(y) for y in yticks],
    #     ylims=[min(yticks), max(yticks)],
    #     xlims=[min(xticks), max(xticks)],
    # )

    # fig.ax[0].set_title(zvar)

    # fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
    # fig.save_fig(name=fig_save)
