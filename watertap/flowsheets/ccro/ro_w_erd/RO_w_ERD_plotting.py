import os
import pprint
import numpy as np
from collections import defaultdict

from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.fig_generator import FigureGenerator

here = os.path.dirname(os.path.abspath(__file__))


def plot_ro_w_erd_sweep():

    n_pump = 1
    # sweep_file = f"{here}/output/ro_w_erd_analysisType_two_stage_1_pump_salinity_recovery_sweep.h5"

    # n_pump = 2
    sweep_file = f"{here}/output/ro_w_erd_analysisType_two_stage_{n_pump}_pump_salinity_recovery_sweep.h5"

    dm = PsDataManager(sweep_file)

    dm.register_data_key("fs.system_recovery", "System Recovery (%)", "%")
    dm.register_data_key("fs.costing.LCOW", "LCOW")
    dm.register_data_key("fs.costing.SEC", "SEC")

    for n in range(1, n_pump + 1):
        dm.register_data_key(
            f"fs.stage[{n}].flux",
            f"Stage {n} Flux",
            "L/m^2/h",
        )
        dm.register_data_key(
            f"fs.stage[{n}].RO.area",
            f"Stage {n} Area",
            "m^2",
        )
        try:
            dm.register_data_key(
                f"fs.stage[{n}].pump.control_volume.properties_out[0.0].pressure",
                f"Pump {n} Pressure",
                "bar",
            )
        except KeyError:
            print(f"Pump {n} Pressure not found in data manager, skipping.")

    dm.load_data()
    # dm.display()

    xvar = "System Recovery (%)"
    yvar = "Salinity (g/L)"

    fig_init = dict(width=3.25, height=3.25, nrows=1, ncols=1)

    input_maps = defaultdict(list)
    xticks = list()
    yticks = list()

    for zvar in [
        "LCOW",
        "SEC",
        "Pump 1 Pressure",
        "Stage 1 Area",
        "Stage 1 Flux",
        "Pump 2 Pressure",
        "Stage 2 Area",
        "Stage 2 Flux",
    ]:
        try:
            for d, key in dm.keys():
                if key == xvar:
                    input_maps[zvar].append(dm[d, zvar].data)
                    xticks = dm[d, xvar].data
                    yticks.append(d[-1] * 1000)

        except KeyError:
            print(f"Key {zvar} not found in data manager")
            continue

        fig = FigureGenerator(save_data=True)
        fig.init_figure(**fig_init)

        input_map = np.array(input_maps[zvar])
        # pprint.pprint(input_map[5])
        # print(yticks)
        # assert False

        xdata = np.array(xticks)
        ydata = np.array(yticks)

        fig.plot_map(
            xdata=xdata, ydata=ydata, zdata=input_map, ax_idx=0, fix_nans=False
        )

        fig.set_axis_ticklabels(
            xlabel=xvar,
            ylabel=yvar,
            ax_idx=0,
            xticklabels=[int(x) for x in xticks],
            yticklabels=[int(y) for y in yticks],
            ylims=[min(yticks), max(yticks)],
            xlims=[min(xticks), max(xticks)],
        )

        fig.ax[0].set_title(f"Two Stage RO with {n_pump} Pumps\n{zvar}")

        fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
        fig.save_fig(name=fig_save)


if __name__ == "__main__":
    plot_ro_w_erd_sweep()
