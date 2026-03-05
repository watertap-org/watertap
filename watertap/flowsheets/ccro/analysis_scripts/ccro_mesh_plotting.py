import os
import pprint
from collections import defaultdict
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator

here = os.path.dirname(os.path.abspath(__file__))

def ccro_mesh_plotting():
    # sweep_file = f"{here}/output/ccro_flow_sweep_analysisType_study_BGW_mesh_study_optimization_lcow.h5"
    sweep_file = f"{here}/output/ccro_bw_mesh_analysisType_BW_mesh_study.h5"
    dm = PsDataManager(sweep_file)

    # for t in range(40):
    #     dm.register_data_key(
    #         f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
    #         f"RO recovery",
    #         "%",
    #         # directory=t,
    #     )

    dm.register_data_key("costing.LCOW", "LCOW")
    dm.register_data_key("costing.SEC", "SEC")
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.load_data()
    dm.display()

    xvar = "Water recovery"
    yvar = "Time Steps"
    zvar = "LCOW"

    input_maps = defaultdict(list)
    xticks = list()
    yticks = list()
    for d, key in dm.keys():
        if d[-1] == 50:
            continue
        print(d, key)
        if key == xvar:
            print("xticks", xticks)
            print("yticks", yticks)
            input_maps[zvar].append(dm[d, zvar].data)
            xticks = dm[d, xvar].data
            yticks.append(d[-1])
        # if key == yvar:
        #     print("yticks", dm[d, yvar].data)


    fig_init = dict(width=3.25, height=3.25, nrows=1, ncols=1, save_data=True)

    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)
    # import matplotlib.pyplot as plt

    # fig, ax = plt.subplots()

    input_map = np.array(input_maps[zvar])
    pprint.pprint(input_map)
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
    fig.show()

    fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
    fig.save_fig(name=fig_save)


if __name__ == "__main__":
    ccro_mesh_plotting()