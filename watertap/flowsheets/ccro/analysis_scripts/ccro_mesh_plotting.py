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
    zdata = list()
    for d, key in dm.keys():
        if d[-1] == 50:
            continue
        if key == xvar:
            for da in dm[d, zvar].data:
                zdata.append(da)
            for da in dm[d, xvar].data:
                xticks.append(da)
                yticks.append(d[-1])
        # if key == yvar:
        #     print("yticks", dm[d, yvar].data)
    print(input_maps)
    fig_init = dict(width=3.25, height=3.25, nrows=1, ncols=1, save_data=True)

    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    fig.plot_map(
        xdata=np.array(xticks),
        ydata=np.array(yticks),
        zdata=np.array(zdata),
        ax_idx=0,
        fix_nans=False,
    )

    fig.set_axis_ticklabels(
        xlabel=xvar,
        ylabel=yvar,
        ax_idx=0,
        xticklabels=np.unique(np.array(xticks)),
        yticklabels=np.unique(np.array(yticks)),
    )
    fig.show()

    fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
    fig.save_fig(name=fig_save)


if __name__ == "__main__":
    ccro_mesh_plotting()
