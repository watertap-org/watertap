import os
import pprint
from collections import defaultdict
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator

here = os.path.dirname(os.path.abspath(__file__))


def ccro_mesh_plotting(water_type="BW"):
    # sweep_file = f"{here}/output/ccro_flow_sweep_analysisType_study_BGW_mesh_study_optimization_lcow.h5"
    sweep_file = f"{here}/output/ccro_bw_mesh_analysisType_{water_type}_mesh_study.h5"
    dm = PsDataManager(sweep_file)

    for t in range(40):
        dm.register_data_key(
            f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
            f"{t} RO recovery",
            "%",
            # directory=t,
        )

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
        if key == xvar:
            for da in dm[d, zvar].data:
                zdata.append(da)
            for da in dm[d, xvar].data:
                xticks.append(da)
                yticks.append(d[-1])
        # if key == yvar:
        #     print("yticks", dm[d, yvar].data)
    print(zdata)
    fig = FigureGenerator(save_data=True)
    fig.init_figure()

    fig.plot_map(
        xdata=np.array(xticks),
        ydata=np.array(yticks),
        zdata=np.array(zdata),
        fix_nans=False,
    )

    fig.set_axis_ticklabels(
        xlabel=xvar,
        ylabel=yvar,
        xticklabels=np.unique(np.array(xticks)),
        yticklabels=np.unique(np.array(yticks)),
    )

    fig_save = sweep_file.replace(".h5", f"_{zvar}.png")
    fig.save_fig(name=fig_save)

    fig.show()

def ccro_line_plotting(
    n_filt_time_steps=10,
    n_flush_time_steps=5,
    n_time_steps=None,
    xdict={},
    ydict={},
    fig_init={
        "width": 4,
        "height": 4,
    },
    linedict={},
    title=None,):


    sweep_file = f"{here}/output/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5"
    if n_time_steps is None:

        n_time_steps = n_filt_time_steps + n_flush_time_steps

    if "units" not in ydict.keys():
        ydict["units"] = None

    cycle_time = list(range(n_time_steps))

    dm = PsDataManager(sweep_file)

    dm.register_data_key(f"{xdict['var']}", xdict["label"], xdict["units"])
    dm.register_data_key(f"{ydict['var']}", ydict["label"], ydict["units"])

    dm.load_data()
    dm.display()

    ydata = dm[(ydict["label"])].data
    xdata = dm[(xdict["label"])].data

    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    fig.plot_line(
        # xdata=np.array(xdict["ticks"]),
        xdata=np.array(xdata),
        # ydata=np.array(ydict["ticks"]),
        ydata=np.array(ydata),
        # zdata=np.array(zdata),
        # fix_nans=False,
        **linedict,
    )
    fig.set_axis_ticklabels(
        # xticklabels=np.array(xdata),
        # yticklabels=np.array(ydata),
        xlabel=xdict["label"],
        ylabel=ydict["label"],
        ax_idx=0,
    )

    fig_save = sweep_file.replace(".h5", f"_LINE_{ydict['label']}.png")
    fig.save_fig(name=fig_save)
    # fig.show()
    # for d in dm.keys():
    #     print(d, dm[d].data)
    # for d, key in dm.keys():
        # if key == zdict["label"]:
        #     z = list()
        #     for da in dm[d, key].data:
        #         z.append(float(da))
        #         print(d, key, da)
            # ydata.append(z)



def ccro_map_plotting(
    n_filt_time_steps=10,
    n_flush_time_steps=5,
    n_time_steps=None,
    xdict={},
    ydict={},
    zdict={},
    fig_init={
        "width": 4,
        "height": 4,
    },
    mapdict={},
    title=None,
):
    """
    Plotting CCRO maps.
    xdict: Sweep variable for x-axis.
    ydict: Build variable for y-axis.
    zdict: Variable for z-axis (color).
    """

    sweep_file = f"{here}/output/ccro_flush_eff_recovery_analysisType_SW_flushing_eff_and_recovery_sweep.h5"
    if n_time_steps is None:

        n_time_steps = n_filt_time_steps + n_flush_time_steps

    if "units" not in zdict.keys():
        zdict["units"] = None

    cycle_time = list(range(n_time_steps))

    dm = PsDataManager(sweep_file)

    dm.register_data_key(f"{zdict['var']}", zdict["label"], zdict["units"])
    dm.register_data_key(f"{xdict['var']}", xdict["label"], xdict["units"])
    dm.register_data_key(f"{ydict['var']}", ydict["label"], ydict["units"])

    # if zdict["var"] is not None:
    #     for t in range(n_time_steps):
    #         dm.register_data_key(
    #             f"blocks[{t}].process.{zdict['var']}",
    #             f"{t} {zdict['label']}",
    #             zdict["units"],
    #         )

    dm.load_data()
    dm.display()

    zdata = list()
    for d, key in dm.keys():
        if key == zdict["label"]:
            z = list()
            for da in dm[d, key].data:
                z.append(float(da))
            zdata.append(z)

    fig = FigureGenerator(save_data=True)
    fig.init_figure(**fig_init)

    fig.plot_map(
        xdata=np.array(xdict["ticks"]),
        ydata=np.array(ydict["ticks"]),
        zdata=np.array(zdata),
        # fix_nans=False,
        **mapdict,
    )

    fig.set_axis_ticklabels(
        xticklabels=xdict["ticks"],
        yticklabels=ydict["ticks"],
        xlabel=xdict["label"],
        ylabel=ydict["label"],
        ax_idx=0,
    )

    if title is not None:
        ax = fig.get_axis(0)
        ax.set_title(title)

    fig_save = sweep_file.replace(".h5", f"_MAP_{zdict['label']}.png")
    fig.save_fig(name=fig_save)

    # fig.show()


if __name__ == "__main__":

    xdict = {
        "var": "flushing.flushing_efficiency",
        "ticks": [20, 27, 34, 41, 48, 55, 62, 69, 76, 83, 90],
        "label": "Flushing Efficiency (%)",
        "units": "%",
    }

    ydict = {
        "var": "overall_recovery",
        "ticks": [40, 45, 50, 55, 60],
        "label": "Water Recovery (%)",
        "units": "%",
    }

    zdict = {
        "var": "costing.LCOW",
        "label": "LCOW",
    }

    zdict = {
        "var": "costing.SEC",
        "label": "SEC",
    }
# f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
# f"blocks[{t}].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
    zdict = {
        "var": "blocks[0].process.fs.RO.area",
        "label": "RO Area",
        "units": "m^2",
    }
    zdict = {
        "var": "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]",
        "label": "RO Recovery",
        "units": "%",
    }
    zdict = {
        "var": "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "label": "Recycle Rate",
        "units": "L/s",
    }
    zdict = {
        "var": "filtration_ramp_rate",
        "label": "Filtration Ramp Rate",
        "units": "bar/min",
    }
    zdict = {
        "var": "total_cycle_time",
        "label": "Total Cycle Time",
        "units": "min",
    }

        # auto_sig_0_1=2,
        # auto_sig_1_10=1,
        # auto_sig_10_inf=0,
    mapdict = {"textfontsize": 8, "auto_sig_1_10": 2}

    fig_init = {
        "width": 4,
        "height": 3,
    }

    title = "Seawater CCRO\nFlushing Efficiency and Recovery Sweep"

    # ccro_map_plotting(zdict=zdict, xdict=xdict, ydict=ydict, mapdict=mapdict, fig_init=fig_init, title=title)

    ########################################################


    xdict = {
        "var": "overall_recovery",
        "ticks": [40, 45, 50, 55, 60],
        "label": "Water Recovery (%)",
        "units": "%",
    }
    ydict = {
        "var": "flushing.flushing_efficiency",
        "ticks": [20, 27, 34, 41, 48, 55, 62, 69, 76, 83, 90],
        "label": "Flushing Efficiency (%)",
        "units": "%",
    }
    ydict = {
        "var": "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "label": "Recycle Rate",
        "units": "L/s",
    }
    ydict = {
        "var": "filtration_ramp_rate",
        "label": "Filtration Ramp Rate",
        "units": "bar/min",
    }
    ydict = {
        "var": "costing.LCOW",
        "label": "LCOW",
    }
    ydict = {
        "var": "costing.SEC",
        "label": "SEC",
    }
    ydict = {
        "var": "total_cycle_time",
        "label": "Total Cycle Time",
        "units": "min",
    }
    ydict = {
        "var": "final_concentration",
        "label": "Permeate Concentration",
        "units": "mg/L",
    }

    linedict = {
        "marker": "o",
        "color": "blue",
    }


    ccro_line_plotting(xdict=xdict, ydict=ydict, fig_init=fig_init, title=title, linedict=linedict)

    # ccro_mesh_plotting("SW")

    # ccro_mesh_plotting("BW")
