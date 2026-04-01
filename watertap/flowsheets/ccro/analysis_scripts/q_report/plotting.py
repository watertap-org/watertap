import itertools
import pprint
import os
import numpy as np
from collections import OrderedDict
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)
from matplotlib.patches import Patch
from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter
from watertap.flowsheets.ccro.analysis_scripts.q_report import get_data as gd

here = os.path.dirname(os.path.abspath(__file__))
par_dir = os.path.dirname(here)
cases = ["Brackish water", "Seawater", "Produced water"]

wc_ac = {
    "Brackish water": "BW",
    "Seawater": "SW",
    "Produced water": "PW",
}

label_dict = {
    "LCOW": "LCOW ($\$$/m$^3$)",
    "SEC": "SEC (kWh/m$^3$)",
    "Permeate concentration": "Perm conc (g/L)",
    "Flux": "Flux (LMH)",
    # "CAPEX": "CAPEX (\\$10$^3$)",
    # "OPEX": "OPEX (\\$10$^3$ yr$^{-1}$)",
    "CAPEX": "CAPEX (kUSD)",
    "OPEX": "OPEX (kUSD/yr)",
    "Area": "RO area (m$^2$)",
    "Membrane Area": "Membrane Area (m$^2$)",
    "Specific area": "Spec area (m$^2$/(L/hr))",
    "Pressure": "RO Pressure (bar)",
    "Pump work": "Pump Work (kW)",
    "Pump size": "Pump Size (kW)",
    "Inlet concentration": "Inlet concentration (g/L)",
    "Initial concentration": "Initial concentration (g/L)",
    ("costing", "total", "LCOW_opex"): "LCOW OPEX ($\$$/m$^3$)",
    ("costing", "total", "LCOW_capex"): "LCOW CAPEX ($\$$/m$^3$)",
    "OPEX/CAPEX Ratio": "OPEX/CAPEX Ratio",
    "OPEX LCOW Fraction": "LCOW OPEX Fraction",
    "CAPEX LCOW Fraction": "LCOW CAPEX Fraction",
    "Specific Membrane Area": "Specific\nMembrane Area\n(m$^2$/(L/hr))",
    "Total cycle time": "Total Cycle Time (min)",
    "Cycle time ratio": "Cycle Time Ratio",
    "Recycle rate": "Recycle Rate (L/s)",
    "Ramp rate": "Ramp Rate (bar/min)",
    "Recycle ratio": "Recycle Ratio (%)",
    "Water recovery": "Water Recovery (%)",
    "Single Pass Recovery": "Single Pass Recovery (%)",
    "Recycle loop concentration": "Recycle Loop Concentration (g/L)",
    "Flushing efficiency": "Flushing Efficiency (%)",
}

yticks_dict = {
    "Brackish water": {
        "LCOW": [0, 0.1, 0.2, 0.3, 0.4],
        "SEC": [0, 0.5, 1, 1.5, 2, 2.5],
        "Permeate concentration": [0, 0.2, 0.4, 0.6],
        "Flux": [0, 15, 30, 45, 60, 75],
        "CAPEX": [0, 20, 40, 60],
        "OPEX": [0, 2, 4, 6, 8],
        "Pump work": [0, 2, 4, 6, 8, 10],
        "Pump size": [0, 2, 4, 6, 8, 10, 12],
        "Area": [0, 50, 100, 150],
        "Membrane Area": [0, 75, 150, 225, 300],
        "Pressure": [0, 25, 50, 75, 100],
        "Inlet concentration": [0, 10, 20, 30],
        ("costing", "total", "LCOW_opex"): [0, 0.05, 0.1, 0.15, 0.2, 0.25],
        ("costing", "total", "LCOW_capex"): [0, 0.05, 0.1, 0.15, 0.2, 0.25],
        "OPEX/CAPEX Ratio": [0, 0.05, 0.1, 0.15, 0.2, 0.25],
        "OPEX LCOW Fraction": [0, 0.25, 0.5, 0.75],
        "CAPEX LCOW Fraction": [0, 0.15, 0.3, 0.45, 0.6],
        "Specific Membrane Area": [0, 0.02, 0.04, 0.06, 0.08, 0.1],
        "Single Pass Recovery": [0, 15, 30, 45, 60],
        "Initial concentration": [0, 10, 20, 30],
        "Ramp rate": [0, 1, 2, 3, 4, 5],
        "Recycle ratio": [0, 1, 2, 3],
        "Recycle rate": [0, 1, 2, 3],
        "Cycle time ratio": [0, 25, 50, 75, 100],
        "Final concentration": [0, 15, 30, 45],
        "Final Pressure": [0, 25, 50, 75, 100],
        "Recycle loop concentration": [0, 30, 60, 90],
        "Flushing efficiency": [0, 25, 50, 75, 100],
    },
    "Seawater": {
        "LCOW": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
        # "SEC": [0, 1, 2, 3],
        "SEC": [0, 1, 2, 3, 4],
        "Permeate concentration": [0, 0.2, 0.4, 0.6],
        "Flux": [
            0,
            10,
            20,
            30,
            40,
        ],
        "CAPEX": [0, 25, 50, 75, 100],
        "OPEX": [0, 2, 4, 6, 8, 10],
        "Pump work": [0, 2, 4, 6, 8, 10],
        "Pump size": [0, 4, 8, 12, 16],
        "Area": [0, 75, 150, 225, 300],
        "Membrane Area": [0, 75, 150, 225, 300],
        "Pressure": [0, 40, 80, 120],
        ("costing", "total", "LCOW_opex"): [0, 0.2, 0.4, 0.6],
        ("costing", "total", "LCOW_capex"): [0, 0.2, 0.4, 0.6],
        "OPEX/CAPEX Ratio": [0, 0.05, 0.1, 0.15, 0.2, 0.25],
        "OPEX LCOW Fraction": [0, 0.25, 0.5, 0.75],
        "CAPEX LCOW Fraction": [0, 0.15, 0.3, 0.45, 0.6],
        "Specific Membrane Area": [0, 0.02, 0.04, 0.06, 0.08, 0.1],
        "Single Pass Recovery": [0, 15, 30, 45, 60],
        "Inlet concentration": [0, 10, 20, 30, 40, 50],
        "Initial concentration": [0, 10, 20, 30, 40, 50],
        "Ramp rate": [0, 2, 4, 6, 8, 10],
        "Recycle ratio": [0, 2.5, 5, 7.5, 10],
        "Recycle rate": [0, 2.5, 5, 7.5, 10],
        "Cycle time ratio": [0, 25, 50, 75, 100],
        "Final concentration": [0, 25, 50, 75, 100],
        "Final Pressure": [0, 40, 80, 120],
        "Recycle loop concentration": [0, 30, 60, 90, 120],
        "Flushing efficiency": [0, 25, 50, 75, 100],
    },
    "Produced water": {
        "LCOW": [
            0,
            0.5,
            1.0,
            1.5,
            2.0,
        ],
        "SEC": [0, 2, 4, 6, 8],
        "Permeate concentration": [0, 0.2, 0.4, 0.6],
        "Flux": [
            0,
            10,
            20,
            30,
            40,
        ],
        "CAPEX": [0, 30, 60, 90, 120],
        "OPEX": [0, 4, 8, 12, 16, 20],
        "Pump work": [0, 4, 8, 12, 16, 20],
        "Pump size": [0, 5, 10, 15, 20, 25],
        "Area": [0, 40, 80, 120, 160],
        "Membrane Area": [0, 40, 80, 120, 160, 200],
        "Pressure": [0, 50, 100, 150, 200],
        "Inlet concentration": [0, 40, 80, 120, 160],
        "Initial concentration": [0, 40, 80, 120, 160],
        ("costing", "total", "LCOW_opex"): [0, 0.25, 0.5, 0.75, 1],
        ("costing", "total", "LCOW_capex"): [0, 0.25, 0.5, 0.75, 1],
        "OPEX/CAPEX Ratio": [0, 0.05, 0.1, 0.15, 0.2, 0.25],
        "OPEX LCOW Fraction": [0, 0.25, 0.5, 0.75],
        "CAPEX LCOW Fraction": [0, 0.15, 0.3, 0.45, 0.6],
        "Specific Membrane Area": [0, 0.02, 0.04, 0.06, 0.08, 0.1],
        "Single Pass Recovery": [0, 15, 30, 45, 60],
        "Ramp rate": [0, 20, 40, 60],
        "Recycle ratio": [0, 5, 10, 15, 20, 25],
        "Recycle rate": [0, 5, 10, 15, 20, 25],
        "Cycle time ratio": [0, 25, 50, 75, 100],
        "Final concentration": [0, 40, 80, 120, 160],
        "Final Pressure": [0, 50, 100, 150, 200],
        "Recycle loop concentration": [0, 40, 80, 120, 160, 200],
        "Flushing efficiency": [0, 25, 50, 75, 100],
    },
}


xticks_dict = {
    "Brackish water": {
        "Water recovery": [75, 80, 85, 90, 95],
        "Recycle rate": [5, 20, 35],
    },
    "Seawater": {"Water recovery": [40, 45, 50, 55, 60], "Recycle rate": [5, 20, 35]},
    "Produced water": {
        "Water recovery": [20, 30, 40, 50, 55],
        "Recycle rate": [5, 20, 35],
    },
}

line_colors = [
    # "#a6cee3",
    "#1f78b4",
    # "#b2df8a",
    "#33a02c",
    # "#fb9a99",
    "#e31a1c",
    # "#fdbf6f",
    "#ff7f00",
    # "#cab2d6",
    "#6a3d9a",
    "k",
]
units = [
    "Feed pump",
    "Recycle pump",
    "RO",
    "Conduit",
    "Pumps + ERD",
    "Stage 1",
    "Stage 2",
]

color_dict = dict(zip(units, line_colors))


def panel_3_by_2_stack(fraction=False):

    # dm_ccro = gd.get_ccro_data()

    # dm_ss = gd.get_ss_data()
    nrows = 3
    ncols = 2
    w = 2.25
    width = ncols * w
    height = nrows * w

    init_figure = dict(
        width=width,
        height=height,
        nrows=nrows,
        ncols=ncols,
        sharex=False,
        sharey=False,
    )

    fig = FigureGenerator()
    fig.init_figure(**init_figure)

    bdp = BreakDownPlotter(
        save_folder=f"{here}/figs",
        save_name=f"test",
        fig=fig,
        color_dict=color_dict,
    )
    # if fraction:

    #     bdp.define_hatch_groups(
    #         {
    #             "LCOW_opex_fraction": {
    #                 "label": "OPEX",
    #                 "hatch": "",
    #                 "color": "white",
    #                 "alpha": 0.85,
    #             },
    #             "LCOW_capex_fraction": {"label": "CAPEX", "hatch": "\\\\", "color": "white"},
    #         }
    #     )
    # else:

    bdp.define_hatch_groups(
        {
            "LCOW_opex": {
                "label": "OPEX",
                "hatch": "",
                "color": "white",
                "alpha": 0.85,
            },
            "LCOW_capex": {"label": "CAPEX", "hatch": "\\\\", "color": "white"},
        }
    )

    for i, water_case in enumerate(cases):

        axis_options = {"xticks": xticks_dict[water_case]}

        axs = list()

        ######################################## CCRO

        ax_idx = (i, 0)  # CCRO on left
        axs.append(bdp.fig.get_axis(ax_idx))
        # axis_options["ax_idx"] = ax_idx

        dm_recov.select_data(water_case)
        wr = dm_recov.get_selected_data()
        bdp.set_selected_data(wr)

        bdp.define_area_groups(
            {
                "Feed pump": {},
                "Recycle pump": {},
                "RO": {},
                "Conduit": {},
            }
        )
        if i != len(cases) - 1:
            axis_options["ylabel"] = label_dict["LCOW"]
            axis_options["yticks"] = yticks_dict[water_case]["LCOW"]
            axis_options["xlabel"] = ""
        else:
            axis_options["ylabel"] = label_dict["LCOW"]
            axis_options["yticks"] = yticks_dict[water_case]["LCOW"]
            axis_options["xlabel"] = "Water recovery (%)"

        bdp.plotbreakdown(
            xdata="Water recovery",
            ydata="costing",
            axis_options=axis_options,
            ax_idx=ax_idx,
        )

        ######################################## STEADY STATE

        ax_idx = (i, 1)  # steady state on right
        axs.append(bdp.fig.get_axis(ax_idx))

        axis_options["yticks"] = yticks_dict[water_case]["LCOW"]

        if water_case == "Brackish water":
            ro_config = "2_stage_1_pump"
            add_erd = "False"
        elif water_case == "Seawater":
            ro_config = "2_stage_1_pump"
            add_erd = "True"
        elif water_case == "Produced water":
            ro_config = "2_stage_2_pump"
            add_erd = "True"

        ro_args = [("add_erd", add_erd), ("stage_sim_cases", ro_config)]
        dm_ss.select_data(
            (water_case, *ro_args),
        )

        wr = dm_ss.get_selected_data()
        bdp.set_selected_data(wr)

        bdp.define_area_groups(
            {
                "Pumps + ERD": {},
                "Stage 1": {},
                "Stage 2": {},
            }
        )
        axis_options["ylabel"] = ""
        axis_options["yticks"] = yticks_dict[water_case]["LCOW"]

        bdp.plotbreakdown(
            xdata="Water recovery",
            ydata="costing",
            axis_options=axis_options,
            ax_idx=ax_idx,
        )

        ys = [ax.get_ylim() for ax in axs]
        y_min = min([y[0] for y in ys])
        y_max = max([y[1] for y in ys]) * 1.05

        for ax in axs:
            ax.set_ylim(y_min, y_max)

    fig.fig.tight_layout()

    # add legend

    ax1 = fig.ax[2][0]
    ax2 = fig.ax[2][1]

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    combined = OrderedDict(zip(labels1 + labels2, handles1 + handles2))
    combined = {k: v for k, v in combined.items() if k not in ["CAPEX", "OPEX"]}

    fig.fig.legend(
        combined.values(),
        combined.keys(),
        loc="upper center",
        ncol=4,
        fontsize=8,
        bbox_to_anchor=(0.54, 1.08),  # position at top center
    )

    ycol = "LCOW"
    if fraction:
        fig.save_fig(name=f"{here}/figs/{ycol}_fraction_tile_plot_stack")
    else:
        fig.save_fig(name=f"{here}/figs/{ycol}_tile_plot_stack")

    fig.show()

    return fig, bdp


def plot_ccro_line(
    dm,
    fig=None,
    ax_idx=0,
    water_case="Brackish water",
    xcol="Water recovery",
    ycol="LCOW",
    plot_options={},
    axis_options={},
):
    if fig is None:
        init_figure = {
            "width": 2,
            "height": 2,
            "nrows": 1,
            "ncols": 1,
        }

        fig = FigureGenerator()
        fig.init_figure(**init_figure)

    fig.plot_line(
        dm[water_case, xcol], dm[water_case, ycol], ax_idx=ax_idx, **plot_options
    )

    if axis_options == {}:
        axis_options = {
            "xlabel": label_dict.get(xcol, None),
            "ylabel": label_dict.get(ycol, None),
            "xticks": xticks_dict[water_case].get(xcol, None),
            "yticks": yticks_dict[water_case].get(ycol, None),
            "ax_idx": ax_idx,
        }

    fig.set_axis(**axis_options)

    # fig.save_fig(name=f"{here}/figs/{ycol}_vs_{xcol}_line_plot_CCRO")

    return fig


def plot_ss_line(
    dm,
    add_erd="False",
    stage_sim_cases="2_stage_1_pump",
    fig=None,
    ax_idx=0,
    water_case="Brackish water",
    xcol="Water recovery",
    ycol="LCOW",
    plot_by_stage=False,
    plot_options={},
    plot2_options={},
    axis_options={},
    leg_labels=[],
):
    if fig is None:
        init_figure = {
            "width": 2,
            "height": 2,
            "nrows": 1,
            "ncols": 1,
        }

        fig = FigureGenerator()
        fig.init_figure(**init_figure)

    if leg_labels == []:
        leg_labels = [water_case]

    ro_args = [("add_erd", add_erd), ("stage_sim_cases", stage_sim_cases)]

    if plot_by_stage:
        # for stage in [1, 2]:
        fig.plot_line(
            dm[water_case, *ro_args, xcol],
            dm[water_case, *ro_args, ("RO 1", ycol)],
            ax_idx=ax_idx,
            **plot_options,
        )
        fig.plot_line(
            dm[water_case, *ro_args, xcol],
            dm[water_case, *ro_args, ("RO 2", ycol)],
            ax_idx=ax_idx,
            **plot2_options,
        )
    else:
        fig.plot_line(
            dm[water_case, *ro_args, xcol],
            dm[water_case, *ro_args, ycol],
            ax_idx=ax_idx,
            **plot_options,
        )

    if axis_options == {}:
        axis_options = {
            "xlabel": label_dict.get(xcol, None),
            "ylabel": label_dict.get(ycol, None),
            "xticks": xticks_dict[water_case].get(xcol, None),
            "yticks": yticks_dict[water_case].get(ycol, None),
            "ax_idx": ax_idx,
        }

    fig.set_axis(**axis_options)

    # make_legend

    # fig.save_fig(name=f"{here}/figs/{ycol}_vs_{xcol}_line_plot_SS")

    return fig



def make_legend(fig, ax_idx=0, labels=[], leg_kwargs={}):
    ax = fig.get_axis(ax_idx)
    lines = ax.lines
    # print(*lines)
    ax.legend([*lines], labels, frameon=False, **leg_kwargs)


if __name__ == "__main__":

    dm_recov = gd.get_ccro_recov_data()
    dm_ss = gd.get_ss_data()
    # dm_sec = gd.get_ccro_SEC_data()
    dm_rr = gd.get_ccro_rr_data()

    colors = itertools.cycle(line_colors)
    plot_options = {
        "marker": "o",
        "markersize": 5,
        "linestyle": "-",
        "color": next(colors),

    }
    # plot2_options = {
    #     "marker": "o",
    #     "markersize": 5,
    #     "linestyle": "-",
    #     "color": next(colors),

    # }
    xcol = "Water recovery"
    ycol = "LCOW"
    init_figure = {
        "width": 2,
        "height": 2,
        "nrows": 1,
        "ncols": 1,
    }

    fig = FigureGenerator()
    fig.init_figure(**init_figure)
    for case in cases:
        fig = plot_ss_line(dm_ss, water_case=case, fig=fig, xcol=xcol, ycol=ycol, plot_by_stage=False, plot_options=plot_options, )
        plot_options["color"] = next(colors)
        fig.show()
    make_legend(fig, labels=cases)
    fig.save_fig(name=f"{here}/figs/{ycol}_vs_{xcol}_line_plot_SS")

    # ycol = "Flux"
    # fig = plot_ss_line(dm_ss, xcol=xcol, ycol=ycol, plot_by_stage=True, plot_options=plot_options, plot2_options=plot2_options)
    # make_legend(fig, labels=["Stage 1", "Stage 2"])
    # fig.save_fig(name=f"{here}/figs/{ycol}_vs_{xcol}_line_plot_SS")
    # plot_ss_line(dm_ss, plot_by_stage=False)

    # ycols = ["LCOW", "SEC", "Flushing efficiency", "Membrane Area", "Flux", "Pump work"]
    # plot_options = {
    #     "marker": "o",
    #     "markersize": 5,
    #     "linestyle": "-",
    # }

    # xcol = "Water recovery"
    # dm = dm_rr

    # for ycol in ycols:
    #     # for dm, xcol in [(dm_recov, "Water recovery"), (dm_rr, "Recycle rate")]:
    #     colors = itertools.cycle(line_colors)
    #     for i, case in enumerate(cases):
    #         color = next(colors)
    #         plot_options["color"] = color
    #         print(ycol, case, color)
    #         if i == 0:
    #             fig = plot_ccro_line(
    #                 dm,
    #                 fig=None,
    #                 water_case=case,
    #                 xcol=xcol,
    #                 ycol=ycol,
    #                 plot_options=plot_options,
    #             )
    #         else:
    #             fig = plot_ccro_line(
    #                 dm,
    #                 fig=fig,
    #                 water_case=case,
    #                 xcol=xcol,
    #                 ycol=ycol,
    #                 plot_options=plot_options,
    #             )

            # fig.show()
        # break
