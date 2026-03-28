import itertools
import pprint
import os
import numpy as np
from collections import OrderedDict
from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)
from psPlotKit.data_manager.costing_packages.watertap_costing import (
    WaterTapCostingPackage,
)
from psPlotKit.data_manager.ps_costing import (
    PsCostingGroup,
    PsCostingManager,
)
from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter
from watertap.flowsheets.ccro.analysis_scripts.q_report import get_data as gd

here = os.path.dirname(os.path.abspath(__file__))
par_dir = os.path.dirname(here)


line_colors = [
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
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
# cases = {
#     "Brackish water": {
#         "xticks": [75, 80, 85, 90, 95],
#         "yticks": [0, 0.1, 0.2, 0.3, 0.4],
#     },
#     "Seawater": {
#         "xticks": [40, 45, 50, 55, 60],
#         "yticks": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
#     },
#     "Produced water": {
#         "xticks": [20, 30, 40, 50, 55],
#         "yticks": [
#             0,
#             0.5,
#             1.0,
#             1.5,
#             2.0,
#         ],
#     },
# }
cases = ["Brackish water", "Seawater", "Produced water"]
ylabel_dict = {
    "LCOW": "LCOW ($\$$/m$^3$)",
    "SEC": "SEC (kWh/m$^3$)",
    "Permeate concentration": "Perm conc (g/L)",
    "Flux": "Flux (LMH)",
    # "CAPEX": "CAPEX (\\$10$^3$)",
    # "OPEX": "OPEX (\\$10$^3$ yr$^{-1}$)",
    "CAPEX": "CAPEX (kUSD)",
    "OPEX": "OPEX (kUSD/yr)",
    "Area": "RO area (m$^2$)",
    "Membrane Area": "RO area (m$^2$)",
    "Specific area": "Spec area (m$^2$/(L/hr))",
    "Pressure": "RO Pressure (bar)",
    "Pump work": "Pump Work (kW)",
    "Pump size": "Pump Size (kW)",
    "Inlet concentration": "Inlet concentration (g/L)",
        ('costing', 'total', 'LCOW_opex'): "LCOW OPEX ($\$$/m$^3$)",
        ('costing', 'total', 'LCOW_capex'): "LCOW CAPEX ($\$$/m$^3$)",
}
yticks_dict = {
    "Brackish water": {
        "LCOW": [0, 0.1, 0.2, 0.3, 0.4],
        "SEC": [0, 0.5, 1, 1.5, 2, 2.5, 3],
        "Permeate concentration": [0, 0.2, 0.4, 0.6],
        "Flux": [0, 15, 30, 45, 60, 75],
        "CAPEX": [0, 20, 40, 60],
        "OPEX": [0, 2, 4, 6, 8],
        "Pump work": [0, 2, 4, 6, 8, 10],
        "Pump size": [0, 2, 4, 6, 8, 10, 12],
        "Area": [0, 50, 100, 150],
        "Membrane Area": [0, 50, 100, 150],
        "Pressure": [0, 25, 50, 75, 100],
        "Inlet concentration": [0, 10, 20, 30],
        ('costing', 'total', 'LCOW_opex'): [0, 0.25, 0.5, 0.75, 1],
        ('costing', 'total', 'LCOW_capex'): [0, 0.25, 0.5, 0.75, 1],
    },
    "Seawater": {
        "LCOW": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
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
        "Inlet concentration": [0, 15, 30, 45, 60, 75],
        ('costing', 'total', 'LCOW_opex'): [0, 0.25, 0.5, 0.75, 1],
        ('costing', 'total', 'LCOW_capex'): [0, 0.25, 0.5, 0.75, 1],
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
        "Membrane Area": [0, 40, 80, 120, 160],
        "Pressure": [0, 50, 100, 150, 200],
        "Inlet concentration": [0, 40, 80, 120, 160],
        ('costing', 'total', 'LCOW_opex'): [0, 0.25, 0.5, 0.75, 1],
        ('costing', 'total', 'LCOW_capex'): [0, 0.25, 0.5, 0.75, 1],
    },
}
xticks_dict = {
    "Brackish water": [75, 80, 85, 90, 95],
    "Seawater": [40, 45, 50, 55, 60],
    "Produced water": [20, 30, 40, 50, 55],
}


def panel_3_by_2_stack(init_figure={}, leg_kwargs=dict()):
    # dm_ccro = gd.get_ccro_data()

    # dm_ss = gd.get_ss_data()

    fig = FigureGenerator()
    fig.init_figure(**init_figure)

    bdp = BreakDownPlotter(
        save_folder=f"{here}/figs",
        save_name=f"test",
        fig=fig,
        color_dict=color_dict,
    )

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

        dm_ccro.select_data(water_case)
        wr = dm_ccro.get_selected_data()
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
            axis_options["ylabel"] = ylabel_dict["LCOW"]
            axis_options["yticks"] = yticks_dict[water_case]["LCOW"]
            axis_options["xlabel"] = ""
        else:
            axis_options["ylabel"] = ylabel_dict["LCOW"]
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
        # axis_options["ax_idx"] = ax_idx
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
    # letters = itertools.cycle(["A", "B", "C", "D", "E", "F"])
    # for ax in fig.ax:
    #     for a in ax:
    #         a.text(
    #                 0.25,
    #                 0.95,
    #                 next(letters),
    #                 fontsize=10,
    #                 ha="center",
    #                 va="center",
    #                 transform=a.transAxes,
    #             )
    ycol = "LCOW"
    fig.save_fig(name=f"{here}/figs/{ycol}_tile_plot_stack")

    fig.show()

    return fig, bdp


def panel_3_by_2_line(
    ycol="SEC", init_figure={}, ccro_kwargs=dict(), ro_kwargs=dict(), row_leg=0
):

    # dm_ccro = get_ccro_data()

    # dm_ss = getssdata()

    fig = FigureGenerator()
    fig.init_figure(**init_figure)

    for i, water_case in enumerate(cases):

        axs = list()

        ax_idx = (i, 0)
        axs.append(fig.get_axis(ax_idx))

        dm_ccro.select_data(water_case)

        if i != len(cases) - 1:
            labels = {"ylabel": ylabel_dict.get(ycol, ""), "xlabel": ""}
        else:
            labels = {
                "ylabel": ylabel_dict.get(ycol, ""),
                "xlabel": "Water recovery (%)",
            }

        fig.plot_line(
            dm_ccro[water_case, "Water recovery"],
            dm_ccro[water_case, ycol],
            ax_idx=ax_idx,
            **ccro_kwargs,
        )

        fig.set_axis(
            xlabel=labels["xlabel"],
            ylabel=labels["ylabel"],
            xticks=xticks_dict[water_case],
            yticks=yticks_dict[water_case].get(ycol, None),
            ax_idx=ax_idx,
        )

        ax_idx = (i, 1)
        axs.append(fig.get_axis(ax_idx))

        if water_case == "Brackish water":
            ro_config = "1_stage_1_pump"
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

        try:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ("RO 1", ycol)],
                ax_idx=ax_idx,
                **ro_kwargs,
            )
        except:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ycol],
                ax_idx=ax_idx,
                **ro_kwargs,
            )
        try:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ("RO 2", ycol)],
                ax_idx=ax_idx,
                marker=ro_kwargs["marker"],
                color=line_colors[3],
                label="Stage 2",
            )
        except:
            pass

        fig.set_axis(
            xlabel=labels["xlabel"],
            ylabel="",
            xticks=xticks_dict[water_case],
            yticks=yticks_dict[water_case].get(ycol, None),
            ax_idx=ax_idx,
        )

    ax1 = fig.ax[0][0]
    ax2 = fig.ax[0][1]
    ax1.legend()
    ax2.legend()

    letters = itertools.cycle(["A", "B", "C", "D", "E", "F"])
    for ax in fig.ax:
        for a in ax:
            a.grid(visible=True)
            # a.text(
            #         0.2,
            #         0.9,
            #         next(letters),
            #         fontsize=10,
            #         ha="center",
            #         va="center",
            #         transform=a.transAxes,
            #     )
    fig.fig.tight_layout()
    fig.save_fig(name=f"{here}/figs/{ycol}_tile_plot_line")

    fig.show()

    return fig

def plot_single_compare(ycol="LCOW", init_figure={}):


    fig = FigureGenerator()
    fig.init_figure(**init_figure)

    for i, water_case in enumerate(cases):

        axs = list()

        # ax_idx = (i, 0)
        # axs.append(fig.get_axis(ax_idx))

        dm_ccro.select_data(water_case)

        if i != len(cases) - 1:
            labels = {"ylabel": ylabel_dict.get(ycol, ""), "xlabel": ""}
        else:
            labels = {
                "ylabel": ylabel_dict.get(ycol, ""),
                "xlabel": "Water recovery (%)",
            }

        fig.plot_line(
            dm_ccro[water_case, "Water recovery"],
            dm_ccro[water_case, ycol],
            # ax_idx=ax_idx,
            **ccro_kwargs,
        )

        fig.set_axis(
            xlabel=labels["xlabel"],
            ylabel=labels["ylabel"],
            xticks=xticks_dict[water_case],
            yticks=yticks_dict[water_case].get(ycol, None),
            # ax_idx=ax_idx,
        )

        # ax_idx = (i, 1)
        # axs.append(fig.get_axis(ax_idx))

        if water_case == "Brackish water":
            ro_config = "1_stage_1_pump"
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

        try:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ("RO 1", ycol)],
                # ax_idx=ax_idx,
                **ro_kwargs,
            )
        except:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ycol],
                # ax_idx=ax_idx,
                **ro_kwargs,
            )
        try:
            fig.plot_line(
                dm_ss[water_case, *ro_args, "Water recovery"],
                dm_ss[water_case, *ro_args, ("RO 2", ycol)],
                # ax_idx=ax_idx,
                marker=ro_kwargs["marker"],
                color=line_colors[3],
                label="Stage 2",
            )
        except:
            pass

        fig.set_axis(
            xlabel=labels["xlabel"],
            ylabel="",
            xticks=xticks_dict[water_case],
            yticks=yticks_dict[water_case].get(ycol, None),
            # ax_idx=ax_idx,
        )

    # ax1 = fig.ax[0][0]
    # ax2 = fig.ax[0][1]
    # ax1.legend()
    # ax2.legend()

if __name__ == "__main__":

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

    dm_ccro = gd.get_ccro_data()

    dm_ss = gd.get_ss_data()
    # panel_3_by_2_stack(init_figure=init_figure)

    ccro_kwargs = {"marker": "o", "color": line_colors[1], "label": "CCRO"}
    ro_kwargs = {"marker": "o", "color": line_colors[4], "label": "RO"}

    plot_single_compare(ycol="LCOW", init_figure={}, ccro_kwargs=ccro_kwargs, ro_kwargs=ro_kwargs)

    # panel_3_by_2_line(
    #     ycol="SEC",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="LCOW",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="CAPEX",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="OPEX",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="Permeate concentration",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="Permeate flow rate",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )
    # panel_3_by_2_line(
    #     ycol="Permeate volume",
    #     init_figure=init_figure,
    #     ccro_kwargs=ccro_kwargs,
    #     ro_kwargs=ro_kwargs,
    # )

#     ro_kwargs = {"marker": "o", "color": line_colors[4], "label": "Stage 1"}
#     panel_3_by_2_line(
#         ycol="Inlet concentration",
#         init_figure=init_figure,
#         ccro_kwargs=ccro_kwargs,
#         ro_kwargs=ro_kwargs,
#         row_leg=2
#     )
#     panel_3_by_2_line(
#         ycol="Flux",
#         init_figure=init_figure,
#         ccro_kwargs=ccro_kwargs,
#         ro_kwargs=ro_kwargs,
#         row_leg=2
#     )
#     panel_3_by_2_line(
#         ycol="Area",
#         init_figure=init_figure,
#         ccro_kwargs=ccro_kwargs,
#         ro_kwargs=ro_kwargs,
#         row_leg=2
#     )
#     panel_3_by_2_line(
#         ycol="Pump work", init_figure=init_figure, ccro_kwargs=ccro_kwargs, ro_kwargs=ro_kwargs, row_leg=2
# )
#     panel_3_by_2_line(
#         ycol="Pump size", init_figure=init_figure, ccro_kwargs=ccro_kwargs, ro_kwargs=ro_kwargs, row_leg=2
# )
#     panel_3_by_2_line(
#         ycol="Pressure", init_figure=init_figure, ccro_kwargs=ccro_kwargs, ro_kwargs=ro_kwargs, row_leg=2
# )
