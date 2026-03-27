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
cases = {
    "Brackish water": {
        "xticks": [75, 80, 85, 90, 95],
        "yticks": [0, 0.1, 0.2, 0.3, 0.4],
    },
    "Seawater": {
        "xticks": [40, 45, 50, 55, 60],
        "yticks": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
    },
    "Produced water": {
        "xticks": [20, 30, 40, 50, 55],
        "yticks": [
            0,
            0.5,
            1.0,
            1.5,
            2.0,
        ],
    },
}

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
    "Specific area": "Spec area (m$^2$/(L/hr))",
    "Pressure": "RO Pressure (bar)",
    "Pump work": "Pump Work (kW)",
    "Pump size": "Pump Size (kW)",
    "Inlet concentration": "Inlet concentration (g/L)",
}
yticks_dict = {
    "Brackish water": {
        "LCOW": [0, 0.1, 0.2, 0.3, 0.4],
        "SEC": [0, 0.5, 1, 1.5, 2, 2.5, 3],
        "Permeate concentration": [0, 0.2, 0.4, 0.6],
        "Flux": [
            0,
            15,
            30,
            45,
            60,
            75
        ],
        "CAPEX": [0, 20, 40, 60],
        "OPEX": [0, 2, 4, 6, 8],
        "Pump work": [0, 2, 4, 6, 8, 10],
        "Pump size": [0, 2, 4, 6, 8, 10, 12],
        "Area" :[0, 50, 100, 150],
        "Pressure" :[0, 25, 50, 75, 100],
        "Inlet concentration": [0, 10, 20, 30],
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
        "Pressure": [0, 40, 80, 120],
        "Inlet concentration": [0, 15, 30, 45, 60, 75],
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
        "Pressure": [0, 50, 100, 150, 200],
        "Inlet concentration": [0, 40, 80, 120, 160],
    },
}
xticks_dict = {
    "Brackish water": [75, 80, 85, 90, 95],
    "Seawater": [40, 45, 50, 55, 60],
    "Produced water": [20, 30, 40, 50, 55],
}



def get_ccro_data():

    dm = PsDataManager()

    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )

    # COSTING
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")

    dm.register_data_key(
        "costing.total_capital_cost",
        "CAPEX",
        assign_units="kUSD",
        conversion_factor=1e-3,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "OPEX",
        assign_units="kUSD/year",
        conversion_factor=1e-3,
    )
    # SYSTEM

    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    # dm.register_data_key("recycle_flowrate", "Recycle rate", "L/s")
    dm.register_data_key(
        "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "Recycle rate",
        "L/s",
    )
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "m^3/s")
    # dm.register_data_key("avg_product_flow_rate", "Permeate flow rate", "m^3/s")
    dm.register_data_key("filtration_ramp_rate", "Ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    dm.register_data_key("overall_rejection", "Overall rejection", "%")
    dm.register_data_key("cycle_time_ratio", "Cycle time ratio", "%")
    dm.register_data_key("permeate_concentration", "Permeate concentration", "g/L")
    # dm.register_data_key(
    #     "blocks[0].process.fs.P2.control_volume.properties_out[0.0].pressure",
    #     "Recycle Pump Pressure",
    #     "bar",
    # )

    # RO
    dm.register_data_key("blocks[0].process.fs.RO.area", "Area", "m^2")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]",
        "Single Pass Recovery",
        "%",
    )

    dm.register_data_key(
        "blocks[19].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Pressure",
        "bar",
    )
    # dm.register_data_key(
    #     "blocks[19].process.fs.P1.control_volume.pressure",
    #     "Pressure",
    #     "bar",
    # )
    dm.register_data_key(
        "blocks[19].process.fs.P1.control_volume.work[0.0]",
        "Pump size",
        "kW",
    )
    dm.register_data_key(
        "blocks[19].process.fs.P1.total_power",
        "Pump work",
        "kW",
    )
    dm.register_data_key(
        f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Inlet concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[19].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Final concentration",
        "g/L",
    )
    # dm.register_data_key(
    #     f"blocks[0].process.fs.product.properties[0.0].flow_vol_phase[Liq]",
    #     "Permeate flow rate",
    #     "m^3/s",
    # )
    dm.register_data_key(
        f"total_permeate_vol",
        "Permeate flow rate",
        # "m^3/s",
    )

    dm.load_data()

    fs = dm.get_expression_keys()
    dm.register_expression(
        fs.Pump_size / fs.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Pump_work / fs.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Average_product_flow_rate / fs.Area,
        "Flux",
        "L/(m^2*hr)",
    )
    dm.register_expression(
        fs.Area/fs.Average_product_flow_rate,
        "Specific area",
        "m^2/(L/hr)",
    )
    dm.register_expression(fs.OPEX / fs.CAPEX, "OPEX/CAPEX Ratio")

    dm.evaluate_expressions()

    # lets create costing pacakage
    package_ccro = WaterTapCostingPackage(
        costing_block="costing", validation_key="costing.LCOW"
    )
    package_ccro.register_product_flow("avg_product_flow_rate")

    # Lets create our groups
    RO = PsCostingGroup("RO")
    RO.add_unit(
        "blocks[0].process.fs.RO",  # only adding "block[0] to specify wher capex is, normally acn just say "RO
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    feed_pump = PsCostingGroup("Feed pump")
    feed_pump.add_unit(
        "blocks[19].process.fs.P1",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    recycle_pump = PsCostingGroup("Recycle pump")
    recycle_pump.add_unit(
        "blocks[19].process.fs.P2",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    conduit = PsCostingGroup("Conduit")
    conduit.add_unit(
        "conduit",
        capex_keys="capital_cost",
    )

    cm_recov = PsCostingManager(
        dm, package_ccro, [RO, feed_pump, recycle_pump, conduit]
    )
    cm_recov.build()
    return dm


def panel_3_by_2_stack(init_figure={}, leg_kwargs=dict()):
    # dm_ccro = get_ccro_data()

    # dm_ss = getssdata()

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

    for i, (water_case, axis_options) in enumerate(cases.items()):

        axs = list()

        ax_idx = (i, 0)
        axs.append(bdp.fig.get_axis(ax_idx))
        axis_options["ax_idx"] = ax_idx

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
            labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": ""}
        else:
            labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": "Water recovery (%)"}
        axis_options.update(labels)

        bdp.plotbreakdown(
            xdata="Water recovery",
            ydata="costing",
            axis_options=axis_options,
            ax_idx=ax_idx,
        )

        ########################################
        ######################################## STEADY STATE
        ########################################

        ax_idx = (i, 1)
        axs.append(bdp.fig.get_axis(ax_idx))
        axis_options["ax_idx"] = ax_idx
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
        labels = {
            "ylabel": "",
        }
        axis_options.update(labels)
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

    ax1 = fig.ax[2][0]
    ax2 = fig.ax[2][1]

    # Collect handles and labels from both axes
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
    # import itertools
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
    # fig.fig.savefig("tes")
    fig.show()

    # for ax in fig.ax:
    #     # print(ax.get_ylims())
    #     for a in ax:
    #         print(a.get_ylim())

    return fig, bdp


def panel_3_by_2_line(ycol="SEC", init_figure={}, ccro_kwargs=dict(), ro_kwargs=dict(), row_leg=0):

    # dm_ccro = get_ccro_data()

    # dm_ss = getssdata()

    fig = FigureGenerator()
    fig.init_figure(**init_figure)

    for i, (water_case, axis_options) in enumerate(cases.items()):

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

    dm_ss = gd.getssdata()
    panel_3_by_2_stack(init_figure=init_figure)

    ccro_kwargs = {"marker": "o", "color": line_colors[1], "label": "CCRO"}
    ro_kwargs = {"marker": "o", "color": line_colors[4], "label": "RO"}

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
    panel_3_by_2_line(
        ycol="Permeate flow rate",
        init_figure=init_figure,
        ccro_kwargs=ccro_kwargs,
        ro_kwargs=ro_kwargs,
    )
    panel_3_by_2_line(
        ycol="Permeate volume",
        init_figure=init_figure,
        ccro_kwargs=ccro_kwargs,
        ro_kwargs=ro_kwargs,
    )

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
