import pprint
import os
import numpy as np
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

here = os.path.dirname(os.path.abspath(__file__))
par_dir = os.path.dirname(here)

from matplotlib.patches import Patch

def getssdata():

    dm_ss = PsDataManager()
    dm_ss.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_BW_sweep.h5",
        directory="Brackish water",
    )
    dm_ss.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_PW_sweep.h5",
        directory="Produced water",
    )
    dm_ss.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_SW_sweep.h5",
        directory="Seawater",
    )
    for stage in [1, 2]:
        dm_ss.register_data_key(
            f"fs.stage[{stage}].RO.flux_mass_phase_comp_avg[0.0,Liq,H2O]",
            (f"RO {stage}", "RO flux"),
            assign_units="L/(m^2*hr)",
            conversion_factor=3600,
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.properties_out[0.0].pressure",
            (f"RO {stage}", "feed pressure"),
            "bar",
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].RO.area",
            (f"RO {stage}", "RO area"),
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.work[0.0]",
            (f"RO {stage}", "pump work"),
            "kW",
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.work[0.0]",
            (f"RO {stage}", "pump work"),
            "kW",
        )
    # dm_ss.register_data_key("fs.ERD.control_volume.work[0.0]", "ERD work", "kW")
    # dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")

    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm_ss.register_data_key("fs.costing.LCOW", "LCOW")
    dm_ss.register_data_key("fs.costing.SEC", "SEC")

    dm_ss.register_data_key(
        "fs.product.properties[0.0].flow_vol_phase[Liq]",
        "Average product flow rate",
        "L/hr",
    )
    dm_ss.register_data_key(
        "fs.product.properties[0.0].flow_vol_phase[Liq]",
        "Average product flow rate",
        "L/hr",
    )

    dm_ss.register_data_key("fs.costing.total_capital_cost", "capex")
    dm_ss.register_data_key("fs.costing.total_operating_cost", "opex")
    dm_ss.load_data()
    ek = dm_ss.get_expression_keys()
    dm_ss.register_expression(
        ek.RO_1_pump_work / ek.Average_product_flow_rate,
        "RO 1 specific pump work",
        "kW/(m^3/s)",
    )
    dm_ss.register_expression(
        ek.RO_2_pump_work / ek.Average_product_flow_rate,
        "RO 2 specific pump work",
        "kW/(m^3/s)",
    )
    dm_ss.register_expression(ek.opex / ek.capex, "opex_capex_ratio")
    dm_ss.evaluate_expressions()
    dm_ss.display()
    # dm_ss.reduce_data(
    #     stack_keys=("Produced water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Produced water", "optimal_design"),
    # )

    # dm_ss.display()
    # dm_ss.reduce_data(
    #     stack_keys=("Brackish water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Brackish water", "optimal_design"),
    # )
    # dm_ss.reduce_data(
    #     stack_keys="stage_sim_cases",
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory="optimal_design",
    # )

    # dm_ss.display()

    package_ss = WaterTapCostingPackage(
        # costing_block="fs.costing", validation_key="fs.costing.LCOW"
    )
    # package_ss.add_flow_cost("electricity", "electricity_cost")
    package_ss.register_product_flow()

    # RO1 = PsCostingGroup("RO1")
    # RO1.add_unit(
    #     "RO",
    #     capex_keys="capital_cost",
    #     fixed_opex_keys="fixed_operating_cost",
    # )
    # # RO2 = PsCostingGroup("RO2")
    # # RO2.add_unit(
    # #     "RO",
    # #     capex_keys="capital_cost",
    # #     fixed_opex_keys="fixed_operating_cost",
    # # )
    # ERD = PsCostingGroup("ERD")
    # ERD.add_unit(
    #     "ERD",
    #     capex_keys="capital_cost",
    #     # fixed_opex_keys="fixed_operating_cost",
    #     flow_keys={"electricity": ["control_volume.work"]},
    # )
    # pump = PsCostingGroup("pump")
    # pump.add_unit(
    #     "fs.stage[1].pump",
    #     capex_keys="capital_cost",
    #     # fixed_opex_keys="fixed_operating_cost",
    #     flow_keys={"electricity": ["control_volume.work"]},
    # )
    RO1 = PsCostingGroup("RO1")
    RO1.add_unit(
        "stage[1].RO",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    RO2 = PsCostingGroup("RO2")
    RO2.add_unit(
        "stage[2].RO",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    pumps = PsCostingGroup("Pumps")
    pumps.add_unit(
        "pump",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    # pump1 = PsCostingGroup("Pump 1")
    # pump1.add_unit(
    #     "stage[1].pump",
    #     capex_keys="capital_cost",
    #     flow_keys={"electricity": "control_volume.work"},
    # )
    # pump2 = PsCostingGroup("Pump 2")
    # pump2.add_unit(
    #     "stage[2].pump",
    #     capex_keys="capital_cost",
    #     flow_keys={"electricity": "control_volume.work"},
    # )
    # ERD = PsCostingGroup("ERD")
    # ERD.add_unit(
    pumps.add_unit(
        "ERD",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    cm_ss = PsCostingManager(
        dm_ss,
        package_ss,
        # [RO1, RO2, pump1, pump2, ERD],
        [RO1, RO2, pumps],
    )
    cm_ss.build()

    return dm_ss


def get_recovery_data():

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
        "capex",
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "opex",
        assign_units="MUSD/year",
        conversion_factor=1e-6,
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
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "L/hr")
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
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
    dm.register_data_key("blocks[0].process.fs.RO.area", "RO area", "m^2")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]", "RO recovery", "%"
    )
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane area", "m^2")
    # dm.register_data_key(
    #     "blocks[19].process.fs.P1.control_volume.properties_out[0.0].pressure",
    #     "RO feed pressure",
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

    dm.load_data()

    fs = dm.get_expression_keys()
    dm.register_expression(
        fs.Pump_size / fs.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Pump_work / fs.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Average_product_flow_rate / fs.RO_area,
        "RO flux",
        "L/(m^2*hr)",
    )
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


def panel_3_by_2(init_figure={}, leg_kwargs=dict()):
    cases = {
        "Brackish water": {
            "xticks": [75, 80, 85, 90, 95],
            # "yticks": [0, 0.1, 0.2, 0.3, 0.4],
            # "yticks_lcow": [0, 0.1, 0.2, 0.3, 0.4],
            # "yticks_sec": [0, 0.5, 1, 1.5, 2, 2.5, 3],
            # "yticks_pressure": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            # "yticks_flux": [0, 10, 20, 30, 40, 50, 60],
        },
        "Seawater": {
            "xticks": [45, 50, 55, 60],
            # "yticks": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            # "yticks_lcow": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            # "yticks_sec": [0, 2, 4, 6, 8],
            # "yticks_pressure": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            # "yticks_flux": [0, 10, 20, 30, 40, 50],
        },
        "Produced water": {
            "xticks": [20, 30, 40, 50, 55],
            # "yticks": [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            # "yticks_lcow": [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            # "yticks_sec": [0, 5, 10, 15, 20, 25],
            # "yticks_pressure": [0, 25, 50, 75, 100, 125, 150, 175, 200],
            # "yticks_flux": [0, 10, 20, 30, 40, 50],
        },
    }

    # water_cases = FigureGenerator.get_plot_options_manager()
    # water_cases.add("Brackish water")
    # water_cases.add("Seawater")
    # water_cases.add("Produced water")

    ########################################
    ######################################## CCRO
    ########################################
    dm_recov = get_recovery_data()

    ########################################
    ######################################## STEADY STATE
    ########################################

    dm_ss = getssdata()

    ########################################
    ########################################
    ########################################

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
        "Pumps",
        "RO1",
        "RO2",
    ]
    color_dict = dict(zip(units, line_colors))

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
            "LCOW_opex": {"label": "OPEX", "hatch": "", "color": "white", "alpha": 0.85},
            "LCOW_capex": {"label": "CAPEX", "hatch": "\\\\", "color": "white"},
        }
    )

    for i, (water_case, axis_options) in enumerate(cases.items()):

        ########################################
        ######################################## CCRO
        ########################################
        axs = list()

        ax_idx = (i, 0)
        axs.append(bdp.fig.get_axis(ax_idx))
        axis_options["ax_idx"] = ax_idx
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
            labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": ""}
        else:
            labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": "Water recovery (%)"}
        axis_options.update(labels)
        # labels.update(init_figure)
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
        
        ro_args = [("add_erd", "True"), ("stage_sim_cases", "2_stage_2_pump")]
        dm_ss.select_data(
            (water_case, *ro_args),
        )

        wr = dm_ss.get_selected_data()
        bdp.set_selected_data(wr)

        bdp.define_area_groups(
            {
                "Pumps": {},
                "RO1": {},
                "RO2": {},
            }
        )
        labels = {"ylabel": "",}
        axis_options.update(labels)

        bdp.plotbreakdown(
            xdata="Water recovery",
            ydata="costing",
            axis_options=axis_options,
            ax_idx=ax_idx,
        )

        ys = [ax.get_ylim() for ax in axs]
        y_min = min([y[0] for y in ys])
        y_max = max([y[1] for y in ys]) * 1.1

        for ax in axs:
            ax.set_ylim(y_min, y_max)

    fig.fig.tight_layout()

    ax1 = fig.ax[2][0]
    ax2 = fig.ax[2][1]

    from collections import OrderedDict

    # Collect handles and labels from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    
    combined = OrderedDict(zip(labels1 + labels2, handles1 + handles2))
    combined = {k: v for k, v in combined.items() if k not in ["CAPEX", "OPEX"]} 
    print( combined)

    fig.fig.legend(
        combined.values(),
        combined.keys(),
        loc="upper center",
        ncol=4, 
        fontsize=8,
        bbox_to_anchor=(0.54, 1.08) # position at top center
    )
    fig.save_fig(name="test_tile_plot.png")
    # fig.show()


    return fig, bdp


if __name__ == "__main__":

    nrows = 3
    ncols = 2
    w = 2
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

    panel_3_by_2(init_figure=init_figure)
