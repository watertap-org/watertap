from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter
import numpy as np

if __name__ == "__main__":

    ##################################################################################
    # Data processing workflow follows plotting or getting cost breakdowns:
    # 1) Load data from h5 file using PsDataManager
    # 2) Createa  dictionary of cost groups to use for combing costs in breakdown,
    #    this is required, each group should have general name, and unit keys, if they are indexed all indexed units will be combined
    # 3) Register data keys to get general LCOW and water recovery which we sweeped over for plotting and data reduction

    # 4) Use breakdown plotter to plot cost breakdowns
    ###################################################################################

    costing_data = PsDataManager(
        [
            "output/lsrro_stage_sweep_analysisType_case_sweep.h5",
        ]
    )

    device_groups = {
        "Pumps and ERDs": {
            "units": ["PrimaryPumps", "BoosterPumps", "EnergyRecoveryDevices"],
        },
        "Membranes": {
            "units": ["ROUnits"],
        },
    }

    costing_data.register_data_key("fs.water_recovery", "Water recovery", "%")
    costing_data.register_data_key("fs.costing.LCOW", "LCOW")
    costing_data.get_costing(
        device_groups,
        default_flow="fs.product.properties[0.0].flow_vol_phase[Liq]",
    )
    costing_data.display()
    costing_data.reduce_data(
        stack_keys="number_of_stages",
        data_key="LCOW",
        reduction_type="min",
    )
    cost_breakdowns = {"CAPEX": {}, "OPEX": {"hatch": "//"}}
    markers_num_stages = {
        2: "s",
        3: "^",
        4: "D",
        5: "v",
        6: "P",
        7: "*",
        8: "X",
        9: "<",
    }

    cases = {
        "case_a": {"xticks": [30, 40, 50, 60, 70, 80]},
        "case_b": {"xticks": [30, 35, 40, 45, 50]},
    }
    for case in cases:
        costing_data.select_data(("sim_cases", case), True)
        wr = costing_data.get_selected_data()

        wr.select_data("stacked_data", True)
        wr = wr.get_selected_data()
        wr.display()
        cost_plotter = BreakDownPlotter(
            wr,
            save_name="Cost breakdown for {}".format(case),
            save_folder="figures",
        )
        cost_plotter.define_area_groups(
            [
                {"Pumps and ERDs": {"label": None, "color": "#d9f0d3"}},
                {"Membranes": {"label": None, "color": "#a6cee3"}},
            ]
        )
        cost_plotter.define_hatch_groups(cost_breakdowns)

        cost_plotter.plotbreakdown(
            xdata="Water recovery",
            ydata=["cost_breakdown", "levelized"],
            axis_options={
                "yticks": [0, 2, 4, 6, 8, 10, 12],
                "xticks": cases[case]["xticks"],
            },
            legend_loc="upper left",
            generate_figure=False,
        )
        wrs = costing_data[("stacked_data", ("sim_cases", case), "Water recovery")].data
        lcow = costing_data[("stacked_data", ("sim_cases", case), "LCOW")].data
        stages = costing_data[
            ("stacked_data", ("sim_cases", case), "number_of_stages")
        ].data
        plotted_labels = []
        for i, wr in enumerate(wrs):
            if stages[i] == stages[i]:  # only plot optimal stages
                if stages[i] not in plotted_labels:
                    label = f"{int(stages[i])} stages"
                    plotted_labels.append(stages[i])
                else:
                    label = ""
                cost_plotter.fig.plot_line(
                    [wrs[i]],
                    [lcow[i]],
                    marker=markers_num_stages[stages[i]],
                    label=label,
                    color="black",
                    ls="",
                    markersize=5,
                )
        cost_plotter.fig.add_legend()
        cost_plotter.generate_figure()
