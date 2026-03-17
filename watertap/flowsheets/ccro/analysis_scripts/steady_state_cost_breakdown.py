from psPlotKit.data_manager.ps_data_manager import PsDataManager

from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter

if __name__ == "__main__":
    water_case = "BW"
    dm = PsDataManager()
    dm.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_PW_sweep.h5",
        directory="Produced water",
    )
    dm.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_SW_sweep.h5",
        directory="Seawater",
    )
    # this is confusing for users
    # the groups specifies desired unit operation, that includes single or multiple untis, flows and so forth
    # the units can be specified to be a singel string, list or a dict, where each unit is married to specific block
    device_groups = {
        "ERD": {"units": "ERD"},
        "Stage 1 RO": {
            "units": {
                "RO": "stage[1]",
                "pump": "stage[1]",
            },
        },
        "Stage 2 RO": {
            "units": {
                "RO": "stage[2]",
                "pump": "stage[2]",
            },
        },
    }

    dm.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm.register_data_key("fs.costing.LCOW", "LCOW")

    dm.load_data()
    dm.get_costing(
        device_groups,
        default_flow="fs.product.properties[0.0].flow_vol_phase[Liq]",
    )
    dm.display()
    dm.reduce_data(
        stack_keys="stage_sim_cases",
        data_key="LCOW",
        reduction_type="min",
        directory="optimal_design",
    )
    dm.display()
    dm.select_data("optimal_design")
    wr = dm.get_selected_data()
    cost_plotter = BreakDownPlotter(
        wr,
        save_folder="figs",
        save_name=f"cost_breakdown_{water_case}",
        show_fig=True,
    )
    cost_plotter.define_area_groups(
        [
            {"ERD": {"label": None, "color": "#f0dcd3"}},
            {"Stage 1 RO": {"label": None, "color": "#a6cee3"}},
            {"Stage 2 RO": {"label": None, "color": "#1f78b4"}},
        ]
    )
    cost_plotter.define_hatch_groups(
        {"TOTAL": {}}
    )  # {"CAPEX": {}, "OPEX": {"hatch": "//"}})
    cost_plotter.plotbreakdown(
        xdata="Water recovery",
        ydata=["cost_breakdown", "levelized"],
        axis_options={
            "yticks": [0, 0.5, 1, 1.5, 2.0],
            "xticks": [50, 60, 70, 80, 90],
        },
        legend_loc="upper left",
        generate_figure=False,
    )
    cost_plotter.fig.plot_line(
        dm["optimal_design", "Water recovery"],
        dm["optimal_design", "LCOW"],
        label="Optimal design",
        color="red",
        linestyle="--",
    )
    cost_plotter.generate_figure()
