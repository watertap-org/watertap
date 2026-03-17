from psPlotKit.data_manager.ps_data_manager import PsDataManager

from psPlotKit.data_plotter.ps_break_down_plotter import BreakDownPlotter
from psPlotKit.data_manager.costing_packages.watertap_costing import (
    WaterTapCostingPackage,
)
from psPlotKit.data_manager.ps_costing import (
    PsCostingGroup,
    PsCostingManager,
)
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)

if __name__ == "__main__":
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

    dm.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm.register_data_key("fs.costing.LCOW", "LCOW")

    dm.load_data()
    package = WaterTapCostingPackage()
    # costing_block="fs.costing", validation_key="fs.costing.LCOW"
    # )
    package.register_product_flow()  # "fs.product.properties[0.0].flow_vol_phase[Liq]")
    # Lets create our groups
    RO = PsCostingGroup("RO")
    RO.add_unit(
        "RO",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    pumps = PsCostingGroup("Pumps & ERD")
    pumps.add_unit(
        "pump",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    pumps.add_unit(
        "ERD",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )

    cm = PsCostingManager(dm, package, [RO, pumps])
    cm.build()
    dm.display()

    cases = {
        "Brackish water": {
            "xticks": [75, 80, 85, 90, 95],
            "yticks": [0, 0.1, 0.2, 0.3, 0.4],
        },
        "Seawater": {
            "xticks": [45, 50, 55, 60],
            "yticks": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        },
        "Produced water": {
            "xticks": [20, 30, 40, 50, 55],
            "yticks": [
                0,
                0.2,
                0.4,
                0.6,
                0.8,
                1.0,
                1.2,
                1.4,
                1.6,
            ],
        },
    }
    for water_case, axis_options in cases.items():
        dm.select_data(water_case)
        wr_water = dm.get_selected_data()
        for erd in ["True", "False"]:
            wr_water.select_data(("add_erd", erd))
            wr = wr_water.get_selected_data()
            wr.display()
            cost_plotter = BreakDownPlotter(
                wr,
                save_folder="figs",
                save_name=f"steady_state_cost_breakdown_{water_case}_ERD_{erd}",
                show_fig=True,
            )
            if erd == "True":
                areas = [
                    # {"ERD": {"label": None, "color": "#f0dcd3"}},
                    {"Pumps & ERD": {"label": "Pumps & ERD", "color": "#1f78b4"}},
                    {"RO": {"label": "RO", "color": "#a6cee3"}},
                ]
            else:
                areas = [
                    # {"ERD": {"label": None, "color": "#f0dcd3"}},
                    {"Pumps & ERD": {"label": "Pumps", "color": "#1f78b4"}},
                    {"RO": {"label": "RO", "color": "#a6cee3"}},
                ]
            cost_plotter.define_area_groups(areas)
            cost_plotter.define_hatch_groups(
                {
                    "LCOW_opex": {"label": "OPEX", "hatch": "", "color": "none"},
                    "LCOW_capex": {"label": "CAPEX", "hatch": "\\\\", "color": "none"},
                }
            )
            labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": "Water recovery (%)"}
            labels.update(axis_options)
            cost_plotter.plotbreakdown(
                xdata="Water recovery",
                ydata="costing",
                axis_options=labels,
                legend_loc="upper left",
                generate_figure=True,
            )
            # cost_plotter.fig.plot_line(
            #     dm["optimal_design", "Water recovery"],
            #     dm["optimal_design", "LCOW"],
            #     label="Optimal design",
            #     color="red",
            #     linestyle="--",
            # )
            # cost_plotter.generate_figure()
