from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
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

if __name__ == "__main__":
    dm = PsDataManager()
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "L/hr")
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")

    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.load_data()

    fs_keys = dm.get_expression_keys()
    fs_keys.print_mapping()
    # this creates a simple class that enables one to store labels and plot options, with automatic color assignment. This also should expose standard options for lines for autocomplete via .add method
    water_cases = FigureGenerator.get_plot_options_manager()
    water_cases.add("Brackish water")
    water_cases.add("Seawater")
    water_cases.add("Produced water")

    # lets create costing pacakage
    package = WaterTapCostingPackage(
        costing_block="costing", validation_key="costing.LCOW"
    )
    package.register_product_flow("avg_product_flow_rate")

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

    cm = PsCostingManager(dm, package, [RO, feed_pump, recycle_pump, conduit])
    cm.build()
    # dm.display()
    # # assert False
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
            "yticks": [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
        },
    }
    for water_case, axis_options in cases.items():
        dm.select_data(water_case)
        wr = dm.get_selected_data()
        cost_plotter = BreakDownPlotter(
            wr,
            save_folder="ccro_cost_breakdown_figs",
            save_name=f"ccro_cost_breakdown_{water_case}",
        )
        cost_plotter.define_area_groups(
            {"Feed pump": {}, "Recycle pump": {}, "RO": {}, "Conduit": {}}
        )
        labels = {"ylabel": "LCOW ($\$$/m$^3$)", "xlabel": "Water recovery (%)"}
        labels.update(axis_options)
        cost_plotter.define_hatch_groups(
            {
                "LCOW_opex": {"label": "OPEX", "hatch": "", "color": "none"},
                "LCOW_capex": {"label": "CAPEX", "hatch": "\\\\", "color": "none"},
            }
        )
        cost_plotter.plotbreakdown(
            xdata="Water recovery",
            ydata="costing",
            axis_options=labels,
            generate_figure=True,
            fig_options={"width": 2, "height": 2},
        )
