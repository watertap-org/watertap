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
    # dm.register_data_file(
    #     "output_13/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5",
    #     directory="Seawater",
    # )
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "L/hr")
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")

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
    # # dm.display()
    # # assert False
    # for water_case in water_cases:
    #     dm.select_data(water_case)
    #     wr = dm.get_selected_data()
    #     print(water_case)
    #     print("RO")
    #     wr[water_case, ("costing", "RO", "total_capital_cost")].display()

    #     wr[water_case, ("costing", "RO", "total_operating_cost")].display()
    #     wr[water_case, ("costing", "RO", "LCOW")].display()
    #     print("Feed pump")
    #     wr[water_case, ("costing", "Feed pump", "total_capital_cost")].display()

    #     wr[water_case, ("costing", "Feed pump", "total_operating_cost")].display()
    #     wr[water_case, ("costing", "Feed pump", "LCOW")].display()
    #     print("Recycle pump")
    #     wr[water_case, ("costing", "Recycle pump", "total_capital_cost")].display()

    #     wr[water_case, ("costing", "Recycle pump", "total_operating_cost")].display()
    #     wr[water_case, ("costing", "Recycle pump", "LCOW")].display()
    #     print("Conduit")
    #     wr[water_case, ("costing", "Conduit", "total_capital_cost")].display()

    #     wr[water_case, ("costing", "Conduit", "total_operating_cost")].display()
    #     wr[water_case, ("costing", "Conduit", "LCOW")].display()
    #     print("calc LCOW")
    #     wr[water_case, ("costing", "total", "total_capital_cost")].display()

    #     wr[water_case, ("costing", "total", "total_operating_cost")].display()
    #     wr[water_case, ("costing", "total", "LCOW")].display()
    #     wr[water_case, ("costing", "aggregate_flow_cost")].display()
    #     print("actual lcow")
    #     wr[water_case, "LCOW"].display()
    #     assert False
    # cost_plotter = BreakDownPlotter(wr)
    # cost_plotter.define_area_groups(
    #     {"RO": {}, "Feed pump": {}, "Recycle pump": {}, "Conduit": {}}
    # )
    # cost_plotter.plotbreakdown(
    #     xdata="Water recovery",
    #     ydata=["costing", "LCOW"],
    #     axis_options={
    #         "yticks": [0, 0.5, 1.0, 1.5, 2.0],
    #         "xticks": "auto",
    #     },
    #     generate_figure=False,
    # )
