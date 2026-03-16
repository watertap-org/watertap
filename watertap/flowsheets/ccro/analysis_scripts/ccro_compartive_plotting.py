from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)


if __name__ == "__main__":
    dm = PsDataManager()
    dm.register_data_file(
        "output_13/ccro_recovery_sweep_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        "output_13/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        "output_13/ccro_recovery_sweep_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )

    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")
    dm.register_data_key("blocks[0].process.fs.RO.area", "RO area", "m^2")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]", "RO recovery", "%"
    )
    dm.register_data_key(
        "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "Recycle rate",
        "L/s",
    )
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    dm.load_data()
    dm.display()
    # this creates a simple class that enables one to store labels and plot options, with automatic color assignment. This also should expose standard options for lines for autocomplete via .add method
    water_cases = FigureGenerator.get_plot_options_manager()
    water_cases.add("Brackish water")
    water_cases.add("Seawater")
    water_cases.add("Produced water")
    for d in dm.data_keys:
        if d != "Water recovery":
            fig = FigureGenerator().init_figure()
            for c, opts in water_cases.items():
                fig.plot_line(dm[c, "Water recovery"], dm[c, d], **opts)
            fig.add_legend()
            fig.set_axis(xlabel="auto", ylabel="auto")
    fig.show()
