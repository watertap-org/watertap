from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)


def getssdata():

    dm_ss = PsDataManager()
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_1_stage_1_pump_recovery_sweep.h5",
        directory=("Brackish water", "system_design", "1 stage 1 pump"),
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_2_stage_1_pump_recovery_sweep.h5",
        directory=("Brackish water", "system_design", "2 stage 1 pump"),
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_2_stage_2_pump_recovery_sweep.h5",
        directory=("Brackish water", "system_design", "2 stage 2 pump"),
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_SW_2_stage_2_pump_recovery_sweep.h5",
        directory=("Seawater", "optimal_design"),
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_PW_2_stage_1_pump_recovery_sweep.h5",
        directory=("Produced water", "system_design", "2 stage 1 pump"),
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_PW_2_stage_2_pump_recovery_sweep.h5",
        directory=("Produced water", "system_design", "2 stage 2 pump"),
    )
    for stage in [1, 2]:
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.properties_out[0.0].pressure",
            (f"RO {stage}", "feed pressure"),
            "bar",
        )
    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")

    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm_ss.register_data_key("fs.costing.LCOW", "LCOW")
    dm_ss.register_data_key("fs.costing.SEC", "SEC")
    dm_ss.load_data()
    dm_ss.display()
    dm_ss.reduce_data(
        stack_keys=("Produced water", "system_design"),
        data_key="LCOW",
        reduction_type="min",
        directory=("Produced water", "optimal_design"),
    )

    dm_ss.display()
    dm_ss.reduce_data(
        stack_keys=("Brackish water", "system_design"),
        data_key="LCOW",
        reduction_type="min",
        directory=("Brackish water", "optimal_design"),
    )
    # dm_ss.reduce_data(
    #     stack_keys=("Seawater", "system_design"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Seawater", "optimal_design"),
    # )

    dm_ss.display()
    return dm_ss


if __name__ == "__main__":
    dm_ss = getssdata()
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
    dm.register_data_key(
        "blocks[19].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "RO feed pressure",
        "bar",
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
                if d == "LCOW" or d == "SEC":
                    fig.plot_line(
                        dm_ss[(c, "optimal_design"), "Water recovery"],
                        dm_ss[(c, "optimal_design"), d],
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True
                if d == "RO feed pressure":
                    pressure_data = dm_ss[
                        (c, "optimal_design"), ("RO 1", "feed pressure")
                    ]
                    if ((c, "optimal_design"), ("RO 2", "feed pressure")) in dm_ss:
                        pressure_data_2 = dm_ss[
                            (c, "optimal_design"), ("RO 2", "feed pressure")
                        ]
                        pressure_data_2.data[pressure_data_2.data == 0] = (
                            pressure_data.data[pressure_data_2.data == 0]
                        )
                        pressure_data = pressure_data_2
                    fig.plot_line(
                        dm_ss[(c, "optimal_design"), "Water recovery"],
                        pressure_data,
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True
            if other_d_plotted:
                fig.plot_line([], [], label="CCRO", color="k", marker="o", ls="-")
                fig.plot_line(
                    [],
                    [],
                    label="Steady state RO",
                    color="k",
                    marker="d",
                    ls="--",
                )
            fig.add_legend()
            fig.set_axis(xlabel="auto", ylabel="auto")
            fig.save_fig(f"figs/{d}_vs_recovery.png")
    fig.show()
