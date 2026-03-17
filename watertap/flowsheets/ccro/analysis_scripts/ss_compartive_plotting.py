from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)


def getssdata():

    dm_ss = PsDataManager()
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_sweep.h5",
        directory="Brackish water",
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_PW_sweep.h5",
        directory="Produced water",
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_SW_sweep.h5",
        directory="Seawater",
    )
    for stage in [1, 2]:
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
    dm_ss.register_data_key("fs.ERD.control_volume.work[0.0]", "ERD work", "kW")
    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")

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
    dm_ss.load_data()
    dm_ss.display()
    # dm_ss.display()
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

    return dm_ss


if __name__ == "__main__":
    dm = getssdata()
    water_cases = FigureGenerator.get_plot_options_manager()
    water_cases.add("Brackish water")
    water_cases.add("Seawater")
    water_cases.add("Produced water")
    design_cases = FigureGenerator.get_plot_options_manager()
    design_cases.add("1 stage 1 pump", color="black", marker="o", ls="-")
    design_cases.add("2 stage 1 pump", color="black", marker="s", ls="-")
    design_cases.add("2 stage 2 pump", color="black", marker="d", ls="-")
    for data in dm.data_keys:
        if data != "Water recovery":
            fig = FigureGenerator().init_figure()
            for c, opts in water_cases.items():
                for d, dopt in design_cases.items():
                    if (
                        c,
                        ("stage_sim_cases", d.replace(" ", "_")),
                        data,
                    ) in dm:
                        opts_combined = opts.copy()
                        opts_combined.ls = dopt.ls
                        opts_combined.marker = dopt.marker
                        opts_combined.label = None
                        fig.plot_line(
                            dm[
                                c,
                                ("stage_sim_cases", d.replace(" ", "_")),
                                "Water recovery",
                            ],
                            dm[c, ("stage_sim_cases", d.replace(" ", "_")), data],
                            **opts_combined,
                        )
                    else:
                        print(f"Data not found for {c}, {d}")
            for c, dopt in water_cases.items():
                fig.plot_line([], [], **dopt)
            for c, dopt in design_cases.items():
                fig.plot_line([], [], **dopt)
            fig.add_legend()
            fig.set_axis(xlabel="auto", ylabel="auto")
            fig.save_fig(f"figs/ss_{data}_vs_recovery.png")
    fig.show()
