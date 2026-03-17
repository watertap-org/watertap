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
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane area", "m^2")
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

    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "L/hr")
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    for t in range(25):

        dm.register_data_key(
            f"operation_time_points[{t}]",
            (t, "Operation Time Points"),
            conversion_factor=1 / 60,
            assign_units="min",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.P1.control_volume.work[0.0]",
            (
                t,
                "Pump 1 Power",
            ),
            "kW",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.P2.control_volume.work[0.0]",
            (
                t,
                "Pump 2 Power",
            ),
            "kW",
        )
    dm.load_data()

    fs_keys = dm.get_expression_keys()
    fs_keys.print_mapping()
    for t in range(25):
        if t == 0:
            dm.register_expression(
                fs_keys[(0, "Operation Time Points")],
                (t, "fractional_time"),
                assign_units="dimensionless",
            )
        else:
            dm.register_expression(
                (
                    fs_keys[(t, "Operation Time Points")]
                    - fs_keys[(t - 1, "Operation Time Points")]
                )
                / fs_keys["Total cycle time"],
                (t, "fractional_time"),
                assign_units="dimensionless",
            )
    dm.evaluate_expressions()
    p1_power = 0
    p2_power = 0
    for t in range(25):
        p1_power += fs_keys[(t, "Pump 1 Power")] * fs_keys[(t, "fractional_time")]
        p2_power += fs_keys[(t, "Pump 2 Power")] * fs_keys[(t, "fractional_time")]
    dm.register_expression(p1_power, "Average Pump 1 Power", assign_units="kW")
    dm.register_expression(p2_power, "Average Pump 2 Power", assign_units="kW")
    dm.evaluate_expressions()
    dm.display()
    # assert False
    # this creates a simple class that enables one to store labels and plot options, with automatic color assignment. This also should expose standard options for lines for autocomplete via .add method
    water_cases = FigureGenerator.get_plot_options_manager()
    water_cases.add("Brackish water")
    water_cases.add("Seawater")
    water_cases.add("Produced water")

    erd_included = False
    other_d_plotted = False
    for d in dm.data_keys:
        if (
            d != "Water recovery"
            and "flux" not in str(d)
            and isinstance(d, str)
            and d != "Average Pump 2 Power"
        ):
            erd_included = False
            other_d_plotted = False
            fig = FigureGenerator().init_figure()
            for c, opts in water_cases.items():
                print(d)
                if d in ["Average Pump 1 Power", "Average Pump 2 Power"]:
                    fig.plot_line(
                        dm[c, "Water recovery"], dm[c, "Average Pump 1 Power"], **opts
                    )
                    optc = opts.copy()
                    optc.ls = "--"
                    optc.label = None
                    optc.marker = "d"
                    fig.plot_line(
                        dm[c, "Water recovery"], dm[c, "Average Pump 2 Power"], **optc
                    )
                    p1 = dm_ss["optimal_design", c, ("RO 1", "pump work")]
                    if ("optimal_design", c, ("RO 2", "RO area")) in dm_ss:
                        p1.data = (
                            p1.data
                            + dm_ss["optimal_design", c, ("RO 2", "pump work")].data
                        )
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        p1,
                        color=opts.color,
                        marker="v",
                        ls="-",
                    )
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        dm_ss["optimal_design", c, "ERD work"],
                        color=opts.color,
                        marker="*",
                        ls="--",
                    )

                    erd_included = True
                    other_d_plotted = True
                # elif d in ["Average product flow rate"]:
                #     ccro_flux = (
                #         dm[c, "Average product flow rate"].data
                #         / dm[c, "Membrane area"].data
                #     )
                #     dm.add_data(c, "CCRO flux", ccro_flux, units="m^2/L/hr")
                #     fig.plot_line(dm[c, "Water recovery"], dm[c, "CCRO flux"], **opts)

                #     area_data = dm_ss["optimal_design", c, ("RO 1", "RO area")]
                #     if ("optimal_design", c, ("RO 2", "RO area")) in dm_ss:
                #         area_data_2 = dm_ss["optimal_design", c, ("RO 2", "RO area")]
                #         area_data_2.data[area_data_2.data == 0] = area_data.data[
                #             area_data_2.data == 0
                #         ]
                #         area_data = area_data_2.data + area_data.data
                #     ro_product_flow = dm_ss[
                #         "optimal_design", c, "Average product flow rate"
                #     ].data
                #     ro_flux = ro_product_flow / area_data.data
                #     dm.add_data(c, "RO flux", ro_flux, units="m^2/L/hr")
                #     fig.plot_line(
                #         dm[c, "Water recovery"],
                #         dm[c, "RO flux"],
                #         color=opts.color,
                #         marker="d",
                #         ls="--",
                #     )
                #     other_d_plotted = True
                else:
                    fig.plot_line(dm[c, "Water recovery"], dm[c, d], **opts)
                if d == "LCOW" or d == "SEC":
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        dm_ss["optimal_design", c, d],
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True
                if d == "RO feed pressure":
                    pressure_data = dm_ss[
                        "optimal_design", c, ("RO 1", "feed pressure")
                    ]
                    if ("optimal_design", c, ("RO 2", "feed pressure")) in dm_ss:
                        pressure_data_2 = dm_ss[
                            "optimal_design", c, ("RO 2", "feed pressure")
                        ]
                        pressure_data_2.data[pressure_data_2.data == 0] = (
                            pressure_data.data[pressure_data_2.data == 0]
                        )
                        pressure_data = pressure_data_2
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        pressure_data,
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True
                if d == "Membrane area":
                    area_data = dm_ss["optimal_design", c, ("RO 1", "RO area")]
                    if ("optimal_design", c, ("RO 2", "RO area")) in dm_ss:
                        area_data_2 = dm_ss["optimal_design", c, ("RO 2", "RO area")]
                        area_data_2.data[area_data_2.data == 0] = area_data.data[
                            area_data_2.data == 0
                        ]
                        area_data = area_data_2.data + area_data.data
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        area_data,
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True
                if d == "Average product flow rate":
                    fig.plot_line(
                        dm_ss["optimal_design", c, "Water recovery"],
                        dm_ss["optimal_design", c, "Average product flow rate"],
                        color=opts.color,
                        marker="d",
                        ls="--",
                    )
                    other_d_plotted = True

            if other_d_plotted:
                if erd_included:
                    fig.plot_line(
                        [],
                        [],
                        label="CCRO pump 1",
                        color="k",
                        marker="o",
                        ls="--",
                    )
                    fig.plot_line(
                        [],
                        [],
                        label="CCRO pump 2",
                        color="k",
                        marker="d",
                        ls="--",
                    )
                    fig.plot_line(
                        [],
                        [],
                        label="Steady state RO",
                        color="k",
                        marker="s",
                        ls="--",
                    )

                    fig.plot_line(
                        [],
                        [],
                        label="Steady state ERD",
                        color="k",
                        marker="*",
                        ls="--",
                    )
                else:
                    fig.plot_line(
                        [],
                        [],
                        label="CCRO",
                        color="k",
                        marker="o",
                        ls="-",
                    )
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
    # fig.show()
