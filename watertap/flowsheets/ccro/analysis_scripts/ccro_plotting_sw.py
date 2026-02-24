from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
import numpy as np
from psPlotKit.data_manager.ps_data import psData


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        if (dir, (t, key)) in data_manager:
            sequence.append(data_manager[dir, (t, key)].data)
    sequence = np.array(sequence)

    return sequence.T


if __name__ == "__main__":
    data_manager = psDataManager(
        [
            "output/ccro_brine_sweep_analysisType_study_SW_optimization_lcow.h5",
        ]
    )
    import_keys = [
        {
            "filekey": "overall_recovery",
            "return_key": "Water recovery",
            "units": "%",
        },
        {
            "filekey": "final_concentration",
            "return_key": "Final TDS",
            "units": "g/L",
        },
        {
            "filekey": "total_cycle_time",
            "return_key": "Total cycle time",
            "units": "min",
        },
        {
            "filekey": "costing.LCOW",
            "return_key": "LCOW",
        },
        {
            "filekey": "costing.specific_energy_consumption[0]",
            "return_key": "SEC",
        },
        {
            "filekey": f"blocks[0].process.fs.RO.area",
            "return_key": "RO area",
        },
        {
            "filekey": f"blocks[0].process.fs.RO.length",
            "return_key": "RO length",
        },
        {
            "filekey": f"blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
            "return_key": "Recycle flow",
        },
        # {
        #     "filekey": f"blocks[0].process.fs.RO.area",
        #     "return_key": "Flushing efficiency",
        # },
    ]
    time_periods = range(200)

    for t in time_periods:
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
                "return_key": (t, "RO feed TDS"),
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,1.0].conc_mass_phase_comp[Liq,NaCl]",
                "return_key": (t, "RO retentate TDS"),
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].pressure",
                "return_key": (t, "RO pressure"),
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.flushing.flushing_efficiency",
                "return_key": "Flushing efficiency",
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
                "return_key": "RO SP recovery",
            }
        )
        # import_keys.append(
        #     {
        #         "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].pressure",
        #         "return_key": (t, "Recycle TDS"),
        #     }
        # )

    data_manager.load_data(import_keys, exact_keys=True)
    data_manager.display()
    line_plots_options = {
        "True": {"color": "#de2d26", "label": "With Hold Up"},
        "False": {"color": "#3182bd", "label": "No Hold Up"},
    }
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "LCOW"].data,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="LCOW ($\$$/m$^3$)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 0.5, 1, 1.5, 2, 2.5],
    )

    fig.add_legend()

    fig.save(file_name="SW lcow", save_location="figs")
    # fig = figureGenerator()
    # fig.init_figure()
    # for lp in line_plots_options:
    #     fig.plot_line(
    #         data_manager[("use_hold_up", lp), "Water recovery"].data,
    #         data_manager[("use_hold_up", lp), "SEC"].data,
    #         color=line_plots_options[lp]["color"],
    #     )
    # fig.set_axis(
    #     xlabel="Water recovery (%)",
    #     ylabel="SEC (kWh/m$^3$)",
    #     xticks=[40,50,60,70],
    #     yticks=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    # )
    # fig.show()
    fig.save(file_name="SW lcow", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        data_manager[("use_hold_up", lp), "Recycle flow"].to_units("L/min")
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "Recycle flow"].data,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Recycle flow rate (L/min)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 50, 100, 150, 200],
    )
    fig.add_legend()
    # fig.show()
    fig.save(file_name="SW Recycle", save_location="figs")

    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        data_manager[("use_hold_up", lp), "Recycle flow"].to_units("L/min")
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "RO area"].data,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Area (m$^2$)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 50, 100, 150, 200, 250, 300, 350, 400],
    )
    # fig.show()
    fig.add_legend()
    fig.save(file_name="SW Area", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        data_manager[("use_hold_up", lp), "Recycle flow"].to_units("L/min")
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "RO length"].data,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Length (m)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 2, 4, 6, 8, 10, 12],
    )
    # fig.show()
    fig.add_legend()
    fig.save(file_name="SW Length", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    for lp in line_plots_options:
        data_manager[("use_hold_up", lp), "Recycle flow"].to_units("L/min")
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "Total cycle time"].data,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Total cycle time (min)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 10, 20, 30, 40, 50, 60],
    )
    # fig.show()
    fig.add_legend(location="upper right")
    fig.save(file_name="SW Total cycle time", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    for lp in line_plots_options:
        data_manager[("use_hold_up", lp), "Recycle flow"].to_units("L/min")
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "Flushing efficiency"].data * 100,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Flushing efficiency (%)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 25, 50, 75, 100],
    )
    # fig.show()
    fig.add_legend(location="upper right")
    fig.save(file_name="SW Flushing efficiency", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    for lp in line_plots_options:
        fig.plot_line(
            data_manager[("use_hold_up", lp), "Water recovery"].data,
            data_manager[("use_hold_up", lp), "RO SP recovery"].data * 100,
            color=line_plots_options[lp]["color"],
            label=line_plots_options[lp]["label"],
        )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="RO SP recovery (%)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 25, 50, 75, 100],
    )
    # fig.show()
    fig.add_legend(location="upper right")
    fig.save(file_name="SW SP recovery", save_location="figs")
    fig.show()

    # d = get_sequence(
    #     data_manager,
    #     "ccro_sweep_brackish_recovery/overall_recovery",
    #     "RO feed TDS",
    #     time_periods,
    # )

    # fig = figureGenerator()
    # fig.init_figure()
    # map_object = fig.gen_colormap(num_samples=14)
    # print(
    #     data_manager[
    #         "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #     ].data
    # )
    # # assert False
    # for i, tds in enumerate(
    #     data_manager[
    #         "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #     ].data
    # ):
    #     if (
    #         i
    #         != len(
    #             data_manager[
    #                 "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #             ].data
    #         )
    #         - 1
    #     ):
    #         fig.plot_line(
    #             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #             d[i],
    #             color=map_object(i / 14),
    #             label=f"{int(round(tds, 0))} %",
    #             marker="o",
    #         )
    # fig.add_legend()
    # fig.set_axis(
    #     xlabel="Time period (#)",
    #     ylabel="RO Feed TDS (g/L)",
    #     xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #     yticks=[0, 20, 40, 60, 80, 100],
    # )
    # # fig.show()
    # fig.save(file_name="SW feed tds", save_location="figs")

    # d = get_sequence(
    #     data_manager,
    #     "ccro_sweep_brackish_recovery/overall_recovery",
    #     "RO retentate TDS",
    #     time_periods,
    # )

    # fig = figureGenerator()
    # fig.init_figure()
    # map_object = fig.gen_colormap(num_samples=14)
    # for i, tds in enumerate(
    #     data_manager[
    #         "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #     ].data
    # ):
    #     if (
    #         i
    #         != len(
    #             data_manager[
    #                 "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #             ].data
    #         )
    #         - 1
    #     ):
    #         fig.plot_line(
    #             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #             d[i],
    #             color=map_object(i / 14),
    #             label=f"{int(round(tds, 0))} %",
    #             marker="o",
    #         )
    # fig.add_legend()
    # fig.set_axis(
    #     xlabel="Time period (#)",
    #     ylabel="RO Retentate TDS (g/L)",
    #     xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #     yticks=[0, 20, 40, 60, 80, 100],
    # )
    # # fig.show()
    # fig.save(file_name="SW retentate tds", save_location="figs")

    # d = get_sequence(
    #     data_manager,
    #     "ccro_sweep_brackish_recovery/overall_recovery",
    #     "RO pressure",
    #     time_periods,
    # )

    # fig = figureGenerator()
    # fig.init_figure()
    # map_object = fig.gen_colormap(num_samples=14)
    # for i, tds in enumerate(
    #     data_manager[
    #         "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #     ].data
    # ):
    #     if (
    #         i
    #         != len(
    #             data_manager[
    #                 "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
    #             ].data
    #         )
    #         - 1
    #     ):
    #         fig.plot_line(
    #             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #             d[i] / 1e5,
    #             color=map_object(i / 14),
    #             label=f"{int(round(tds, 0))} %",
    #             marker="o",
    #         )
    # fig.add_legend()
    # fig.set_axis(
    #     xlabel="Time period (#)",
    #     ylabel="RO pressure (bar)",
    #     xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #     yticks=[0, 20, 40, 60, 80, 100],
    # )
    # fig.show()
    # fig.save(file_name="SW pressure", save_location="figs")
