from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
import numpy as np
from psPlotKit.data_manager.ps_data import psData


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        print((*dir, (t, key)))
        if (*dir, (t, key)) in data_manager:
            sequence.append(data_manager[(*dir, (t, key))].data)
    sequence = np.array(sequence)

    return sequence.T


def add_legend(fig, cl, ms):
    for c in cl:
        fig.plot_line([], [], color=cl[c], ls="-", marker="s", label=f"Periods: {c}")

    for m in ms:
        fig.plot_line(
            [],
            [],
            color="black",
            ls=ms[m]["ls"],
            marker=ms[m]["marker"],
            label=ms[m]["label"],
        )


if __name__ == "__main__":
    data_manager = psDataManager(
        [
            "output/ccro_brine_sweep_analysisType_study_SW_mesh_study_optimization_lcow.h5",
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
            "units": "L/s",
        },
        # {
        #     "filekey": f"blocks[0].process.fs.RO.area",
        #     "return_key": "Flushing efficiency",
        # },
    ]
    time_periods = range(300)

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

    mesh_steps = {11: "#a6cee3", 21: "#1f78b4", 51: "#b2df8a"}  # , 101: "#33a02c"}
    data_manager.load_data(import_keys, exact_keys=True)
    data_manager.display()

    line_plots_options = {
        "True": {"ls": "--", "marker": "o", "label": "With Hold Up"},
        "False": {"ls": "-", "marker": "d", "label": "No Hold Up"},
    }
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        for ts in mesh_steps:
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[("time_steps", ts), ("use_hold_up", lp), "LCOW"].data,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    add_legend(fig, mesh_steps, line_plots_options)
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="LCOW ($\$$/m$^3$)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 0.5, 1.0, 1.5, 2.0, 2.5],
    )
    fig.add_legend()
    # fig.show()
    fig.save(file_name="SW lcow", save_location="figs")
    # assert False
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        for ts in mesh_steps:
            data_manager[
                ("time_steps", ts), ("use_hold_up", lp), "Recycle flow"
            ].to_units("L/s")
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Recycle flow"
                ].data,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Recycle flow rate (L/s)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 1, 2, 3, 4, 5],
    )
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend()
    # fig.show()
    fig.save(file_name="SW Recycle", save_location="figs")

    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        for ts in mesh_steps:
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[("time_steps", ts), ("use_hold_up", lp), "RO area"].data,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend()
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Area (m$^2$)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 100, 200, 300, 400],
    )
    # fig.show()
    fig.save(file_name="SW Area", save_location="figs")
    # assert False

    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        for ts in mesh_steps:
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[("time_steps", ts), ("use_hold_up", lp), "RO length"].data,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend()
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Length (m)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 2, 4, 6, 8, 10],
    )
    # fig.show()
    fig.save(file_name="SW Length", save_location="figs")

    # fig = figureGenerator()
    # fig.init_figure()
    # for ts in mesh_steps:
    #     # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     #     "l/hr"
    #     # )
    #     fig.plot_line(
    #         data_manager[("time_steps", ts), "Water recovery"].data,
    #         data_manager[("time_steps", ts), "RO width"].data,
    #         label=f"Time periods: {ts}",
    #     )

    # add_legend(fig, mesh_steps, line_plots_options)
    # fig.add_legend()
    # fig.set_axis(
    #     xlabel="Water recovery (%)",
    #     ylabel="width (m)",
    #     xticks=[40, 50, 60, 70],
    #     yticks=[0, 200, 400, 600, 800, 1000],
    # )
    # # fig.show()
    # fig.save(file_name="SW Width", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()

    for lp in line_plots_options:
        for ts in mesh_steps:

            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "RO SP recovery"
                ].data
                * 100,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend()
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="RO Single Pass Recovery (%)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 20, 40, 60, 80, 100],
    )
    fig.save(file_name="SW RO SP", save_location="figs")
    # fig = figureGenerator()
    # fig.init_figure()

    # # assert False
    # for lp in line_plots_options:
    #     for ts in mesh_steps:

    #         residence_time = (
    #             data_manager[("time_steps", ts), ("use_hold_up", lp), "RO length"].data
    #             / data_manager[
    #                 ("time_steps", ts), ("use_hold_up", lp), "RO inlet velocity"
    #             ].data
    #         )
    #         fig.plot_line(
    #             data_manager[
    #                 ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
    #             ].data,
    #             residence_time,
    #             color=mesh_steps[ts],
    #             ls=line_plots_options[lp]["ls"],
    #             marker=line_plots_options[lp]["marker"],
    #         )
    # add_legend(fig, mesh_steps, line_plots_options)
    # fig.add_legend()
    # fig.set_axis(
    #     xlabel="Water recovery (%)",
    #     ylabel="Residence Time (s)",
    #     xticks=[40, 50, 60, 70],
    #     yticks=[0, 20, 40, 60],
    # )
    # fig.save(file_name="SW RO residence time", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    for lp in line_plots_options:
        for ts in mesh_steps:
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Total cycle time"
                ].data,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Total cycle time (min)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 10, 20, 30, 40, 50, 60],
    )
    # fig.show()
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend(location="upper right")
    fig.save(file_name="SW Total cycle time", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    for lp in line_plots_options:
        for ts in mesh_steps:
            fig.plot_line(
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Water recovery"
                ].data,
                data_manager[
                    ("time_steps", ts), ("use_hold_up", lp), "Flushing efficiency"
                ].data
                * 100,
                color=mesh_steps[ts],
                ls=line_plots_options[lp]["ls"],
                marker=line_plots_options[lp]["marker"],
            )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Flushing efficiency (%)",
        xticks=[40, 50, 60, 70],
        yticks=[0, 25, 50, 75, 100],
    )
    # fig.show()
    add_legend(fig, mesh_steps, line_plots_options)
    fig.add_legend(location="upper right")
    fig.save(file_name="SW Flushing efficiency", save_location="figs")
    fig.show()
    assert False
    fig = figureGenerator()
    fig.init_figure()
    map_object = fig.gen_colormap(num_samples=14)
    print(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data
    )
    # assert False
    for i, tds in enumerate(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data
    ):
        if (
            i
            != len(
                data_manager[
                    "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
                ].data
            )
            - 1
        ):
            fig.plot_line(
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                d[i],
                color=map_object(i / 14),
                label=f"{int(round(tds, 0))} %",
                marker="o",
            )
    fig.add_legend()
    fig.set_axis(
        xlabel="Time period (#)",
        ylabel="RO Feed TDS (g/L)",
        xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        yticks=[0, 20, 40, 60, 80, 100],
    )
    fig.show()
    fig.save(file_name="SW feed tds", save_location="figs")

    d = get_sequence(
        data_manager,
        "ccro_sweep_brackish_recovery/overall_recovery",
        "RO retentate TDS",
        time_periods,
    )

    fig = figureGenerator()
    fig.init_figure()
    map_object = fig.gen_colormap(num_samples=14)
    for i, tds in enumerate(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data
    ):
        if (
            i
            != len(
                data_manager[
                    "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
                ].data
            )
            - 1
        ):
            fig.plot_line(
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                d[i],
                color=map_object(i / 14),
                label=f"{int(round(tds, 0))} %",
                marker="o",
            )
    fig.add_legend()
    fig.set_axis(
        xlabel="Time period (#)",
        ylabel="RO Retentate TDS (g/L)",
        xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        yticks=[0, 20, 40, 60, 80, 100],
    )
    # fig.show()
    fig.save(file_name="SW retentate tds", save_location="figs")

    d = get_sequence(
        data_manager,
        "ccro_sweep_brackish_recovery/overall_recovery",
        "RO pressure",
        time_periods,
    )

    fig = figureGenerator()
    fig.init_figure()
    map_object = fig.gen_colormap(num_samples=14)
    for i, tds in enumerate(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data
    ):
        if (
            i
            != len(
                data_manager[
                    "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
                ].data
            )
            - 1
        ):
            fig.plot_line(
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                d[i] / 1e5,
                color=map_object(i / 14),
                label=f"{int(round(tds, 0))} %",
                marker="o",
            )
    fig.add_legend()
    fig.set_axis(
        xlabel="Time period (#)",
        ylabel="RO pressure (bar)",
        xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        yticks=[0, 20, 40, 60, 80, 100],
    )
    fig.show()
    fig.save(file_name="SW pressure", save_location="figs")
