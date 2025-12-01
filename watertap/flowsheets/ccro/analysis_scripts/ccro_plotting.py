from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
from psPlotKit.data_plotter.ps_map_plotter import mapPlotter
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
            "output/ccro_brine_sweep_analysisType_ccro_sweep_brackish_recovery.h5",
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
            "filekey": "blocks[0].process.fs.costing.LCOW",
            "return_key": "LCOW",
        },
        {
            "filekey": f"blocks[0].process.fs.RO.area",
            "return_key": "RO area",
        },
        {
            "filekey": f"blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
            "return_key": "Recycle flow",
        },
    ]
    time_periods = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

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
        # import_keys.append(
        #     {
        #         "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].pressure",
        #         "return_key": (t, "Recycle TDS"),
        #     }
        # )

    data_manager.load_data(import_keys, exact_keys=True)
    data_manager.display()

    fig = figureGenerator()
    fig.init_figure()
    fig.plot_line(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data,
        data_manager["ccro_sweep_brackish_recovery/overall_recovery", "LCOW"].data,
        color="black",
    )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="LCOW ($\$$/m$^3$)",
        xticks=[50, 60, 70, 80, 90],
        yticks=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    )
    # fig.show()
    fig.save(file_name="lcow", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    data_manager[
        "ccro_sweep_brackish_recovery/overall_recovery", "Recycle flow"
    ].to_units("L/min")
    fig.plot_line(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data,
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Recycle flow"
        ].data,
        color="black",
    )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Recycle flow rate (L/min)",
        xticks=[50, 60, 70, 80, 90],
        yticks=[0, 10, 20, 30],
    )
    # fig.show()
    fig.save(file_name="Recycle", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    # data_manager["ccro_sweep_brackish_recovery/overall_recovery", "Area flow"].to_units(
    #     "l/hr"
    # )
    fig.plot_line(
        data_manager[
            "ccro_sweep_brackish_recovery/overall_recovery", "Water recovery"
        ].data,
        data_manager["ccro_sweep_brackish_recovery/overall_recovery", "RO area"].data,
        color="black",
    )
    fig.set_axis(
        xlabel="Water recovery (%)",
        ylabel="Area (m$^2$)",
        xticks=[50, 60, 70, 80, 90],
        yticks=[0, 10, 20, 30, 40],
    )
    # fig.show()
    fig.save(file_name="Area", save_location="figs")

    d = get_sequence(
        data_manager,
        "ccro_sweep_brackish_recovery/overall_recovery",
        "RO feed TDS",
        time_periods,
    )

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
    # fig.show()
    fig.save(file_name="feed tds", save_location="figs")

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
    fig.save(file_name="retentate tds", save_location="figs")

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
    fig.save(file_name="pressure", save_location="figs")
