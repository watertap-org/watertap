from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
import numpy as np
from psPlotKit.data_manager.ps_data import psData
from watertap.flowsheets.ccro import PS_RO


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        if (dir, (t, key)) in data_manager:
            sequence.append(data_manager[dir, (t, key)].data)
    sequence = np.array(sequence)
    print(sequence.T)
    return sequence.T[0]


if __name__ == "__main__":
    val_data = PS_RO.load_validation_data_into_model(
        file_path="../validation_data/sine_700-900psi_60s_period.csv",
    )
    data_manager = psDataManager(
        [
            "output/val_sweep_analysisType_val_sweep.h5",
        ]
    )
    import_keys = [
        {
            "filekey": f"blocks[0].process.fs.RO.area",
            "return_key": "RO area",
        },
    ]
    time_periods = range(200)
    for t in time_periods:
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.time_point",
                "return_key": (t, "Simulated time point"),
                "units": "s",
            },
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.velocity[0.0,0.1]",
                "return_key": (t, "RO inlet velocity"),
            },
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.mixed_permeate[0.0].flow_vol_phase[Liq]",
                "return_key": (t, "RO mixed permeate flow"),
                "units": "L/hr",
            },
        )

        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
                "return_key": (t, "RO SP"),
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,1.0].conc_mass_phase_comp[Liq,NaCl]",
                "return_key": (t, "Retentante TDS"),
            }
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.mixed_permeate[0.0].conc_mass_phase_comp[Liq,NaCl]",
                "return_key": (t, "Permeate TDS"),
            }
        )
    data_manager.load_data(import_keys, exact_keys=True)
    data_manager.display()

    line_plots_options = {
        "True": {"color": "#de2d26", "label": "With Hold Up"},
        "False": {"color": "#3182bd", "label": "No Hold Up"},
    }
    sim_time_points = get_sequence(
        data_manager,
        ("use_hold_up", "True"),
        "Simulated time point",
        time_periods,
    )
    # sim_time_points = sim_time_points * list(range(len(sim_time_points)))

    exp_time_points = val_data["Runtime (min)"].to_numpy() * 60
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = get_sequence(
            data_manager,
            ("use_hold_up", lp),
            "RO inlet velocity",
            time_periods,
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    # fig.plot_line(
    #     exp_time_points,
    #     val_data["Permeate 1 Flowrate (L/min)"].to_numpy()
    #     / val_data["Feed Flowrate (L/min)"].to_numpy()
    #     * 100,
    #     label="Experimental Data",
    #     color="black",
    #     marker="o",
    #     ls="",
    #     zorder=5,
    # )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Cross flow velocity (%)",
        xticks=[0, 50, 100, 150, 200],
        yticks=[0, 0.1, 0.2, 0.3, 0.4],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="validation_water_flux", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = (
            get_sequence(
                data_manager,
                ("use_hold_up", lp),
                "RO SP",
                time_periods,
            )
            * 100
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    fig.plot_line(
        exp_time_points,
        val_data["Permeate 1 Flowrate (L/min)"].to_numpy()
        / val_data["Feed Flowrate (L/min)"].to_numpy()
        * 100,
        label="Experimental Data",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Water recovery (%)",
        xticks=[0, 50, 100, 150, 200],
        yticks=[0, 5, 10, 15, 20, 25, 30],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="validation_water_flux", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = (
            get_sequence(
                data_manager,
                ("use_hold_up", lp),
                "RO mixed permeate flow",
                time_periods,
            )
            / 7.2
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    fig.plot_line(
        exp_time_points,
        val_data["Water Flux (LMH)"].to_numpy(),
        label="Experimental Data",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Water flux (LMH)",
        xticks=[0, 50, 100, 150, 200],
        yticks=[0, 10, 20, 30, 40, 50],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="validation_water_flux", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = get_sequence(
            data_manager,
            ("use_hold_up", lp),
            "Retentante TDS",
            time_periods,
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    fig.plot_line(
        exp_time_points - 12,
        val_data["Reject Conductivity (mS/cm)"].to_numpy() * 0.6638 * 1.08,
        label="Experimental Data (offset -11 s)",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Retentante (g/L)",
        xticks=[0, 50, 100, 150, 200],
        yticks=[36, 38, 40, 42, 44],
    )
    fig.add_legend(location="upper right")
    # fig.show()
    fig.save(file_name="validation_retentate_tds", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = get_sequence(
            data_manager,
            ("use_hold_up", lp),
            "Permeate TDS",
            time_periods,
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    fig.plot_line(
        exp_time_points - 12,
        val_data["Permeate 1 Conductivity (uS/cm)"].to_numpy() / 1000 * 0.6638,
        label="Experimental Data",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Permeate (g/L)",
        xticks=[0, 50, 100, 150, 200],
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
    )
    fig.add_legend(location="upper right")
    # fig.show()
    fig.save(file_name="validation_permeate_tds", save_location="figs")
    fig.show()
