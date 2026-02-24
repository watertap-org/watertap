from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
import numpy as np
from psPlotKit.data_manager.ps_data import psData
from watertap.flowsheets.ccro import PS_RO
import pandas as pd


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        if (dir, (t, key)) in data_manager:
            sequence.append(data_manager[dir, (t, key)].data)
    sequence = np.array(sequence)
    print(sequence.T)
    return sequence.T[0]


def get_data(df, key):
    comp_time = df["Runtime (hr)"].to_numpy() * 60

    time_frame = 25
    time_start = 1.75
    comp_filter = np.where(
        (comp_time > time_start) & (comp_time < (time_start + time_frame))
    )

    comp_filt_time = comp_time[comp_filter] - comp_time[comp_filter][0]
    print(comp_time, comp_filt_time)
    return comp_filt_time, df[key].to_numpy()[comp_filter]


def get_val_data():
    data = pd.read_csv("../validation_data/Approx_3hr_15LMH_50Recirc.csv")

    return data


if __name__ == "__main__":
    data_manager = psDataManager(
        [
            "output/validation_sweep_ccro_analysisType_val_sweep.h5",
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
                "filekey": f"blocks[{t}].process.fs.dead_volume.accumulation_time[0.0]",
                "return_key": (t, "cycle time"),
                "units": "s",
            },
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].pressure",
                "return_key": (t, "Feed pressure"),
                "units": "bar",
            },
        )
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,1.0].pressure",
                "return_key": (t, "Retentate pressure"),
                "units": "bar",
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
    # data_manager.display_loaded_contents()
    data = get_val_data()
    line_plots_options = {
        "True": {"color": "#de2d26", "label": "With Hold Up"},
        "False": {"color": "#3182bd", "label": "No Hold Up"},
    }
    d = get_sequence(
        data_manager,
        ("use_ro_with_hold_up", "True"),
        "cycle time",
        time_periods,
    )
    sim_time_points = np.array(list(range(len(d) - 1))) * d[0] / 60
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = get_sequence(
            data_manager,
            ("use_ro_with_hold_up", lp),
            "Feed pressure",
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
    t, pressure = get_data(data, "Feed Pressure (psi)")
    fig.plot_line(
        t,
        pressure * 0.06894757,
        label="Experimental Data",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (min)",
        ylabel="Feed pressure (bar)",
        xticks=[0, 5, 10, 15, 20],
        yticks=[0, 5, 10, 15, 20, 25, 30, 35, 40],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="ccro_val_feed_pressure", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        d = get_sequence(
            data_manager,
            ("use_ro_with_hold_up", lp),
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
    t, pressure = get_data(data, "Reject Concetration (g/L)")
    fig.plot_line(
        t,
        pressure,
        label="Experimental Data",
        color="black",
        marker="o",
        ls="",
        zorder=5,
    )
    fig.set_axis(
        xlabel="Time (min)",
        ylabel="Retentate TDS (g/L)",
        xticks=[0, 5, 10, 15, 20],
        yticks=[0, 5, 10, 15, 20, 25, 30, 35, 40],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="ccro_val_retentante_tds", save_location="figs")
    fig.show()
