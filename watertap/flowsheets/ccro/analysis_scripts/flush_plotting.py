from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
import numpy as np
from psPlotKit.data_manager.ps_data import psData
from watertap.flowsheets.ccro import PS_RO


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        d = []
        if isinstance(dir, list):
            for _d in dir:
                d.append(_d)
            d.append((t, key))
        else:
            d = [dir, (t, key)]
        if tuple(d) in data_manager:
            sequence.append(data_manager[tuple(d)].data)
    sequence = np.array(sequence)
    print(sequence.T)
    return sequence.T[0]


if __name__ == "__main__":

    data_manager = psDataManager(
        [
            "output/flush_sweep_analysisType_val_sweep.h5",
        ]
    )
    import_keys = [
        {
            "filekey": f"blocks[0].process.fs.RO.area",
            "return_key": "RO area",
        },
    ]
    import_keys = [
        {
            "filekey": f"flushing.flushing_time",
            "return_key": "Flushing time",
        },
    ]
    time_periods = range(200)
    for t in time_periods:
        import_keys.append(
            {
                "filekey": f"blocks[{t}].process.fs.RO.feed_side.accumulation_time[0.0]",
                "return_key": (t, "accumulation_time"),
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
                "filekey": f"blocks[{t}].process.fs.P1.control_volume.properties_out[0.0].pressure",
                "return_key": (t, "RO inlet pressure"),
                "units": "bar",
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
        5: {"color": "#3f3f3f", "label": "5 steps"},
        10: {"color": "#de2d26", "label": "10 steps"},
        20: {"color": "#3182bd", "label": "20 steps"},
        50: {"color": "#48bd31", "label": "50 steps"},
        # 100: {"color": "#bd31b1", "label": "100 steps"},
    }

    # sim_time_points = sim_time_points * list(range(len(sim_time_points)))
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        sim_time_points = get_sequence(
            data_manager,
            [("add_flushing", "True"), ("time_steps", lp)],
            "accumulation_time",
            time_periods,
        )
        flush_time = data_manager[
            ("add_flushing", "True"), ("time_steps", lp), "Flushing time"
        ].data
        flush_delta = flush_time / len(sim_time_points)
        sim_time_points[:] = flush_delta
        sim_time_points = [
            sum(sim_time_points[: i + 1]) for i in range(len(sim_time_points))
        ]
        d = (
            get_sequence(
                data_manager,
                [("add_flushing", "True"), ("time_steps", lp)],
                "RO mixed permeate flow",
                time_periods,
            )
            / 100
        )
        fig.plot_line(
            sim_time_points,
            d,
            label=line_plots_options[lp]["label"],
            color=line_plots_options[lp]["color"],
            marker="o",
            zorder=10,
        )
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Water flux (LMH)",
        xticks=[0, 10, 20, 30, 40, 50, 60],
        yticks=[0, 10, 20, 30, 40, 50],
    )
    fig.add_legend(location="upper right")
    fig.save(file_name="flush_water_flux", save_location="figs")

    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        sim_time_points = get_sequence(
            data_manager,
            [("add_flushing", "True"), ("time_steps", lp)],
            "accumulation_time",
            time_periods,
        )
        flush_time = data_manager[
            ("add_flushing", "True"), ("time_steps", lp), "Flushing time"
        ].data
        flush_delta = flush_time / len(sim_time_points)
        sim_time_points[:] = flush_delta
        sim_time_points = [
            sum(sim_time_points[: i + 1]) for i in range(len(sim_time_points))
        ]
        d = get_sequence(
            data_manager,
            [("add_flushing", "True"), ("time_steps", lp)],
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
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="Retentante (g/L)",
        xticks=[0, 10, 20, 30, 40, 50, 60],
        yticks=[0, 5, 10, 15, 20, 25, 30, 35, 40],
    )
    fig.add_legend(location="upper right")
    # fig.show()
    fig.save(file_name="flushing_retentate_tds", save_location="figs")
    fig = figureGenerator()
    fig.init_figure()
    for lp in line_plots_options:
        sim_time_points = get_sequence(
            data_manager,
            [("add_flushing", "True"), ("time_steps", lp)],
            "accumulation_time",
            time_periods,
        )
        flush_time = data_manager[
            ("add_flushing", "True"), ("time_steps", lp), "Flushing time"
        ].data
        flush_delta = flush_time / len(sim_time_points)
        sim_time_points[:] = flush_delta
        sim_time_points = [
            sum(sim_time_points[: i + 1]) for i in range(len(sim_time_points))
        ]
        d = get_sequence(
            data_manager,
            [("add_flushing", "True"), ("time_steps", lp)],
            "RO inlet pressure",
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
    fig.set_axis(
        xlabel="Time (s)",
        ylabel="RO Inlet Pressure (Bar)",
        xticks=[0, 10, 20, 30, 40, 50, 60],
        yticks=[0, 10, 20, 30, 40, 50, 60],
    )
    fig.add_legend(location="upper right")
    # fig.show()
    fig.save(file_name="flushing_outlet_pressure", save_location="figs")
    fig.show()
