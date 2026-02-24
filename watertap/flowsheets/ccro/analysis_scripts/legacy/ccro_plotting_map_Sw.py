from psPlotKit.data_manager.ps_data_manager import psDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator
from psPlotKit.data_plotter.ps_map_plotter import MapPlotter
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
            "output/ccro_map_sweep_analysisType_study_SW_map.h5",
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
                "units": "%",
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
    for lp in line_plots_options:
        mp = MapPlotter(
            data_manager, save_folder="map_figs", save_name=f"SW lcow map {lp}"
        )
        mp.plot_map(
            [("use_hold_up", lp)],
            "Water recovery",
            "Recycle flow",
            "LCOW",
            zlevels=[1.8, 2, 2.2, 2.4, 2.6],
            axis_options={
                "xlabel": "Water recovery (%)",
                "ylabel": "Recycle flow rate (L/s)",
                "zlabel": "LCOW ($\$$/m$^3$)",
                "xticklabels": [40, 50, 60, 70],
                "yticklabels": [1, 5, 10, 15, 20],
                "zticks": [1.8, 2, 2.2, 2.4, 2.6],
            },
        )
        mp = MapPlotter(
            data_manager, save_folder="map_figs", save_name=f"SW area map {lp}"
        )
        mp.plot_map(
            [("use_hold_up", lp)],
            "Water recovery",
            "Recycle flow",
            "RO area",
            zlevels=[0, 100, 200, 300, 400, 500, 600, 700],
            axis_options={
                "xlabel": "Water recovery (%)",
                "ylabel": "Recycle flow rate (L/s)",
                "zlabel": "Area (m$^2$)",
                "xticklabels": [40, 50, 60, 70],
                "yticklabels": [1, 5, 10, 15, 20],
                "zticks": [0, 100, 200, 300, 400, 500, 600, 700],
            },
        )
        mp = MapPlotter(
            data_manager, save_folder="map_figs", save_name=f"SW cycle time map {lp}"
        )
        mp.plot_map(
            [("use_hold_up", lp)],
            "Water recovery",
            "Recycle flow",
            "Total cycle time",
            zlevels=[0, 30, 60, 90, 120, 150, 180],
            axis_options={
                "xlabel": "Water recovery (%)",
                "ylabel": "Recycle flow rate (L/s)",
                "zlabel": "Total cycle time (min)",
                "xticklabels": [40, 50, 60, 70],
                "yticklabels": [1, 5, 10, 15, 20],
                "zticks": [0, 30, 60, 90, 120, 150, 180],
            },
        )
        mp = MapPlotter(
            data_manager,
            save_folder="map_figs",
            save_name=f"SW flushing efficiency map {lp}",
        )
        mp.plot_map(
            [("use_hold_up", lp)],
            "Water recovery",
            "Recycle flow",
            "Flushing efficiency",
            zlevels=[0, 20, 40, 60, 80, 100],
            axis_options={
                "xlabel": "Water recovery (%)",
                "ylabel": "Recycle flow rate (L/s)",
                "zlabel": "Flushing efficiency (%)",
                "xticklabels": [40, 50, 60, 70],
                "yticklabels": [1, 5, 10, 15, 20],
                "zticks": [0, 20, 40, 60, 80, 100],
            },
        )
