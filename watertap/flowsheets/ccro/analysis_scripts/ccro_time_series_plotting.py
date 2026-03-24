import pandas as pd
import numpy as np
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator


def get_sequence(data_manager, dir, key, time_periods):
    sequence = []
    for t in time_periods:
        if isinstance(dir, str):
            dir = (dir,)
        if (*dir, (t, key)) in data_manager:
            sequence.append(data_manager[(*dir, (t, key))].data)
    sequence = np.array(sequence)
    return sequence.T


if __name__ == "__main__":
    dm = PsDataManager()
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        "output/ccro_recovery_sweep_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )

    n_filt_time_steps = 20
    n_flush_time_steps = 5
    n_time_periods = n_filt_time_steps + n_flush_time_steps

    time_periods = list(range(n_time_periods))
    t0 = 0
    tf = n_time_periods - 1
    xvar = "Water Recovery"
    xvar = "Flushing Efficiency"
    yvar = "LCOW"
    yvar = "Total Filtration Time"
    dm.register_data_key("total_cycle_time", "Total Cycle Time", "s")
    dm.register_data_key("total_filtration_time", "Total Filtration Time")
    dm.register_data_key("total_flush_volume", "Total Flushing Volume")
    dm.register_data_key("avg_product_flow_rate", "Average Product Flow Rate", "L/s")
    dm.register_data_key("cycle_time_ratio", "Cycle Time Ratio", "%")
    dm.register_data_key("flushing.flushing_efficiency", "Flushing Efficiency")
    dm.register_data_key("overall_recovery", "Water Recovery", "%")
    dm.register_data_key(
        "costing.aggregate_flow_electricity", "Total Power Required", "kW"
    )
    dm.register_data_key(
        "costing.aggregate_flow_costs[electricity]", "Total Electricity Cost"
    )
    dm.register_data_key("costing.total_capital_cost", "Total Capital Cost")
    dm.register_data_key("costing.total_operating_cost", "Total Operating Cost")
    dm.register_data_key("costing.capital_recovery_factor", "CRF")
    dm.register_data_key("costing.utilization_factor", "UF")
    dm.register_data_key("costing.maintenance_labor_chemical_factor", "MLCF")
    dm.register_data_key("costing.wacc", "WACC")
    dm.register_data_key("costing.wacc", "WACC")

    units = ["Conduit", "Pump", "ReverseOsmosis1DwithHoldUp"]
    unit_names = ["Conduit", "Pump", "RO"]
    for u, n in zip(units, unit_names):
        dm.register_data_key(
            f"costing.LCOW_aggregate_direct_capex[{u}]", f"{n} Direct CAPEX"
        )

    for t in time_periods:
        print(t)

        dm.register_data_key(
            f"operation_time_points[{t}]",
            (t, "Operation Time Points"),
            conversion_factor=1 / 60,
            assign_units="min",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.dead_volume.accumulation_time[0.0]",
            (
                t,
                "Accumulation Time",
            ),
            # "kW",
        )
        dm.register_data_key(
            f"flushing.flushing_time",
            (
                t,
                "Flushing Time",
            ),
            # "kW",
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
        dm.register_data_key(
            f"blocks[{t}].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
            (
                t,
                "Recycle Rate",
            ),
            "L/s",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.P1.control_volume.properties_out[0.0].pressure",
            (
                t,
                "Pump 1 Pressure",
            ),
            "bar",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.P2.control_volume.properties_out[0.0].pressure",
            (
                t,
                "Pump 2 Pressure",
            ),
            "bar",
        )

        dm.register_data_key(
            f"blocks[{t}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
            (
                t,
                "RO recovery",
            ),
            "%",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,1.0].conc_mass_phase_comp[Liq,NaCl]",
            (
                t,
                "RO outlet TDS",
            ),
            "g/L",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
            (
                t,
                "RO inlet TDS",
            ),
            "g/L",
        )
        dm.register_data_key(
            f"ramp_rate[{t}]",
            (
                t,
                "Ramp Rate",
            ),
        )

    dm.load_data()
    dm.display()
    cases = {
        "Brackish water": {
            "recoveries": [95],  # [75, 85, 95],
            "xticks": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "yticks_pressure": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            "yticks_flux": [0, 10, 20, 30, 40, 50, 60],
            "yticks_tds": [0, 20, 40, 60, 80, 100],
        },
        "Seawater": {
            "recoveries": [45, 51, 60],
            "xticks": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "yticks_pressure": [0, 20, 40, 60, 80, 100, 120],
            "yticks_flux": [0, 10, 20, 30, 40, 50],
            "yticks_tds": [0, 20, 40, 60, 80, 100, 120],
        },
        "Produced water": {
            "recoveries": [20, 40, 55],
            "xticks": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
            "yticks_pressure": [0, 25, 50, 75, 100, 125, 150, 175, 200],
            "yticks_flux": [0, 10, 20, 30, 40, 50],
            "yticks_tds": [0, 50, 100, 150],
        },
    }
    for case, opts in cases.items():
        fig = FigureGenerator()
        fig.init_figure()
        for r in opts["recoveries"]:
            time = get_sequence(dm, case, "Operation Time Points", time_periods)
            pressure = get_sequence(dm, case, "Pump 1 Pressure", time_periods)
            print(
                r,
                dm[(case, "Water Recovery")].data,
            )
            idx = np.where(abs(dm[(case, "Water Recovery")].data - float(r)) <= 1e-8)[
                0
            ][0]
            time = [0] + list(time[idx])
            pressure = [pressure[idx][-1]] + list(pressure[idx])
            fig.plot_line(
                xdata=time[:21],
                ydata=pressure[:21],
                label=f"Standard operation",  # f"Recovery of {r}%",
                marker="o",
                color="black",
                zorder=10,
            )
            fig.plot_line(
                xdata=time[20:],
                ydata=pressure[20:],
                label=f"Flushing operation",  # f"Recovery of {r}%",
                marker="o",
                color="blue",
            )
        fig.set_axis(
            xlabel="Cycle Time (min)",
            ylabel="Pressure (bar)",
            xticks=opts["xticks"],
            yticks=opts["yticks_pressure"],
        )
        fig.add_legend()
        fig.save("ccro_time_series_figs", f"{case}_pressure_time_series.png")
        fig = FigureGenerator()
        fig.init_figure()
        for i, r in enumerate(opts["recoveries"]):
            time = get_sequence(dm, case, "Operation Time Points", time_periods)
            tds_inlet = get_sequence(dm, case, "RO inlet TDS", time_periods)
            tds_outlet = get_sequence(dm, case, "RO outlet TDS", time_periods)
            print(
                r,
                dm[(case, "Water Recovery")].data,
            )
            idx = np.where(abs(dm[(case, "Water Recovery")].data - float(r)) <= 1e-8)[
                0
            ][0]
            # print(
            #     len(time[idx]),
            #     len(pressure[idx]),
            #     r,
            #     dm[(case, "Water Recovery")].data,
            #     idx,
            # )
            time = [0] + list(time[idx])
            tds_inlet = [tds_inlet[idx][-1]] + list(tds_inlet[idx])
            tds_outlet = [tds_outlet[idx][-1]] + list(tds_outlet[idx])
            # fig.plot_line(
            #     xdata=time,
            #     ydata=tds_inlet,
            #     color=i,
            #     marker="o",
            # )
            fig.plot_line(
                xdata=time[:21],
                ydata=tds_outlet[:21],
                color="black",
                marker="o",
                label="Standard operation",  # f"Recovery of {r}%",
                zorder=10,
            )
            fig.plot_line(
                xdata=time[20:],
                ydata=tds_outlet[20:],
                label=f"Flushing operation",  # f"Recovery of {r}%",
                color="blue",
                marker="o",
            )
            # fig.plot_line([], [], label=f"Recovery of {r}% ", marker="o", color=i)
        # fig.plot_line([], [], label=f"Inlet", marker="o", color="black")
        # fig.plot_line([], [], label=f"Outlet", marker="d", color="black")
        fig.set_axis(
            xlabel="Cycle Time (min)",
            ylabel="RO outlet TDS (g/L)",
            xticks=opts["xticks"],
            yticks=opts["yticks_tds"],
        )
        fig.add_legend()
        fig.save("ccro_time_series_figs", f"{case}_tds_time_series")
    # t_sequence = list()
    # y_sequence = list()
    # # yvar = "Recycle Rate"
    # yvar = "Pump 1 Pressure"
    # # yvar = "Ramp Rate"
    # # yvar = "Pump 2 dP"
    # tvar = "Operation Time Points"
    # for t in time_periods:
    #     t_sequence.append(dm[(t, tvar)].data)
    #     y_sequence.append(dm[(t, yvar)].data)
    # y_sequence = [list(group) for group in zip(*y_sequence)]
    # # t_sequence.append(dm[(t, tvar)].data)
    # t_sequence = [list(group) for group in zip(*t_sequence)]
    # fig = FigureGenerator()
    # fig.init_figure()

    # for t, s, r in zip(t_sequence, y_sequence, dm[xvar].data):

    #     _r = f"{int(r)}%"
    #     fig.plot_line(xdata=t, ydata=s, label=_r)

    # fig.set_axis_ticklabels(
    #     xlabel="Cycle Time (s)",
    #     ylabel=yvar,
    #     ax_idx=0,
    # )
    # fig.add_legend()
    # # fig.show()

    # fig_save = sweep_file.replace(".h5", f"_{yvar}.png")
    # fig.save_fig(name=fig_save)
