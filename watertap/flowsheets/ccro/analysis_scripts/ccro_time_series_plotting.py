import pandas as pd
import numpy as np
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator
if __name__ == "__main__":
    sweep_file = "/Users/ksitterl/Documents/Python/watertap/watertap/watertap/flowsheets/ccro/analysis_scripts/output/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5"
    dm = PsDataManager(sweep_file)
    n_filt_time_steps = 20
    n_flush_time_steps = 5
    n_time_periods = n_filt_time_steps + n_flush_time_steps

    time_periods = list(range(n_time_periods))
    t0 = 0
    tf = n_time_periods - 1
    xvar = "Water Recovery"
    yvar = "LCOW"
    yvar = "Total Filtration Time"
    dm.register_data_key("costing.LCOW", "LCOW")
    dm.register_data_key("costing.SEC", "SEC")
    dm.register_data_key("total_cycle_time", "Total Cycle Time", "s")
    dm.register_data_key("total_filtration_time", "Total Filtration Time")
    dm.register_data_key("total_flush_volume", "Total Flushing Volume")
    dm.register_data_key("avg_product_flow_rate", "Average Product Flow Rate", "L/s")
    dm.register_data_key("cycle_time_ratio", "Cycle Time Ratio", "%")
    dm.register_data_key("overall_recovery", "Water Recovery", "%")
    dm.register_data_key("costing.aggregate_flow_electricity", "Total Power Required", "kW")
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

    for t in range(n_time_periods):
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
            f"blocks[{t}].process.fs.RO.area",
            (
                t,
                "RO area",
            ),
            "m^2",
        )

    dm.load_data()

    filt_t_sequence = list()
    for t in time_periods:
        filt_t_sequence.append(dm[(t, "Accumulation Time")].data)
    filt_t_sequence = [list(s) for s in zip(*filt_t_sequence)]

    flushing_t_sequence = list()
    for t in time_periods:
        flushing_t_sequence.append(dm[(t, "Flushing Time")].data)
    flushing_t_sequence = [list(s) for s in zip(*flushing_t_sequence)]

    y_sequence = list()
    yvar = "Recycle Rate"
    yvar = "Pump 1 Pressure"
    for t in time_periods:
        y_sequence.append(dm[(t, yvar)].data)
    y_sequence = [list(group) for group in zip(*y_sequence)]


    fig = FigureGenerator()
    fig.init_figure()

    # cycle_time_sequence = list()
    for s, r, dt, ft in zip(
        y_sequence, dm[xvar].data, filt_t_sequence, flushing_t_sequence
    ):
        filt_dt = dt[0]
        flush_dt = ft[0] / n_flush_time_steps
        cts_filt = [filt_dt * _t for _t in range(n_filt_time_steps)]
        cts_flush = [(flush_dt * _t) + cts_filt[-1] for _t in range(n_flush_time_steps)]
        cts = cts_filt + cts_flush
        # cycle_time_sequence.append(cts)
        _r = f"{int(r)}%"
        fig.plot_line(xdata=cts, ydata=s, label=_r)

    fig.set_axis_ticklabels(
        xlabel="Cycle Time (s)",
        ylabel=yvar,
        ax_idx=0,
    )
    fig.add_legend()
    fig.show()

    fig_save = sweep_file.replace(".h5", f"_{yvar}.png")
    fig.save_fig(name=fig_save)