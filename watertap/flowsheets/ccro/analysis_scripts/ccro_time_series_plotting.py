import pandas as pd
import numpy as np
from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator

if __name__ == "__main__":
    sweep_file = "output/ccro_recovery_sweep_analysisType_SW_recovery_sweep.h5"
    # sweep_file = "output/ccro_flush_eff_sweep_var_recovery_analysisType_SW_flushing_efficiency_sweep.h5"
    sweep_file = "output/ccro_flush_eff_sweep_var_recovery_analysisType_SW_flushing_eff_sweep_var_recovery.h5"
    dm = PsDataManager(sweep_file)
    n_filt_time_steps = 10
    n_flush_time_steps = 5
    n_time_periods = n_filt_time_steps + n_flush_time_steps

    time_periods = list(range(n_time_periods))
    t0 = 0
    tf = n_time_periods - 1
    xvar = "Water Recovery"
    xvar = "Flushing Efficiency"
    yvar = "LCOW"
    yvar = "Total Filtration Time"
    dm.register_data_key("costing.LCOW", "LCOW")
    dm.register_data_key("costing.SEC", "SEC")
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

    for t in range(n_flush_time_steps):
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
            f"blocks[{t}].process.fs.P1.control_volume.deltaP[0.0]",
            (
                t,
                "Pump 1 dP",
            ),
            "bar",
        )
        dm.register_data_key(
            f"blocks[{t}].process.fs.P2.control_volume.deltaP[0.0]",
            (
                t,
                "Pump 2 dP",
            ),
            "bar",
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
        dm.register_data_key(
            f"ramp_rate[{t}]",
            (
                t,
                "Ramp Rate",
            ),
            # conversion_factor=14.5037737730, assign_units="psi/min",
            # "psi/min",
        )

    dm.load_data()
    dm.display()

    t_sequence = list()
    y_sequence = list()
    # yvar = "Recycle Rate"
    yvar = "Pump 1 Pressure"
    # yvar = "Ramp Rate"
    # yvar = "Pump 2 dP"
    tvar = "Operation Time Points"
    recovery = 0.5
    for t in time_periods:
        t_sequence.append(dm[(("overall_recovery", recovery), (t, tvar))].data)
        y_sequence.append(dm[(("overall_recovery", recovery), (t, yvar))].data)
        # print(dm[(("overall_recovery", 0.6), (t, tvar))].data)
    #     t_sequence.append(dm[(t, tvar)].data)
    #     y_sequence.append(dm[(t, yvar)].data)

    y_sequence = [list(group) for group in zip(*y_sequence)]
    # t_sequence.append(dm[(t, tvar)].data)
    t_sequence = [list(group) for group in zip(*t_sequence)]
    print(t_sequence)

    fig = FigureGenerator()
    fig.init_figure()
    # # print(dm[xvar].data)

    for t, s, r in zip(
        t_sequence, y_sequence, dm[(("overall_recovery", recovery), xvar)].data
    ):

        _r = f"{int(r)}%"
        fig.plot_line(xdata=t, ydata=s, label=r)

    fig.set_axis_ticklabels(
        xlabel="Cycle Time (min)",
        ylabel=yvar,
        ax_idx=0,
    )
    fig.add_legend()
    fig.show()

    fig_save = sweep_file.replace(".h5", f"_{yvar}.png")
    fig.save_fig(name=fig_save)

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
