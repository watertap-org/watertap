from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import FigureGenerator


def get_sequence(data_manager, dir, key, time_periods, new_dir=None, index=0):
    sequence = []

    for t in time_periods:
        if dir == "":  # top-level key
            if (t, key) in data_manager:
                sequence.append(data_manager[t, key].data[index])
        else:
            if (dir, (t, key)) in data_manager:
                sequence.append(data_manager[dir, (t, key)].data[index])
    if len(sequence) == 0:
        raise ValueError(f"No data found for key {key} in directory {dir}")
    sequence = np.array(sequence)
    if new_dir is not None:
        data_manager.add_data(new_dir, "time_periods", time_periods)
        data_manager.add_data(new_dir, key, sequence.T)
    return sequence


if __name__ == "__main__":
    # Load data
    dm = PsDataManager(
        "output/ccro_flow_sweep_analysisType_study_BGW_mesh_study_optimization_lcow.h5"
    )
    for i in range(40):
        dm.register_data_key(
            f"blocks[{i}].process.fs.RO.recovery_vol_phase[0.0,Liq]",
            (
                i,
                "RO recovery",
            ),
            "%",
        )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_out[0.0].flow_vol_phase[Liq]",
    #     "Dead volume outflow volume",
    #     "L/min",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_in[0.0].flow_vol_phase[Liq]",
    #     "Dead volume inflow volume",
    #     "L/min",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_out[0.0].flow_mass_phase_comp[Liq,NaCl]",
    #     ("Dead mass outflow mass", "NaCl"),
    #     "kg/s",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_in[0.0].flow_mass_phase_comp[Liq,NaCl]",
    #     ("Dead mass inflow mass", "NaCl"),
    #     "kg/s",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_out[0.0].flow_mass_phase_comp[Liq,H2O]",
    #     ("Dead mass outflow mass", "H2O"),
    #     "kg/s",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_in[0.0].flow_mass_phase_comp[Liq,H2O]",
    #     ("Dead mass inflow mass", "H2O"),
    #     "kg/s",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.delta_state.conc_mass_phase_comp[Liq,NaCl]",
    #     "Dead volume delta concentration",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.delta_state.mass_phase_comp[0.0,Liq,NaCl]",
    #     ("Dead volume delta mass phase", "NaCl"),
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.delta_state.mass_phase_comp[0.0,Liq,H2O]",
    #     ("Dead volume delta mass phase", "H2O"),
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.mass_phase_comp[0.0,Liq,NaCl]",
    #     ("Dead volume mass phase", "NaCl"),
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.mass_phase_comp[0.0,Liq,H2O]",
    #     ("Dead volume mass phase", "H2O"),
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_out[0.0].conc_mass_phase_comp[Liq,NaCl]",
    #     "Dead volume outlet concentration",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_in[0.0].conc_mass_phase_comp[Liq,NaCl]",
    #     "Dead volume inlet concentration",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.volume[0.0,Liq]",
    #     "Dead volume volume",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.delta_state.volume[0.0,Liq]",
    #     "Dead volume delta volume",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.accumulation_time[0.0]",
    #     "Dead volume accumulation time",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_out[0.0].dens_mass_phase[Liq]",
    #     "Dead volume outlet density",
    #     directory=i,
    # )
    # dm.register_data_key(
    #     f"blocks[{i}].process.fs.dead_volume.dead_volume.properties_in[0.0].dens_mass_phase[Liq]",
    #     "Dead volume inlet density",
    #     directory=i,
    # )

    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.load_data()
    dm.display()
    get_sequence(
        dm,
        dir=("time_steps", 20),
        key="RO recovery",
        time_periods=list(range(20)),
        new_dir="overall_recovery_sequence",
        index=16,
    )
    dm.display()
    dm["overall_recovery_sequence", "RO recovery"].display()
    # # # Plot recovery vs. time for different flushing efficiencies
    # ek = dm.get_expression_keys()
    # for e in ek:
    #     if e != "Water recovery":
    #         get_sequence(
    #             dm,
    #             dir="",
    #             key=e,
    #             time_periods=list(range(24)),
    #             new_dir="overall_recovery_sequence",
    #             index=16,
    #         )
    # dm.select_data("overall_recovery_sequence")
    # dm.display()
    # dms = dm.get_selected_data()

    # dms.display()
    # dks = dms.get_expression_keys()
    # dks.print_mapping()
    # dms.register_expression(
    #     dks.Dead_volume_mass_phase_H2O + dks.Dead_volume_mass_phase_NaCl,
    #     "Dead volume mass",
    #     assign_units="kg",
    # )
    # dms.register_expression(
    #     dks.Dead_volume_delta_mass_phase_H2O + dks.Dead_volume_delta_mass_phase_NaCl,
    #     "Dead volume delta mass",
    #     assign_units="kg",
    # )
    # dms.register_expression(
    #     dks.Dead_mass_outflow_mass_H2O + dks.Dead_mass_outflow_mass_NaCl,
    #     "Dead mass outflow",
    #     assign_units="kg/s",
    # )
    # dms.register_expression(
    #     dks.Dead_mass_inflow_mass_H2O + dks.Dead_mass_inflow_mass_NaCl,
    #     "Dead mass inflow",
    #     assign_units="kg/s",
    # )
    # dms.register_expression(
    #     (dks.Dead_mass_inflow_mass_NaCl - dks.Dead_mass_outflow_mass_NaCl)
    #     * dks.Dead_volume_accumulation_time
    #     + dks.Dead_volume_delta_mass_phase_NaCl,
    #     "Dead mass phase calc NaCl",
    #     assign_units="kg/s",
    # )
    # dms.register_expression(
    #     (-dks.Dead_mass_inflow_mass_H2O + dks.Dead_mass_outflow_mass_H2O)
    #     * dks.Dead_volume_accumulation_time
    #     + dks.Dead_volume_delta_mass_phase_H2O,
    #     "Dead mass phase calc inv H2O",
    #     assign_units="kg/s",
    # )
    # print(dms._registered_key_import_status)
    # dms.evaluate_expressions()
    # dks = dms.get_expression_keys()
    # # dms.register_expression(
    # #     (dks.Dead_mass_phase_calc_H2O + dks.Dead_mass_phase_calc_NaCl),
    # #     "Calc mass",
    # #     assign_units="kg",
    # # )
    # dms.register_expression(
    #     (dks.Dead_volume_mass / dks.Dead_volume_outlet_density),
    #     "Calc outlet vol",
    #     assign_units="L",
    # )
    # dms.register_expression(
    #     (dks.Dead_volume_mass / dks.Dead_volume_inlet_density),
    #     "Calc inlet vol",
    #     assign_units="L",
    # )
    # dms.register_expression(
    #     (dks.Dead_mass_inflow - dks.Dead_mass_outflow),
    #     "mass delta",
    #     assign_units="L",
    # )
    # dms.register_expression(
    #     (dks.Dead_mass_inflow_mass_NaCl - dks.Dead_mass_outflow_mass_NaCl)
    #     * dks.Dead_volume_accumulation_time,
    #     "accumulated nacl",
    #     assign_units="kg",
    # )
    # dms.register_expression(
    #     (dks.Dead_mass_inflow_mass_H2O - dks.Dead_mass_outflow_mass_H2O)
    #     * dks.Dead_volume_accumulation_time,
    #     "accumulated h2o",
    #     assign_units="kg",
    # )
    # dms.register_expression(
    #     (dks.Dead_mass_inflow - dks.Dead_mass_outflow)
    #     * dks.Dead_volume_accumulation_time,
    #     "overall accumulated mass",
    # )

    # dms.register_expression(
    #     (dks.Dead_volume_delta_mass - dks.Dead_volume_mass),
    #     "Change in dead volume mass",
    #     assign_units="kg",
    # )
    # dms.evaluate_expressions()

    # dks = dms.get_expression_keys()
    # dms.register_expression(
    #     (dks.accumulated_h2o - dks.accumulated_h2o) * dks.Dead_volume_accumulation_time,
    #     "accumulated mass",
    #     assign_units="kg",
    # )
    # dms.register_expression(
    #     (dks.Calc_inlet_vol - dks.Calc_outlet_vol)
    #     / dks.Dead_volume_accumulation_time
    #     * dks.Dead_volume_outlet_density,
    #     "mass outlet delta",
    #     assign_units="kg/s",
    # )
    # dms.register_expression(
    #     (dks.Calc_inlet_vol - dks.Calc_outlet_vol)
    #     / dks.Dead_volume_accumulation_time
    #     * dks.Dead_volume_inlet_density,
    #     "mass inlet delta",
    #     assign_units="kg/s",
    # )
    # dms.display()
    # dks2 = dms.get_expression_keys()
    # dks2.print_mapping()
    # print(dms["overall_recovery_sequence", "time_periods"].data.size)
    # print(dms["overall_recovery_sequence", e].data.size)

    # dms.evaluate_expressions()
    # print(dms["overall_recovery_sequence", "Dead volume accumulation time"].data)
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume accumulation time"].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Accumulation time (seconds)",
    # )
    # fig.save_fig(name="figs/acc time.png")
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "RO recovery"].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="RO single pass recovery",
    # )
    # fig.save_fig(name="figs/RO_recovery.png")
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead mass outflow"].data,
    #     label="outflow mass",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead mass inflow"].data,
    #     label="inflow mass",
    # )

    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Dead volume mass flows",
    # )
    # fig.add_legend()

    # fig.save_fig(name="figs/mass flows into dead volume.png")
    # fig = FigureGenerator()
    # fig.init_figure()
    # # fig.plot_line(
    # #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    # #     ydata=dms["overall_recovery_sequence", "accumulated mass"].data,
    # #     label="sum((in_j - out_j) * accumulation time)",
    # # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Change in dead volume mass"].data,
    #     label="Change in dead volume mass (delta_mass - mass)",
    # )
    # # fig.plot_line(
    # #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    # #     ydata=dms["overall_recovery_sequence", "overall accumulated mass"].data,
    # #     label="Overall accumulated mass ((sum(in_j) - sum(out_j)) * accumulation time)",
    # # )

    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Dead volume accumulated mass",
    # )

    # fig.save_fig(name="figs/acc dead volume.png")
    # fig.show()
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Calc mass"].data,
    #     label="Calculated mass",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume mass"].data,
    #     label="Dead volume mass",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume delta mass"].data,
    #     label="Dead volume delta mass",
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Dead volume mass",
    # )
    # fig.add_legend()
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume inlet density"].data,
    #     label="inlet density",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume outlet density"].data,
    #     label="outlet density",
    # )

    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Dead volume density",
    # )
    # fig.add_legend()

    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Calc inlet vol"].data,
    #     label="inlet volume",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Calc outlet vol"].data,
    #     label="outlet volume",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume volume"].data,
    #     label="actual volume",
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Calculated volume",
    # )
    # fig.add_legend()

    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "mass delta"].data,
    #     label="mass delta",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "mass outlet delta"].data,
    #     label="mass outlet delta",
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "mass inlet delta"].data,
    #     label="mass inlet delta",
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="Calculated volume",
    # )
    # fig.add_legend()
    # fig.show()
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead mass inflow"].data,
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead mass outflow"].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="dead mass flow",
    # )
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume accumulation time"].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="dead accumulation time",
    # )
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume delta volume"].data,
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", "Dead volume volume"].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="dead volume",
    # )
    # fig = FigureGenerator()
    # fig.init_figure()
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", ("Dead mass inflow mass", "NaCl")].data,
    # )
    # fig.plot_line(
    #     xdata=dms["overall_recovery_sequence", "time_periods"].data,
    #     ydata=dms["overall_recovery_sequence", ("Dead mass inflow mass", "H2O")].data,
    # )
    # fig.set_axis(
    #     xlabel=dms["overall_recovery_sequence", "time_periods"].data_label,
    #     ylabel="dead mass flow",
    # )
# fig.show()
# lp = line_plotter()
# lp.plot(
#     x=dm.get_column("time_hr"),
#     y=dm.get_column("recovery"),
#     group_by=dm.get_column("flushing_efficiency"),
#     xlabel="Time (hr)",
#     ylabel="Recovery",
#     title="Recovery vs. Time for Different Flushing Efficiencies",
#     legend_title="Flushing Efficiency",
# )
