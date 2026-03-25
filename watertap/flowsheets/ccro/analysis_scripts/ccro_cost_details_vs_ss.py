from psPlotKit.data_manager.ps_data_manager import PsDataManager
import numpy as np
from psPlotKit.data_plotter.ps_line_plotter import LinePlotter
from psPlotKit.data_plotter.fig_generator import (
    FigureGenerator,
)


def getssdata():

    dm_ss = PsDataManager()
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_BW_sweep.h5",
        directory="Brackish water",
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_PW_sweep.h5",
        directory="Produced water",
    )
    dm_ss.register_data_file(
        f"output_steady_state\RO_w_ERD_analysisType_SW_sweep.h5",
        directory="Seawater",
    )
    for stage in [1, 2]:
        dm_ss.register_data_key(
            f"fs.stage[{stage}].RO.flux_mass_phase_comp_avg[0.0,Liq,H2O]",
            (f"RO {stage}", "RO flux"),
            assign_units="L/(m^2*hr)",
            conversion_factor=3600,
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.properties_out[0.0].pressure",
            (f"RO {stage}", "feed pressure"),
            "bar",
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].RO.area",
            (f"RO {stage}", "RO area"),
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.work[0.0]",
            (f"RO {stage}", "pump work"),
            "kW",
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.work[0.0]",
            (f"RO {stage}", "pump work"),
            "kW",
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].pump.costing.capital_cost",
            (f"Pump {stage}", "capex"),
            assign_units="MUSD",
            conversion_factor=1e-6,
        )
        dm_ss.register_data_key(
            f"fs.stage[{stage}].RO.costing.capital_cost",
            (f"RO {stage}", "capex"),
            assign_units="MUSD",
            conversion_factor=1e-6,
        )
    dm_ss.register_data_key("fs.ERD.control_volume.work[0.0]", "ERD work", "kW")
    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")

    dm_ss.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm_ss.register_data_key("fs.costing.LCOW", "LCOW")
    dm_ss.register_data_key("fs.costing.SEC", "SEC")

    dm_ss.register_data_key(
        "fs.product.properties[0.0].flow_vol_phase[Liq]",
        "Average product flow rate",
        "L/hr",
    )
    dm_ss.register_data_key(
        "fs.product.properties[0.0].flow_vol_phase[Liq]",
        "Average product flow rate",
        "L/hr",
    )

    dm_ss.register_data_key(
        "fs.costing.total_capital_cost",
        "capex",
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm_ss.register_data_key(
        "fs.costing.total_operating_cost",
        "opex",
        assign_units="MUSD/year",
        conversion_factor=1e-6,
    )
    dm_ss.load_data()
    ek = dm_ss.get_expression_keys()
    dm_ss.register_expression(
        ek.RO_1_pump_work / ek.Average_product_flow_rate,
        "RO 1 specific pump work",
        "kW/(m^3/s)",
    )
    dm_ss.register_expression(
        ek.RO_2_pump_work / ek.Average_product_flow_rate,
        "RO 2 specific pump work",
        "kW/(m^3/s)",
    )
    dm_ss.register_expression(ek.opex / ek.capex, "opex_capex_ratio")
    dm_ss.evaluate_expressions()
    # dm_ss.display()
    # dm_ss.reduce_data(
    #     stack_keys=("Produced water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Produced water", "optimal_design"),
    # )

    # dm_ss.display()
    # dm_ss.reduce_data(
    #     stack_keys=("Brackish water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Brackish water", "optimal_design"),
    # )
    # dm_ss.reduce_data(
    #     stack_keys="stage_sim_cases",
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory="optimal_design",
    # )

    dm_ss.display()
    return dm_ss


if __name__ == "__main__":
    dm_ss = getssdata()
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

    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key(
        f"blocks[19].process.fs.P1.costing.capital_cost",
        (f"Feed pump", "capex"),
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        f"blocks[19].process.fs.P2.costing.capital_cost",
        (f"Recycle pump", "capex"),
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        f"blocks[19].process.fs.P2.costing.capital_cost",
        (f"Recycle pump", "capex"),
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        f"blocks[0].process.fs.RO.costing.capital_cost",
        (f"RO {stage}", "capex"),
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        "costing.total_capital_cost",
        "capex",
        assign_units="MUSD",
        conversion_factor=1e-6,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "opex",
        assign_units="MUSD/year",
        conversion_factor=1e-6,
    )
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")
    dm.register_data_key("blocks[0].process.fs.RO.area", "RO area", "m^2")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]", "RO recovery", "%"
    )
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane area", "m^2")
    dm.register_data_key(
        "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "Recycle rate",
        "L/s",
    )
    dm.register_data_key(
        "blocks[19].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "RO feed pressure",
        "bar",
    )
    dm.register_data_key(
        "blocks[19].process.fs.P1.control_volume.work[0.0]",
        "Pump size",
        "kW",
    )
    dm.register_data_key(
        "blocks[19].process.fs.P1.total_power",
        "Pump work",
        "kW",
    )
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "L/hr")
    dm.register_data_key("filtration_ramp_rate", "Filtration ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")

    dm.load_data()
    fs = dm.get_expression_keys()
    dm.register_expression(
        fs.Pump_size / fs.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Pump_work / fs.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        fs.Average_product_flow_rate / fs.RO_area,
        "RO flux",
        "L/(m^2*hr)",
    )
    dm.evaluate_expressions()
    cases = {
        "Brackish water": {
            "xticks": [75, 80, 85, 90, 95],
            "yticks_lcow": [0, 0.1, 0.2, 0.3, 0.4],
            "yticks_sec": [0, 0.5, 1, 1.5, 2, 2.5, 3],
            "yticks_pressure": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            "yticks_flux": [0, 10, 20, 30, 40, 50, 60],
        },
        "Seawater": {
            "xticks": [45, 50, 55, 60],
            "yticks_lcow": [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            "yticks_sec": [0, 2, 4, 6, 8],
            "yticks_pressure": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            "yticks_flux": [0, 10, 20, 30, 40, 50],
        },
        "Produced water": {
            "xticks": [20, 30, 40, 50, 55],
            "yticks_lcow": [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            "yticks_sec": [0, 5, 10, 15, 20, 25],
            "yticks_pressure": [0, 25, 50, 75, 100, 125, 150, 175, 200],
            "yticks_flux": [0, 10, 20, 30, 40, 50],
        },
    }
    system = FigureGenerator.get_plot_options_manager()
    system.add("CCRO", marker="d")

    system.add("1 stage RO")
    system.add("1 stage RO with ERD")
    system.add("2 stage RO with ERD")
    for case, case_opts in cases.items():
        fig = FigureGenerator().init_figure()
        fig.plot_line(dm[case, "Water recovery"], dm[case, "LCOW"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "LCOW",
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "LCOW",
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "LCOW",
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="LCOW ($\$$/m$^3$)",
            xticks=case_opts["xticks"],
            yticks=case_opts["yticks_lcow"],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_LCOW")
        # fig.show()

        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(dm[case, "Water recovery"], dm[case, "SEC"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "SEC",
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "SEC",
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "SEC",
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="SEC (kW/m$^3$)",
            xticks=case_opts["xticks"],
            yticks=case_opts["yticks_sec"],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_SEC")
        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(dm[case, "Water recovery"], dm[case, "SEC"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "SEC",
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "SEC",
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "SEC",
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="SEC (kW/m$^3$)",
            xticks=case_opts["xticks"],
            yticks=case_opts["yticks_sec"],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_SEC")
        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(dm[case, "Water recovery"], dm[case, "capex"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "capex",
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "capex",
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "capex",
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="CAPEX (MUSD)",
            xticks=case_opts["xticks"],
            yticks=[0, 0.02, 0.04, 0.06, 0.08, 0.1],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_capex")
        # assert False
        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(dm[case, "Water recovery"], dm[case, "opex"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "opex",
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "opex",
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "opex",
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="OPEX (MUSD/year)",
            xticks=case_opts["xticks"],
            yticks=[0, 0.01, 0.02, 0.03, 0.04, 0.05],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_opex")
        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(dm[case, "Water recovery"], dm[case, "RO flux"], **system["CCRO"])
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                (f"RO 1", "RO flux"),
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                ("RO 1", "RO flux"),
            ],
            **system["1 stage RO with ERD"],
        )

        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                ("RO 1", "RO flux"),
            ],
            **system["2 stage RO with ERD"],
        )
        sc = system["2 stage RO with ERD"].copy()
        sc.ls = "--"
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                ("RO 2", "RO flux"),
            ],
            **sc,
        )
        fig.plot_line([], [], color="black", ls="-", label="Stage 1 RO flux")
        fig.plot_line([], [], color="black", ls="--", label="Stage 2 RO flux")
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="RO flux (LMH)",
            xticks=case_opts["xticks"],
            yticks=case_opts["yticks_flux"],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_ro_flux")

        fig = FigureGenerator().init_figure(width=2, height=2)
        fig.plot_line(
            dm[case, "Water recovery"], dm[case, "Pump size"], **system["CCRO"]
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                ("RO 1", "pump work"),
            ],
            **system["1 stage RO"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="Feed pump size (kW)",
            xticks=case_opts["xticks"],
            yticks=[0, 5, 10, 15, 20, 25],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")
        fig.save("performance_figs", file_name=f"{case}_pump_size_comparison")
        fig = FigureGenerator().init_figure()  # width=2, height=2)
        fig.plot_line(
            dm[case, "Water recovery"], dm[case, "RO feed pressure"], **system["CCRO"]
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                ("RO 1", "feed pressure"),
            ],
            **system["1 stage RO"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "1_stage_1_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "1_stage_1_pump"),
                ("RO 1", "feed pressure"),
            ],
            **system["1 stage RO with ERD"],
        )
        fig.plot_line(
            dm_ss[
                case,
                ("add_erd", "False"),
                ("stage_sim_cases", "2_stage_2_pump"),
                "Water recovery",
            ],
            dm_ss[
                case,
                ("add_erd", "True"),
                ("stage_sim_cases", "2_stage_2_pump"),
                ("RO 2", "feed pressure"),
            ],
            **system["2 stage RO with ERD"],
        )
        fig.set_axis(
            xlabel="Water recovery (%)",
            ylabel="RO feed pressure (bar)",
            xticks=case_opts["xticks"],
            yticks=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
            # yticks=[0, 100],
        )
        fig.add_legend(loc="upper left")

        fig.save("performance_figs", file_name=f"{case}_pressure")
        # fig = FigureGenerator().init_figure(width=2, height=2)
        # fig.plot_line(
        #     dm[case, "Water recovery"], dm[case, "Pump work"], **system["CCRO"]
        # )
        # fig.plot_line(
        #     dm_ss[
        #         case,
        #         ("add_erd", "False"),
        #         ("stage_sim_cases", "1_stage_1_pump"),
        #         "Water recovery",
        #     ],
        #     dm_ss[
        #         case,
        #         ("add_erd", "False"),
        #         ("stage_sim_cases", "1_stage_1_pump"),
        #         ("RO 1", "pump work"),
        #     ],
        #     **system["1 stage RO"],
        # )
        # fig.set_axis(
        #     xlabel="Water recovery (%)",
        #     ylabel="Feed pump work (kW)",
        #     xticks=case_opts["xticks"],
        #     yticks=[0, 5, 10, 15, 20, 25],
        # )
        # fig.add_legend(loc="upper left")
        # fig.save("performance_figs", file_name=f"{case}_pump_work_comparison")
        # fig.show()
