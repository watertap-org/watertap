from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.fig_generator import figureGenerator


if __name__ == "__main__":
    ##################################################################################
    # Data processing workflow follows plotting or tornda plot:
    # 1) Load data from h5 file using PsDataManager
    # 2) Create a dictionary of baseline values and their names for import
    #    this is required, each group should have general name, and unit keys, if they are indexed all indexed units will be combined
    # 3) Register data keys for sensitivities
    # 4) Normalize data (E.g. get change in LCOW)
    # 5) Build bar plot useinf FigureGenerator
    ###################################################################################

    data_manager = PsDataManager(
        [
            "output/lsrro_sensitivity_sweep_analysisType_tornado_sweep.h5",
        ]
    )
    data_manager.register_data_key("fs.water_recovery", "Water recovery", "%")
    data_manager.register_data_key(
        "fs.feed.feed.properties[0].conc_mass_phase_comp[Liq, TDS]", "Feed TDS", "g/L"
    )
    data_manager.register_data_key("fs.costing.LCOW", "LCOW")

    # the baseline_value here is 1, but if we are using absolute value (e.g. 85 bar, then we would use this for normalizations)
    sensitivity_baselines = {
        "Max. pump pressure": {
            "baseline_val": 1,
            "directory": "pump_pressure_multiplier",
        },
        "Water permeability": {
            "baseline_val": 1,
            "directory": "a_value_multiplier",
        },
        "LSRRO mem. cost": {
            "baseline_val": 1,
            "directory": "lsrro_membrane_cost_multiplier",
        },
        "Pump cost": {
            "baseline_val": 1,
            "directory": "pump_cost_multiplier",
        },
        "Pump efficiency": {
            "baseline_val": 1,
            "directory": "pump_efficiency_multiplier",
        },
    }
    _sens_base = {}
    for sense_name, sense_info in sensitivity_baselines.items():
        data_manager.register_data_key(
            f"fs.sense_manager.{sense_info['directory']}",
            sense_name,
            directories=sense_info["directory"],
        )
        _sens_base[sense_name] = sense_info["baseline_val"]

    data_manager.load_data()
    data_manager.display()

    # This will normalize the data, it will take the keys in the dictory, find specified
    # sensetivity var, and "center point" which is the baseline value, and calculate change of LCOW
    # relative from center point.
    # e.g. if you do mutlplier then Max. pump pressure will be 0.8, 1, and 1.2, this will find LCOW  a 1 and use that
    # to get change for other values.
    data_manager.normalize_data(base_value_dict=_sens_base, related_keys="LCOW")
    data_manager.display()

    treatment_scenarios = {
        "case_a": {
            "label": "Case A (65 g/L TDS to 70% WR)",
            "color": "#7db6df",
        },
        "case_b": {
            "label": "Case B (125 g/L TDS to 50% WR)",
            "color": "#ff5d5d",
        },
    }
    fig = figureGenerator()
    fig.init_figure()

    for j, sense in enumerate(sensitivity_baselines):
        for i, (case_key, case_value) in enumerate(treatment_scenarios.items()):
            pos = j + i * 0.4 - 0.2
            _sense_data = data_manager[
                (
                    ("sim_cases", case_key),
                    sensitivity_baselines[sense]["directory"],
                    "LCOW",
                )
            ]
            sense_data = data_manager[
                (
                    ("sim_cases", case_key),
                    sensitivity_baselines[sense]["directory"],
                    "normalized_data",
                    "LCOW",
                )
            ]
            sense_steps = data_manager[
                (
                    ("sim_cases", case_key),
                    sensitivity_baselines[sense]["directory"],
                    "normalized_data",
                    sense,
                )
            ]
            print(sense, case_key, sense_data.data, sense_steps.data, _sense_data.data)

            # We know first value is [0], center point is index [1] and high value is [2]
            fig.plot_bar(
                [pos],
                [sense_data.data[0]],
                bottom=0,
                hatch="///",
                color=case_value["color"],
                vertical=False,
                width=0.4,
            )
            fig.plot_bar(
                [pos],
                [sense_data.data[2]],
                bottom=0,
                hatch="",
                color=case_value["color"],
                vertical=False,
                width=0.4,
            )
            print(j, pos)
            if j % 2 == 0:
                fig.plot_area(
                    [-120, -120],
                    [pos - 0.3, pos + 0.3],
                    x2data=[60, 60],
                    color="#d6d6d6b2",
                    edgecolor=None,
                    zorder=-10,
                    lw=1,
                    clip_on=False,
                )
    for case, case_value in treatment_scenarios.items():
        fig.plot_bar(
            [-10],
            [0],
            hatch="",
            color=case_value["color"],
            label=case_value["label"],
        )

    fig.plot_bar(
        [-10],
        [0],
        hatch="///",
        color="white",
        label="-20%",
    )
    fig.plot_bar(
        [-10],
        [0],
        hatch="",
        color="white",
        label="+20%",
    )
    fig.set_axis(
        xlims=[-30, 60],
        xlabel="Change in LCOW (%)",
        xticks=[-30, -20, -10, 0, 10, 20, 30, 40, 50, 60],
    )
    fig.set_axis_ticklabels(
        yticks=list(range(len(sensitivity_baselines))),
        ylims=[
            len(sensitivity_baselines) - 0.5,
            -0.5,
        ],
        yticklabels=[key for key, label in sensitivity_baselines.items()],
    )
    fig.add_legend(loc="upper center", bbox_to_anchor=[0.5, 1.3])
    fig.save("figures", "tornado_plot")
    fig.show()
