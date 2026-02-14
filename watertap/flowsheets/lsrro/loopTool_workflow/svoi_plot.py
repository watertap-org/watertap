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
            "output/lsrro_sensitivity_sweep_analysisType_svoi_sweep.h5",
        ]
    )
    data_manager.register_data_key("fs.water_recovery", "Water recovery", "%")
    data_manager.register_data_key(
        "fs.feed.feed.properties[0].conc_mass_phase_comp[Liq, TDS]", "Feed TDS", "g/L"
    )
    data_manager.register_data_key("fs.costing.LCOW", "LCOW")

    # This is for SVOI Specifically
    data_manager.register_data_key("nominal_idx", "nominal_idx")
    data_manager.register_data_key("differential_idx", "differential_idx")

    # Dict of registering keys
    sensitivity_baselines = {
        "Max. pump pressure": {
            "directory": "pump_pressure",
        },
        "Water permeability": {
            "directory": "a_value",
        },
        "LSRRO mem. cost": {
            "directory": "lsrro_membrane_cost",
        },
        "Pump cost": {
            "directory": "pump_cost",
        },
        "Pump efficiency": {
            "directory": "pump_efficiency",
        },
    }

    for sense_name, sense_info in sensitivity_baselines.items():
        data_manager.register_data_key(
            f"fs.sense_manager.{sense_info['directory']}",
            sense_name,
            directories=sense_info["directory"],
        )

    data_manager.load_data()
    data_manager.display()

    # define SVOI function do calculatin,
    # We need to calcualte chnage in LCOW between base optimziation (referneced by nominal_idx), and
    # simulation with improvment (differential_idx)
    # in our SVOI setup file we are doing 1 percentile step, as such change in LCOW is
    # the percent change in LCOW/percentile improvment in parameter
    def svoi_calc(cost, nominal_idx, diff_idx):
        nominal_cost = cost[nominal_idx == nominal_idx]
        diff_cost = cost[diff_idx == diff_idx]

        # invert it so its positive (LCOW improvement!)
        cost_change_percent = -1 * (diff_cost - nominal_cost) / nominal_cost * 100
        # assert False
        return cost_change_percent

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

    # we need to now run calculation on SVOI, this can be done using
    # data manager, by creating a dictionary that conntains the directories to be passed into
    # our function, and then calling eval_function on data manager
    # it will pull out the data for each directory and directly pass it into the function
    for case in treatment_scenarios:
        for sense in sensitivity_baselines:
            voi_dict = {
                "cost": (
                    ("sim_cases", case),
                    sensitivity_baselines[sense]["directory"],
                    "LCOW",
                ),
                "nominal_idx": (
                    ("sim_cases", case),
                    sensitivity_baselines[sense]["directory"],
                    "nominal_idx",
                ),
                "diff_idx": (
                    ("sim_cases", case),
                    sensitivity_baselines[sense]["directory"],
                    "differential_idx",
                ),
            }
            data_manager.eval_function(
                (("sim_cases", case), sensitivity_baselines[sense]["directory"]),
                "voi",
                svoi_calc,
                voi_dict,
            )
            data_manager[
                ("sim_cases", case),
                sensitivity_baselines[sense]["directory"],
                "voi",
            ].display()

    data_manager.display()

    fig = figureGenerator()
    fig.init_figure()

    for j, sense in enumerate(sensitivity_baselines):
        for i, (case_key, case_value) in enumerate(treatment_scenarios.items()):
            pos = j + i * 0.4 - 0.2
            voi_data = data_manager[
                ("sim_cases", case_key),
                sensitivity_baselines[sense]["directory"],
                "voi",
            ]
            fig.plot_box(
                pos,
                voi_data.get_data(exclude_nan_values=True),
                color=case_value["color"],
                vertical=False,
                width=0.4,
            )
            if j % 2 == 0:
                fig.plot_area(
                    [-0.6, -0.6],
                    [pos - 0.3, pos + 0.3],
                    x2data=[1, 1],
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
    fig.set_axis(
        xlims=[0, 1],
        xlabel="SVOI (%$_{\Delta LCOW}$ / %$_{\Delta performance}$)",
        xticks=[0, 0.2, 0.4, 0.6, 0.8, 1],
    )
    fig.set_axis_ticklabels(
        yticks=list(range(len(sensitivity_baselines))),
        ylims=[
            len(sensitivity_baselines) - 0.5,
            -0.5,
        ],
        yticklabels=[key for key, label in sensitivity_baselines.items()],
    )
    fig.add_legend(loc="upper center", bbox_to_anchor=[0.5, 1.2])
    fig.save("figures", "svoi_plot")
    fig.show()
