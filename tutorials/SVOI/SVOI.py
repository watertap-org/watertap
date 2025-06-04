import os
import time
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
from idaes.core.solvers import get_solver
from pyomo.environ import (
    check_optimal_termination,
)
import example_flowsheet as RO
import yaml
import pandas as pd

from psPlotKit.data_manager import (
    data_importer,
)
from psPlotKit.data_manager.ps_data_manager import psDataManager
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def build_sweep(
    inlet_flow_rate=100,
    inlet_concentration=35,
    inlet_pressure=50e5,
    membrane_area=10000,
):
    m = RO.build_system()
    RO.set_operating_conditions(m, Qin=inlet_flow_rate, Cin=inlet_concentration)
    RO.set_ro_system_operating_conditions(
        m, mem_area=membrane_area, RO_pressure=inlet_pressure
    )
    RO.set_scaling(m)
    RO.init_system(m)
    RO.add_costing(m)
    RO.solve(m)
    RO.report_RO(m)
    RO.optimize(
        m,
        water_recovery=0.5,
    )

    return m


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results


cwd = get_working_dir()
yaml_file = os.path.join(cwd, "tutorials", "SVOI", "diff_param_sweep.yaml")
save_name = "sweep"
save_dir = os.path.join(cwd, "tutorials", "SVOI", "sweep")

print(yaml_file)
print(save_dir)

solver = get_solver()
solver.options["max_iter"] = 3000

# lT = loopTool(
#     yaml_file,
#     build_function=build_sweep,
#     optimize_function=solve,
#     solver=solver,
#     save_name=save_name,
#     saving_dir=save_dir,
#     execute_simulations=True,
#     number_of_subprocesses=1,
# )


file = os.path.join(save_dir, "output", "sweep_analysisType_example_diff_analysis.h5")
data_manager = psDataManager(file)

VOI = data_manager.get_voi(yaml_file)

# print(VOI)

frames = [pd.DataFrame({"VOI": VOI[key]["VOI"], "Sweep": key}) for key in VOI.keys()]
df = pd.concat(frames, ignore_index=True)


def plot_voi(df):
    fig, ax = plt.subplots(figsize=(7, 9))
    sns.boxplot(
        df,
        x="VOI",
        y="Sweep",
        gap=0.1,
        fliersize=0,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(df["Sweep"].unique()) - 0.5)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_ylabel("")
    ax.set_xlabel("VOI (%$_{\Delta LCOW}$/%$_{\Delta performance}$)", fontsize=16)
    color_bands = ["#d1d1d1", "white"]
    for i in range(0, len(ax.get_yticks())):
        ax.axhspan(i - 0.5, i + 0.5, facecolor=color_bands[i % 2], alpha=0.75)

    plt.legend(
        bbox_to_anchor=(-0.3, 1.1),
        loc="upper left",
        ncol=len(df["Sweep"].unique()),
        frameon=False,
        fontsize=14,
        alignment="left",
        handlelength=1,
        handleheight=1,
    )

    plt.tight_layout()
    plt.show()


yaml_file = os.path.join(cwd, "tutorials", "SVOI", "diff_param_sweep_full.yaml")

# lT = loopTool(
#         yaml_file,
#         build_function=build_sweep,
#         optimize_function=solve,
#         solver=solver,
#         save_name=save_name,
#         saving_dir=save_dir,
#         execute_simulations=True,
#         number_of_subprocesses=1, # For some reason my cpu is not working > 1
#     )

# Load results file
results_file = os.path.join(
    save_dir, "output", "sweep_analysisType_example_diff_analysis_full.h5"
)

data_manager = psDataManager(results_file)
VOI = data_manager.get_voi(yaml_file)
frames = [pd.DataFrame({"VOI": VOI[key]["VOI"], "Sweep": key}) for key in VOI.keys()]
df = pd.concat(frames, ignore_index=True)

plot_voi(df)
