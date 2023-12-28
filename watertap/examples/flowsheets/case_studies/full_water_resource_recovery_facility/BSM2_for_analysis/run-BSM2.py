from tkinter import Y
from adam_bsm2 import main

import pyomo.environ as pyo

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import numpy as np

import matplotlib.pyplot as plt


# Import python path
import os

# Import idaes model serializer to store initialized model
from idaes.core.util import model_serializer as ms 
from pyomo.opt import TerminationCondition, SolverStatus

solver = get_solver(options={"bound_push": 1e-8})

def solve(blk, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=False)
    pyo.assert_optimal_termination(results)
    return results

# LCOW vs inlet SS
if __name__ == '__main__':  
    m, results = main()
    if not os.path.exists("modelsolve.json.gz"):
        ms.to_json(m, fname="modelsolve.json.gz")
    else:
        ms.from_json(m, fname="modelsolve.json.gz")
    cost = np.zeros((11))
    inlet = np.linspace(70, 50, num=10)
    x=[]
    z=[]
    for j in range(len(inlet)):
        m.fs.FeedWater.conc_mass_comp[0, "S_S"].fix(inlet[j] * pyo.units.g / pyo.units.m**3)
        try:

            # solve the model
            status = solver.solve(m, tee=False)

            if (status.solver.status == SolverStatus.ok) and (
                status.solver.termination_condition == TerminationCondition.optimal):
            
                cost[j] = (pyo.value(m.fs.costing.LCOW))
                x.append(pyo.value(m.fs.FeedWater.conc_mass_comp[0, "S_S"]))
                z.append(cost[j])
        except ValueError:
            pass
    print(cost)

    fig, ax = plt.subplots(figsize=(15, 10))
    plt.plot(x, z)
    ax.set_xlabel("Inlet S$_S$ (kg/m$^3$)", fontsize=25, labelpad=20)
    ax.set_ylabel("LCOW \$/m$^3$", fontsize=25, labelpad=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig("S_S.png", bbox_inches="tight", dpi=300)
    plt.show()

# # LCOW vs inlet SS
# if __name__ == '__main__':  
#     m, results = main()
#     if not os.path.exists("modelsolve.json.gz"):
#         ms.to_json(m, fname="modelsolve.json.gz")
#     else:
#         ms.from_json(m, fname="modelsolve.json.gz")
#     cost = np.zeros((10))
#     inlet = np.linspace(60, 130.5, num=10)
#     x=[]
#     z=[]
#     for j in range(len(inlet)):
#         m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(inlet[j] * pyo.units.g / pyo.units.m**3)
#         try:

#             # solve the model
#             status = solver.solve(m, tee=False)

#             if (status.solver.status == SolverStatus.ok) and (
#                 status.solver.termination_condition == TerminationCondition.optimal):

#                 cost[j] = (pyo.value(m.fs.costing.LCOW))
#                 x.append(pyo.value(m.fs.FeedWater.conc_mass_comp[0, "X_I"]))
#                 z.append(cost[j])
#         except ValueError:
#             pass
#     print(cost)

    # fig, ax = plt.subplots(figsize=(15, 10))
    # plt.plot(x, z)
    # ax.set_xlabel("Inlet X$_I$ (kg/m$^3$)", fontsize=25, labelpad=20)
    # ax.set_ylabel("LCOW \$/m$^3$", fontsize=25, labelpad=20)
    # plt.xticks(fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.savefig("X_I.png", bbox_inches="tight", dpi=300)
    # plt.show()




