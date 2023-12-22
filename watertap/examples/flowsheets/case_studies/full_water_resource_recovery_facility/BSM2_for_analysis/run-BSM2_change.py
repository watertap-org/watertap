from tkinter import Y
from BSM2_changes import main

import pyomo.environ as pyo

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import numpy as np

import matplotlib.pyplot as plt


# Import python path
import os

# Import idaes model serializer to store initialized model
from idaes.core.util import model_serializer as ms 


def solve(blk, solver=None):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=False)
    pyo.assert_optimal_termination(results)
    return results

# Cost breakdown
if __name__ == '__main__':  
    m, results = main()

    labels = 'Primary', 'Secondary', 'Digestion', 'dewatering', 'thickening'


    tss = np.linspace(0.25, 0.15, num=3)
    x1=[]

    for i in range(len(tss)):
        # ms.from_json(m, fname="modelsolve.json.gz")
        m.fs.totalN_max.fix = tss[i]

        results = solve(m)

        x = [pyo.value(m.fs.CL.costing.capital_cost),
            pyo.value(m.fs.R1.costing.capital_cost) +
            pyo.value(m.fs.R2.costing.capital_cost) +
            pyo.value(m.fs.R3.costing.capital_cost) +
            pyo.value(m.fs.R4.costing.capital_cost) +
            pyo.value(m.fs.R5.costing.capital_cost) +
            pyo.value(m.fs.CL1.costing.capital_cost),
            pyo.value(m.fs.RADM.costing.capital_cost),
            pyo.value(m.fs.DU.costing.capital_cost),
            pyo.value(m.fs.TU.costing.capital_cost),
        ]
        x1.append(x)
        print("soooooolved")

    x2=np.zeros((5,5))
    rows = 3
    columns = 5
    print(x1)
    for m in range(rows):
        for i in range(columns):
            print(m)
            print(i)
            x2[m][i] = x1[m][i] /3.154e+7/10*1.08
    print(x2)
    # print(sum( x2[i] for i in range(len(x2))))
    # fig, ax = plt.subplots()
    # ax.pie(x, labels=labels, autopct='%1.1f%%')

    # plt.show()


