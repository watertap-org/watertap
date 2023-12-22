from binascii import b2a_base64
from tkinter import Y
from adam_bsm2_sim_vs_opt import main
from adam_bsm2_sim_vs_opt import main2
import matplotlib.colors as mcolors

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

if __name__ == '__main__':  
    m, results = main()
    if not os.path.exists("modelsolve_flow.json.gz"):
        ms.to_json(m, fname="modelsolve_flow.json.gz")
    else:
        ms.from_json(m, fname="modelsolve_flow.json.gz")
    # cost = np.zeros((1,8))
    inlet = np.linspace(18000, 50000, num=10)
    x=[]
    lco=[]
    cap = []
    anual = []
    oper = []
    unit = []
    unit1 = []
    y = []
    y1 = []


    x.append(pyo.value(m.fs.FeedWater.flow_vol[0])*3600*24)
    lco.append(pyo.value(m.fs.costing.LCOW))
    cap.append(pyo.value(m.fs.costing.total_capital_cost))
    oper.append(pyo.value(m.fs.costing.total_operating_cost))
    anual.append(pyo.value(m.fs.costing.total_annualized_cost))


    unit = [pyo.value(m.fs.CL.costing.capital_cost),
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

    y = [pyo.value(m.fs.CL.electricity_consumption[0]),
        pyo.value(m.fs.R3.electricity_consumption[0]) +
        pyo.value(m.fs.R4.electricity_consumption[0]) +
        pyo.value(m.fs.R5.electricity_consumption[0]) +
        pyo.value(m.fs.CL1.electricity_consumption[0]),
        pyo.value(m.fs.RADM.electricity_consumption[0]),
        pyo.value(m.fs.DU.electricity_consumption[0]),
        pyo.value(m.fs.TU.electricity_consumption[0]),
    ]
    unit1.append(unit)
    y1.append(y)

    m, results = main2()

    x.append(pyo.value(m.fs.FeedWater.flow_vol[0])*3600*24)
    lco.append(pyo.value(m.fs.costing.LCOW))
    cap.append(pyo.value(m.fs.costing.total_capital_cost))
    oper.append(pyo.value(m.fs.costing.total_operating_cost))
    anual.append(pyo.value(m.fs.costing.total_annualized_cost))


    unit = [pyo.value(m.fs.CL.costing.capital_cost),
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

    y = [pyo.value(m.fs.CL.electricity_consumption[0]),
        pyo.value(m.fs.R3.electricity_consumption[0]) +
        pyo.value(m.fs.R4.electricity_consumption[0]) +
        pyo.value(m.fs.R5.electricity_consumption[0]) +
        pyo.value(m.fs.CL1.electricity_consumption[0]),
        pyo.value(m.fs.RADM.electricity_consumption[0]),
        pyo.value(m.fs.DU.electricity_consumption[0]),
        pyo.value(m.fs.TU.electricity_consumption[0]),
    ]
    unit1.append(unit)
    y1.append(y)

    unit2=np.zeros((2,5))
    rows = 2
    columns = 5
    for m in range(rows):
        for i in range(columns):
            unit2[m][i] = unit1[m][i] * 1.08 /365 /10 /inlet[m]

    y2=np.zeros((2,5))
    rows = 2
    columns = 5
    for m in range(rows):
        for i in range(columns):
            y2[m][i] = y1[m][i] * 0.07 /(inlet[m]/24)

    labels = 'Primary', 'Secondary', 'Digestion', 'dewatering', 'thickening'

    Primary = unit2[:,0]
    Secondary = unit2[:,1]
    Digestion = unit2[:,2]
    dewatering = unit2[:,3]
    thickening = unit2[:,4]
    Primary_e = y2[:,0]
    Secondary_e = y2[:,1]
    Digestion_e = y2[:,2]
    dewatering_e = y2[:,3]
    thickening_e = y2[:,4]

    groups = ["sim", "opt"]
    b1 = np.add(Primary, Secondary)
    b2 = np.add(b1, Digestion)
    b3 = np.add(b2, dewatering)
    b4 = np.add(b3, thickening)
    b5 = np.add(b4, Primary_e)
    b6 = np.add(b5, Secondary_e)
    b7 = np.add(b6, Digestion_e)
    b8 = np.add(b7, dewatering_e)


    fig, ax = plt.subplots()
    ax.bar(groups, Primary, color = "navy", label ='Primary treatment capital cost')
    ax.bar(groups, Secondary, bottom = Primary, color = "blue", label ='Secondary treatment capital cost')
    ax.bar(groups, Digestion, bottom = b1, color = "dodgerblue", label = 'Anaerobic Digestion  capital cost')
    ax.bar(groups, dewatering, bottom = b2, color = "skyblue", label ='Dewatering unit capital cost')
    ax.bar(groups, thickening, bottom = b3, color = "mediumaquamarine", label ='Thickener capital cost')
    ax.bar(groups, Primary_e, bottom = b4, color = "mediumseagreen", label ='Primary treatment operating cost')
    ax.bar(groups, Secondary_e, bottom = b5, color = "forestgreen", label ='Secondary treatment operating cost')
    ax.bar(groups, Digestion_e, bottom = b6, color = "gold", label ='Anaerobic Digestion  operating cost')
    ax.bar(groups, dewatering_e, bottom = b7, color = "orange", label ='Dewatering unit operating cost')
    ax.bar(groups, thickening_e, bottom = b8, color = "tomato", label ='Thickener operating cost')
    # plt.xticks(x)
    plt.ylabel('LCOW \$/m$^3$')
    plt.title('Simulation vs. Optimization')
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show() 