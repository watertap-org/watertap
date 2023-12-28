from binascii import b2a_base64
from tkinter import Y
from adam_bsm2 import main
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
    op1 =[]
    for i in range(len(inlet)):
        ms.from_json(m, fname="modelsolve_flow.json.gz")
        m.fs.FeedWater.flow_vol.fix(inlet[i] * pyo.units.m**3 / pyo.units.day)
        results = solve(m)
        x.append(pyo.value(m.fs.FeedWater.flow_vol[0])*3600*24)
        lco.append(pyo.value(m.fs.costing.LCOW))
        cap.append(pyo.value(m.fs.costing.total_capital_cost))
        oper.append(pyo.value(m.fs.costing.total_operating_cost))
        anual.append(pyo.value(m.fs.costing.total_annualized_cost))


        unit = [pyo.value(
            m.fs.costing.factor_capital_annualization
            * m.fs.CL.costing.capital_cost
            / m.fs.costing.annual_water_production),
            pyo.value(m.fs.costing.factor_capital_annualization
            * (m.fs.R1.costing.capital_cost +
            m.fs.R2.costing.capital_cost+
            m.fs.R3.costing.capital_cost+
            m.fs.R4.costing.capital_cost+
            m.fs.R5.costing.capital_cost+
                m.fs.CL1.costing.capital_cost)
                / m.fs.costing.annual_water_production),
            pyo.value(m.fs.costing.factor_capital_annualization
            * m.fs.RADM.costing.capital_cost
            / m.fs.costing.annual_water_production),
            pyo.value(m.fs.costing.factor_capital_annualization
            * m.fs.DU.costing.capital_cost
            / m.fs.costing.annual_water_production),
            pyo.value(m.fs.costing.factor_capital_annualization
            * m.fs.TU.costing.capital_cost
            / m.fs.costing.annual_water_production),
        ]

        # y = [pyo.value(m.fs.CL.electricity_consumption[0]),
        #     pyo.value(m.fs.R3.electricity_consumption[0]) +
        #     pyo.value(m.fs.R4.electricity_consumption[0]) +
        #     pyo.value(m.fs.R5.electricity_consumption[0]) +
        #     pyo.value(m.fs.CL1.electricity_consumption[0]),
        #     pyo.value(m.fs.RADM.electricity_consumption[0]),
        #     pyo.value(m.fs.DU.electricity_consumption[0]),
        #     pyo.value(m.fs.TU.electricity_consumption[0]),
        # ]

        op = [(pyo.value((m.fs.costing.aggregate_fixed_operating_cost
            + m.fs.costing.maintenance_labor_chemical_operating_cost )/ m.fs.costing.annual_water_production)),
        (pyo.value((m.fs.costing.aggregate_variable_operating_cost
                + sum(m.fs.costing.aggregate_flow_costs[f] for f in m.fs.costing.used_flows)
                * m.fs.costing.utilization_factor)/ m.fs.costing.annual_water_production))]
        unit1.append(unit)
        op1.append(op)

    unit2=np.zeros((10,5))
    rows = 10
    columns = 5
    for m in range(rows):
        for i in range(columns):
            unit2[m][i] = unit1[m][i]

    op2=np.zeros((10,5))
    rows = 10
    columns = 2
    for m in range(rows):
        for i in range(columns):
            op2[m][i] = op1[m][i]

    print(op2)

    labels = 'Primary', 'Secondary', 'Digestion', 'dewatering', 'thickening'

    Primary = unit2[:,0]
    Secondary = unit2[:,1]
    Digestion = unit2[:,2]
    dewatering = unit2[:,3]
    thickening = unit2[:,4]
    # Primary_e = y2[:,0]
    # Secondary_e = y2[:,1]
    # Digestion_e = y2[:,2]
    # dewatering_e = y2[:,3]
    # thickening_e = y2[:,4]
    operating_fixed = op2[:,0]
    operating_variable = op2[:,1]


    groups = ["18000", "21555", "25111", "28666", "32222", "35777", "39333", "42888", "46444", "50000"]
    b1 = np.add(Primary, Secondary)
    b2 = np.add(b1, Digestion)
    b3 = np.add(b2, dewatering)
    b4 = np.add(b3, thickening)
    # b5 = np.add(b4, Primary_e)
    # b6 = np.add(b5, Secondary_e)
    # b7 = np.add(b6, Digestion_e)
    # b8 = np.add(b7, dewatering_e)
    b5 = np.add(b4, operating_fixed)
    b6 = np.add(b5, operating_variable)


    fig, ax = plt.subplots(figsize=(15, 10))
    ax.bar(groups, Primary, color = "navy", label ='Primary treatment capital cost')
    ax.bar(groups, Secondary, bottom = Primary, color = "blue", label ='Secondary treatment capital cost')
    ax.bar(groups, Digestion, bottom = b1, color = "dodgerblue", label = 'Anaerobic Digestion  capital cost')
    ax.bar(groups, dewatering, bottom = b2, color = "skyblue", label ='Dewatering unit capital cost')
    ax.bar(groups, thickening, bottom = b3, color = "mediumaquamarine", label ='Thickener capital cost')
    ax.bar(groups, operating_fixed, bottom = b4, color = "forestgreen", label ='Fixed operating cost')
    ax.bar(groups, operating_variable, bottom = b5, color = "gold", label ='Variable operating cost')
    # ax.bar(groups, Primary_e, bottom = b4, color = "mediumseagreen", label ='Primary treatment operating cost')
    # ax.bar(groups, Secondary_e, bottom = b5, color = "forestgreen", label ='Secondary treatment operating cost')
    # ax.bar(groups, Digestion_e, bottom = b6, color = "gold", label ='Anaerobic Digestion  operating cost')
    # ax.bar(groups, dewatering_e, bottom = b7, color = "orange", label ='Dewatering unit operating cost')
    # ax.bar(groups, thickening_e, bottom = b8, color = "tomato", label ='Thickener operating cost')
    # plt.xticks(x)
    ax.set_ylabel('LCOW \$/m$^3$', fontsize=25, labelpad=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.legend(loc = 'lower left', fontsize=20)
    plt.savefig("cost_by_flow.png", bbox_inches="tight", dpi=300)
    plt.show() 
 
    

    oper1 = []
    anual1 =  []

    for i in range(len(oper)):
        oper1.append(oper[i] /365 /inlet[i])

    for i in range(len(anual)):
        anual1.append(anual[i] /365 /inlet[i])

    fig2, ax = plt.subplots(figsize=(15, 10))
    lineslco = ax.plot(x, oper1, label="Operating cost")
    lineslco = ax.plot(x, anual1, label="Capital cost")
    lineslco = ax.plot(x, lco, label="LCOW")
    ax.set_xlabel('Plant flow [m$^3$/day]', fontsize=25, labelpad=20)
    ax.set_ylabel('Costs \$/m$^3$', fontsize=25, labelpad=20)
    # ax.set_title('LCOW', fontsize=34)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.legend(fontsize=20)
    plt.savefig("cost_plot.png", bbox_inches="tight", dpi=300)
    plt.show()  
