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
    inlet = np.linspace(20648, 20648, num=10)
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
    op = []
    op1 = []

    print("clarifier", pyo.value(
        m.fs.costing.factor_capital_annualization
        * m.fs.CL.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume))

    print("ASM", pyo.value(m.fs.costing.factor_capital_annualization
        * (m.fs.R1.costing.capital_cost +
           m.fs.R2.costing.capital_cost+
           m.fs.R3.costing.capital_cost+
           m.fs.R4.costing.capital_cost+
           m.fs.R5.costing.capital_cost)
            / m.fs.costing.annual_feedwater_volume))

    print("anaerobic", pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.RADM.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume))
    
    print("secondary cl", pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.CL1.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume))

    print("thickener", pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.TU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume))

    print("dewaterind", pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.DU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume))

    print("electricity", pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
        * m.fs.costing.utilization_factor
        / m.fs.costing.annual_feedwater_volume
    ))
    print("electricity", pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
    ))

    print("opex",pyo.value(m.fs.costing.maintenance_labor_chemical_operating_cost
        / m.fs.costing.annual_feedwater_volume
    ))

    print("fixed",(pyo.value((m.fs.costing.aggregate_fixed_operating_cost
            + m.fs.costing.maintenance_labor_chemical_operating_cost )/ m.fs.costing.annual_feedwater_volume)))
    print("variable",(pyo.value((m.fs.costing.aggregate_variable_operating_cost
                + sum(m.fs.costing.aggregate_flow_costs[f] for f in m.fs.costing.used_flows)
                * m.fs.costing.utilization_factor)/ m.fs.costing.annual_feedwater_volume)))

    unit = [pyo.value(
        m.fs.costing.factor_capital_annualization
        * m.fs.CL.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * (m.fs.R1.costing.capital_cost +
           m.fs.R2.costing.capital_cost+
           m.fs.R3.costing.capital_cost+
           m.fs.R4.costing.capital_cost+
           m.fs.R5.costing.capital_cost+
            m.fs.CL1.costing.capital_cost)
            / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.RADM.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.DU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.TU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
    ]

    # y = [pyo.value(m.fs.CL.electricity_consumption[0])+o1,
    #     pyo.value(m.fs.R3.electricity_consumption[0]) +
    #     pyo.value(m.fs.R4.electricity_consumption[0]) +
    #     pyo.value(m.fs.R5.electricity_consumption[0]) +
    #     pyo.value(m.fs.CL1.electricity_consumption[0])+o2,
    #     pyo.value(m.fs.RADM.electricity_consumption[0])+o3,
    #     pyo.value(m.fs.DU.electricity_consumption[0])+o4,
    #     pyo.value(m.fs.TU.electricity_consumption[0])+o5,
    # ]

    op = [(pyo.value((m.fs.costing.aggregate_fixed_operating_cost
            + m.fs.costing.maintenance_labor_chemical_operating_cost )/ m.fs.costing.annual_feedwater_volume)),
        (pyo.value((m.fs.costing.aggregate_variable_operating_cost
                + sum(m.fs.costing.aggregate_flow_costs[f] for f in m.fs.costing.used_flows)
                * m.fs.costing.utilization_factor)/ m.fs.costing.annual_feedwater_volume))]
    unit1.append(unit)
    # y1.append(y)
    op1.append(op)

    m, results = main2()

    x.append(pyo.value(m.fs.FeedWater.flow_vol[0])*3600*24)
    lco.append(pyo.value(m.fs.costing.LCOW))
    cap.append(pyo.value(m.fs.costing.total_capital_cost))
    oper.append(pyo.value(m.fs.costing.total_operating_cost))
    anual.append(pyo.value(m.fs.costing.total_annualized_cost))


    unit = [pyo.value(
        m.fs.costing.factor_capital_annualization
        * m.fs.CL.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * (m.fs.R1.costing.capital_cost +
           m.fs.R2.costing.capital_cost+
           m.fs.R3.costing.capital_cost+
           m.fs.R4.costing.capital_cost+
           m.fs.R5.costing.capital_cost+
            m.fs.CL1.costing.capital_cost)
            / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.RADM.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.DU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
        pyo.value(m.fs.costing.factor_capital_annualization
        * m.fs.TU.costing.capital_cost
        / m.fs.costing.annual_feedwater_volume),
    ]

    # o1 = (pyo.value(
    #     m.fs.costing.factor_capital_annualization
    #     * m.fs.CL.costing.capital_cost
    #     / m.fs.costing.annual_feedwater_volume)* 0.03)
    # o2 = (pyo.value(m.fs.costing.factor_capital_annualization
    #     * (m.fs.R1.costing.capital_cost +
    #        m.fs.R2.costing.capital_cost+
    #        m.fs.R3.costing.capital_cost+
    #        m.fs.R4.costing.capital_cost+
    #        m.fs.R5.costing.capital_cost+
    #         m.fs.CL1.costing.capital_cost)
    #         / m.fs.costing.annual_feedwater_volume)* 0.03)
    # o3 = (pyo.value(m.fs.costing.factor_capital_annualization
    #     * m.fs.RADM.costing.capital_cost
    #     / m.fs.costing.annual_feedwater_volume)* 0.03)
    # o4 = (pyo.value(m.fs.costing.factor_capital_annualization
    #     * m.fs.DU.costing.capital_cost
    #     / m.fs.costing.annual_feedwater_volume)* 0.03)
    # o5 = (pyo.value(m.fs.costing.factor_capital_annualization
    #     * m.fs.TU.costing.capital_cost
    #     / m.fs.costing.annual_feedwater_volume)* 0.03)

    # y = [pyo.value(m.fs.CL.electricity_consumption[0])+o1,
    #     pyo.value(m.fs.R3.electricity_consumption[0]) +
    #     pyo.value(m.fs.R4.electricity_consumption[0]) +
    #     pyo.value(m.fs.R5.electricity_consumption[0]) +
    #     pyo.value(m.fs.CL1.electricity_consumption[0])+o2,
    #     pyo.value(m.fs.RADM.electricity_consumption[0])+o3,
    #     pyo.value(m.fs.DU.electricity_consumption[0])+o4,
    #     pyo.value(m.fs.TU.electricity_consumption[0])+o5,
    # ]
    op = [(pyo.value((m.fs.costing.aggregate_fixed_operating_cost
            + m.fs.costing.maintenance_labor_chemical_operating_cost )/ m.fs.costing.annual_feedwater_volume)),
        (pyo.value((m.fs.costing.aggregate_variable_operating_cost
                + sum(m.fs.costing.aggregate_flow_costs[f] for f in m.fs.costing.used_flows)
                * m.fs.costing.utilization_factor)/ m.fs.costing.annual_feedwater_volume))]
    unit1.append(unit)
    # y1.append(y)
    op1.append(op)

    unit1.append(unit)
    # y1.append(y)

    print("LCOW cap",pyo.value(
        # m.fs.costing.primary_pump_capex_lcow = Expression(
        # expr=
        m.fs.costing.factor_capital_annualization *(
        pyo.value(m.fs.CL.costing.capital_cost) +
        pyo.value(m.fs.R1.costing.capital_cost) +
        pyo.value(m.fs.R2.costing.capital_cost) +
        pyo.value(m.fs.R3.costing.capital_cost) +
        pyo.value(m.fs.R4.costing.capital_cost) +
        pyo.value(m.fs.R5.costing.capital_cost) +
        pyo.value(m.fs.CL1.costing.capital_cost) +
        pyo.value(m.fs.RADM.costing.capital_cost) +
        pyo.value(m.fs.DU.costing.capital_cost) +
        pyo.value(m.fs.TU.costing.capital_cost) )
        / m.fs.costing.annual_feedwater_volume
    ))
    print("LCOW direct",pyo.value(
        # m.fs.costing.primary_pump_capex_lcow = Expression(
        # expr=
        m.fs.costing.factor_capital_annualization *(
        pyo.value(m.fs.CL.costing.direct_capital_cost) +
        pyo.value(m.fs.R1.costing.direct_capital_cost) +
        pyo.value(m.fs.R2.costing.direct_capital_cost) +
        pyo.value(m.fs.R3.costing.direct_capital_cost) +
        pyo.value(m.fs.R4.costing.direct_capital_cost) +
        pyo.value(m.fs.R5.costing.direct_capital_cost) +
        pyo.value(m.fs.CL1.costing.direct_capital_cost) +
        pyo.value(m.fs.RADM.costing.direct_capital_cost) +
        pyo.value(m.fs.DU.costing.direct_capital_cost) +
        pyo.value(m.fs.TU.costing.direct_capital_cost) )
        / m.fs.costing.annual_feedwater_volume
    ))
    print("electricity", pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
        * m.fs.costing.utilization_factor
        / m.fs.costing.annual_feedwater_volume
    ))
    print("electricity", pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
    ))

    print("opex",pyo.value(m.fs.costing.maintenance_labor_chemical_operating_cost
        / m.fs.costing.annual_feedwater_volume
    ))

    print("fixed",(pyo.value((m.fs.costing.aggregate_fixed_operating_cost
            + m.fs.costing.maintenance_labor_chemical_operating_cost )/ m.fs.costing.annual_feedwater_volume)))
    print("variable",(pyo.value((m.fs.costing.aggregate_variable_operating_cost
                + sum(m.fs.costing.aggregate_flow_costs[f] for f in m.fs.costing.used_flows)
                * m.fs.costing.utilization_factor)/ m.fs.costing.annual_feedwater_volume)))

    unit2=np.zeros((2,5))
    rows = 2
    columns = 5
    for j in range(rows):
        for i in range(columns):
            unit2[j][i] = unit1[j][i]

    print(unit2)
    # y2=np.zeros((2,5))
    # rows = 2
    # columns = 5
    # for j in range(rows):
    #     for i in range(columns):
    #         y2[j][i] = y1[j][i] * 0.07 /(inlet[j]/24)

    # print(y2)

    print(op1)
    op2=np.zeros((2,2))
    rows = 2
    columns = 2
    for j in range(rows):
        for i in range(columns):
            op2[j][i] = op1[j][i]

    # print(y2)

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

    groups = ["Simulation", "Optimization"]
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
    ax.legend(fontsize=20)
    plt.savefig("cost_sim_vs_opt.png", bbox_inches="tight", dpi=300)
    plt.show() 

    
    # sections = ['Primary', 'Secondary', 'Digestion', 'dewatering', 'thickening']
    # barWidth = 0.25
    # br1 = np.arange(len(reactors)) 
    # br2 = [x + barWidth for x in br1] 
    # br3 = [x + barWidth for x in br2]

    # fig, ax = plt.subplots()
    # ax.bar(br1, V1_1, color = "navy", width = barWidth, label ='Simulation')
    # ax.bar(br2, V1_2, color = "dodgerblue", width = barWidth, label ='Optimization')

    # # plt.xticks(x)
    # plt.ylabel('Volume m$^3$')
    # plt.title('Simulation vs. Optimization')
    # plt.xticks([r + barWidth for r in range(len(reactors))], 
    #     ["Reactor 1", "Reactor 2", "Reactor 3", "Reactor 4", "Reactor 5"])
    # ax.legend()
    # plt.show() 