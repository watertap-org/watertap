import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyomo.environ import (
    Var,
    Constraint,
    ConcreteModel,
    units as pyunits,
    value)

from pyomo.util.check_units import assert_units_consistent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
import watertap.property_models.cryst_prop_pack as props_cryst


def main():
    save_dir = "C:/Users/carso/Documents/MVC/watertap_results/Property comparisons/"
    save_file = "Seawater enthalpy.csv"
    prop_pack = {}
    prop_pack['Parameters'] = props_sw.SeawaterParameterBlock()
    prop_pack['solute'] = 'TDS'
    t = np.arange(30,100,5)
    mass_fractions = [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2]

    wf_rr = [(0.05, 0.5), (0.1, 0.5), (0.05, 0.7), (0.1, 0.7)]
    save_dir = "C:/Users/carso/Documents/MVC/watertap_results/Property comparisons/energy_balance_comparison/"
    # save_q_evap_change(t,wf_rr,save_dir,prop_pack)


    prop_pack['Parameters'] = props_cryst.NaClParameterBlock()
    prop_pack['solute'] = 'NaCl'
    # save_q_evap_change(t,wf_rr,save_dir,prop_pack)
    plot_energy_balance_difference(t,wf_rr,save_dir)

def save_q_evap_change(t,wf_rr,save_dir, prop_pack, t_f_preheated=None):
    """
    t: list of temperatures [C]
    wf_rr: list of (wf, rr) tuples
    save_dir: directory for saving results
    prop_pack: dictionary pointing to parameter block and name of solute
    t_f_preheated: if None, then T_f = T_b. Otherwise T_f is specified value
    saves file of Q_evap, m_f, m_b, m_v, h_f, h_b, h_v, H_f, H_b, H_v
    """

    T = t + 273.15
    solute = prop_pack['solute']

    # model
    m = build_evap_energy_balance_model(prop_pack)
    # initialize
    initialize_evap_model(m,wf_rr[0][0], wf_rr[0][1], solute)

    # iterate over wf, rr
    solver = get_solver()
    optarg = solver.options

    for (wf,rr) in wf_rr:
        data = {}
        data['temperature'] = t
        Q_evap = []
        m_f = []
        m_b = []
        m_v = []
        h_f = []
        h_b = []
        h_v = []
        for temp in T:
            m.fs.feed.properties[0].mass_frac_phase_comp['Liq', solute].fix(wf)
            m.fs.recovery.fix(rr)
            m.fs.feed.properties[0].temperature.fix(temp)
            m.fs.brine.properties[0].temperature.fix(temp)
            results = solver.solve(m, tee=False)
            Q_evap.append(m.fs.heat_transfer.value)
            m_f.append(m.fs.feed.properties[0].flow_mass_phase_comp['Liq',"H2O"].value + m.fs.feed.properties[0].flow_mass_phase_comp['Liq',solute].value)
            m_b.append(m.fs.brine.properties[0].flow_mass_phase_comp['Liq',"H2O"].value + m.fs.brine.properties[0].flow_mass_phase_comp['Liq',solute].value)
            m_v.append(m.fs.vapor.properties[0].flow_mass_phase_comp['Vap',"H2O"].value)
            h_f.append(m.fs.feed.properties[0].enth_mass_phase['Liq'].value)
            h_b.append(m.fs.brine.properties[0].enth_mass_phase['Liq'].value)
            h_v.append(m.fs.vapor.properties[0].enth_mass_phase['Vap'].value)

        # Save results
        data['Q_evap'] = Q_evap
        data['m_f'] = m_f
        data['m_b'] = m_b
        data['m_v'] = m_v
        data['h_f'] = h_f
        data['h_b'] = h_b
        data['h_v'] = h_v
        # Save to csv
        df = pd.DataFrame(data=data)
        filename = solute + '_wf_' + str(wf*1e3) + '_rr_' + str(rr*100) + ".csv"
        df.to_csv(save_dir + filename, index=False)

def build_evap_energy_balance_model(prop_pack):
    solute = prop_pack['solute']
    # Build model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Build feed state block
    m.fs.properties = prop_pack['Parameters']
    m.fs.feed = Feed(default={"property_package": m.fs.properties})

    # Build brine state block
    m.fs.brine = Feed(default={"property_package": m.fs.properties})

    # Build vapor state block
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.vapor = Feed(default={"property_package": m.fs.properties_vapor})

    # Add recovery ratio
    m.fs.recovery = Var(m.fs.config.time, initialize=0.5, bounds=(0, 1))
    m.fs.recovery_equation = Constraint(expr=m.fs.vapor.properties[0].flow_mass_phase_comp["Vap", "H2O"] ==
                                             m.fs.recovery[0] * (
                                                         m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"] +
                                                         m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute]))

    # Add water mass balance
    @m.fs.Constraint()
    def eq_water_mass_balance(blk):
        return (
                blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"] ==
                blk.brine.properties[0].flow_mass_phase_comp["Liq", "H2O"] +
                blk.vapor.properties[0].flow_mass_phase_comp['Vap', "H2O"]
        )
    # Add solute mass balance
    @m.fs.Constraint()
    def eq_solute_mass_balance(blk):
        return (
                blk.feed.properties[0].flow_mass_phase_comp["Liq", solute] ==
                blk.brine.properties[0].flow_mass_phase_comp["Liq", solute]
        )

    # Add energy balance
    m.fs.heat_transfer = Var(initialize=1e4, bounds=(1, 1e10), units=pyunits.J * pyunits.s**-1)

    @m.fs.Constraint()
    def eq_energy_balance(blk):
        return (
                blk.heat_transfer + blk.feed.properties[0].enth_flow
                == blk.brine.properties[0].enth_flow
                + blk.vapor.properties[0].enth_flow_phase["Vap"]
        )

    # Brine pressure constraint
    @m.fs.Constraint()
    def eq_brine_pressure(blk):
        return blk.brine.properties[0].pressure == blk.brine.properties[0].pressure_sat

    # Vapor pressure
    @m.fs.Constraint()
    def eq_vapor_pressure(blk):
        return blk.vapor.properties[0].pressure == blk.brine.properties[0].pressure

    # Vapor temperature
    @m.fs.Constraint()
    def eq_vapor_temperature(blk):
        return (
                blk.vapor.properties[0].temperature == blk.brine.properties[0].temperature
        )

    # Scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", solute)
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )


    return m

def initialize_evap_model(m, wf, rr, solute):
    # set operating conditions for initialization
    # feed
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', solute].fix(wf)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    # recovery
    m.fs.recovery.fix(rr)

    # brine
    m.fs.brine.properties[0].temperature.fix(273.15 + 25)
    if solute == 'NaCl':
        m.fs.feed.properties[0].flow_mass_phase_comp['Sol','NaCl'].fix(0)
        m.fs.feed.properties[0].flow_mass_phase_comp['Vap', 'H2O'].fix(0)
        m.fs.brine.properties[0].flow_mass_phase_comp['Sol','NaCl'].fix(0)
        m.fs.brine.properties[0].flow_mass_phase_comp['Vap', 'H2O'].fix(0)

    # vapor
    m.fs.vapor.properties[0].flow_mass_phase_comp['Liq', "H2O"].fix(0)

    # initialize
    solver = get_solver()
    optarg = solver.options

    print('DOF: ', degrees_of_freedom(m.fs))
    results = solver.solve(m, tee=False)
    print(results.solver.termination_condition)

    if results.solver.termination_condition != 'optimal':
        print('Initialization not optimal')
        assert False

def plot_energy_balance(t, wf_rr, save_dir):
    colors = plt.cm.gnuplot2(np.linspace(0, 1, len(wf_rr) + 1))
    i = 0
    plt.figure()
    for (wf, rr) in wf_rr:
        filename = save_dir + 'TDS_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_tds = pd.read_csv(filename)
        filename = save_dir + 'NaCl_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_nacl = pd.read_csv(filename)
        plt.plot(t,df_tds['Q_evap'], '.',color=colors[i], label = 'Seawater ' + str(wf*1000) + 'g/kg, ' + str(rr*100) + '%')
        plt.plot(t,df_nacl['Q_evap'], color=colors[i], label = 'NaCl Cryst ' + str(wf*1000) + 'g/kg, ' + str(rr*100) + '%')
        i +=1
    plt.xlabel('Temperature [C]')
    plt.ylabel('Heat transfer [J/s]')
    plt.legend()
    plt.show()

    i = 0
    plt.figure()
    for (wf, rr) in wf_rr:
        filename = save_dir + 'TDS_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_tds = pd.read_csv(filename)
        m_v = df_tds['m_v'][0]
        filename = save_dir + 'NaCl_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_nacl = pd.read_csv(filename)
        plt.plot(t,df_tds['Q_evap']/m_v, '.',color=colors[i], label = 'Seawater ' + str(wf*1000) + 'g/kg, ' + str(rr*100) + '%')
        plt.plot(t,df_nacl['Q_evap']/m_v, color=colors[i], label = 'NaCl Cryst ' + str(wf*1000) + 'g/kg, ' + str(rr*100) + '%')
        i +=1
    plt.xlabel('Temperature [C]')
    plt.ylabel('Heat transfer [J/kg-vapor]')
    plt.legend()
    plt.show()

def plot_energy_balance_difference(t, wf_rr, save_dir):
    colors = plt.cm.gnuplot2(np.linspace(0, 1, len(wf_rr) + 1))
    i = 0
    plt.figure()
    for (wf, rr) in wf_rr:
        filename = save_dir + 'TDS_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_tds = pd.read_csv(filename)
        filename = save_dir + 'NaCl_wf_' + str(wf*1000) + "_rr_"+ str(rr*100)+".csv"
        df_nacl = pd.read_csv(filename)
        Q_dif = df_nacl['Q_evap']-df_tds['Q_evap']
        Q_dif_per = Q_dif.div(df_tds['Q_evap'])*100
        plt.plot(t,Q_dif_per, color=colors[i], label = str(wf*1000) + 'g/kg, ' + str(rr*100) + '%')
        i +=1
    plt.xlabel('Temperature [C]')
    plt.ylabel(r'($Q_{NaCl} - Q_{seawater}$)/$Q_{seawater} \times 100$ [%]')
    plt.legend()
    plt.show()

    # i = 0
    # plt.figure()
    # for (wf, rr) in wf_rr:
    #     filename = save_dir + 'TDS_wf_' + str(wf * 1000) + "_rr_" + str(rr * 100) + ".csv"
    #     df_tds = pd.read_csv(filename)
    #     filename = save_dir + 'NaCl_wf_' + str(wf * 1000) + "_rr_" + str(rr * 100) + ".csv"
    #     df_nacl = pd.read_csv(filename)
    #     m_v = df_tds['m_v'][0]
    #     Q_dif = df_nacl['Q_evap']/m_v - df_tds['Q_evap']/m_v
    #     Q_dif_per = Q_dif.div(df_tds['Q_evap']/m_v)
    #     plt.plot(t, Q_dif_per, color=colors[i], label=str(wf * 1000) + 'g/kg, ' + str(rr * 100) + '%')
    #     i += 1
    # plt.xlabel('Temperature [C]')
    # # plt.ylabel(r'$(Q_{NaCl} - Q_{seawater})/m_{vapor}$ [J/kg-vapor]')
    # plt.ylabel(r'$(Q_{NaCl} - Q_{seawater})/m_{vapor}$ [J/kg-vapor]')
    # plt.legend()
    # plt.show()

def molality_to_mass_fraction_NaCl(m):
    mol_weight = 1e-3*58.4428 # kg/mol for NaCl
    return m*mol_weight/(1+m*mol_weight)

def save_h_wT(t, mass_fractions, save_dir, filename, prop=None, prop_pack=None):
    """
    t: temperatures to be evaluated (C)
    w: mass fractions to be evaluated (kg/kg)
    save_dir: directory for saving results
    filename: filename for saving results
    prop: the property to save
    prop_pack: dictionary pointing to parameter block and name of solute
    """

    T = t + 273.15 # K
    solute = prop_pack['solute']

    # Build state block
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = prop_pack['Parameters']
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', solute].fix(mass_fractions[0])
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    # Scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", solute)
    )

    # touch properties
    m.fs.feed.properties[0].enth_mass_phase['Liq']
    m.fs.feed.properties[0].pressure_sat

    # initialize
    solver = get_solver()
    optarg = solver.options

    print('DOF: ', degrees_of_freedom(m.fs.feed))
    results = solver.solve(m, tee=False)
    print(results.solver.termination_condition)

    data = {}
    data['temperature'] = t
    for w in mass_fractions:
        m.fs.feed.properties[0].mass_frac_phase_comp['Liq', solute].fix(w)
        h = []
        for temp in T:
            m.fs.feed.properties[0].temperature.fix(temp)
            results = solver.solve(m, tee=False)
            if results.solver.termination_condition == 'infeasible':
                print('Infeasible for w = ', str(w), 'kg/kg, T = ', str(temp), ' K')
            if prop == 'enthalpy':
                h.append(m.fs.feed.properties[0].enth_mass_phase['Liq'].value*1e3)
            elif prop == 'vapor pressure':
                h.append(m.fs.feed.properties[0].pressure_sat.value)
        data[w] = h

    # save to csv
    df = pd.DataFrame(data=data)
    df.to_csv(save_dir+filename,index=False)

def plot_enthalpy_comparison(file_1, file_2,t, mass_fractions):
    # read in files
    df_1 = pd.read_csv(file_1['name'])
    df_2 = pd.read_csv(file_2['name'])
    # reference if enthalpies
    df_1 = df_1 - df_1['0'][0]
    df_2 = df_2 - df_2['0'][0]
    df_1 = df_1*1e-3
    df_2 = df_2*1e-3

    df_1.drop('temperature', axis=1, inplace=True)
    df_2.drop('temperature', axis=1, inplace=True)

    df_dif = df_2-df_1
    df_per_dif = df_dif.div(df_1)*100

    plt.figure()
    # plt.plot([120,120],[-100,900],'r--') # temperature limit
    colors = plt.cm.gnuplot2(np.linspace(0, 1, len(mass_fractions)+1))
    i = 0
    for col_name in df_per_dif.columns:
        label = 'w = ' + str(round(float(col_name), 3))
        plt.plot(t,df_per_dif[col_name],color=colors[i],label=label)
        i += 1
    plt.xlabel('Temperature')
    #plt.ylabel(r'$(h_{cryst NaCl}-h_{seawater})$ [kJ/kg]')

    plt.ylabel(r'$(h_{cryst NaCl}-h_{seawater})$/$h_{seawater}\times 100$ [%]')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_vapor_pressure_comparison(file_1, file_2,t, mass_fractions):
    # read in files
    df_1 = pd.read_csv(file_1['name'])
    df_2 = pd.read_csv(file_2['name'])
    df_1 = df_1*1e-3
    df_2 = df_2*1e-3

    df_1.drop('temperature', axis=1, inplace=True)
    df_2.drop('temperature', axis=1, inplace=True)



    plt.figure()
    colors = plt.cm.gnuplot2(np.linspace(0, 1, len(mass_fractions)+1))
    i = 0
    for col_name in df_1.columns:
        label = 'Seawater w = ' + str(round(float(col_name), 3))
        plt.plot(t,df_1[col_name],color=colors[i],label=label)
        label = 'NaCl w = ' + str(round(float(col_name), 3))
        plt.plot(t,df_2[col_name],'.',color=colors[i],label=label)
        i += 1
    plt.xlim(25,100)
    plt.ylim(0,100)
    plt.xlabel('Temperature')
    plt.ylabel(r'$P^{sat}$ [kPa]')
    plt.legend(loc='center left')
    plt.tight_layout()
    plt.show()

    df_dif = df_2-df_1
    df_per_dif = df_dif.div(df_1)*100
    plt.figure()
    i = 0
    for col_name in df_dif.columns:
        label = 'w = ' + str(round(float(col_name), 3))
        plt.plot(t,df_dif[col_name],color=colors[i],label=label)
        i += 1
    plt.xlabel('Temperature')
    plt.ylabel(r'$(P^{sat}_{cryst NaCl}-P^{sat}_{seawater})$ [kPa]')
    #plt.ylabel(r'$(P^{sat}_{cryst NaCl}-P^{sat}_{seawater})$/$P^{sat}_{seawater}\times 100$ [%]')
    plt.legend()
    plt.tight_layout()
    plt.show()

def seawater_properties(t, mass_fractions,colors):
    m = ConcreteModel()
    T = t + 273.15

    # Build state block
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties_seawater = props_sw.SeawaterParameterBlock()
    m.fs.feed = Feed(default={"property_package": m.fs.properties_seawater})
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].fix(mass_fractions[0])
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    # Scaling
    m.fs.properties_seawater.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_seawater.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    # touch properties
    m.fs.feed.properties[0].enth_mass_phase['Liq']

    # initialize
    solver = get_solver()
    optarg = solver.options

    print('DOF: ', degrees_of_freedom(m.fs.feed))
    results = solver.solve(m, tee=False)
    print(results.solver.termination_condition)
    h_ref = m.fs.feed.properties[0].enth_mass_phase['Liq'].value*1e-3
    i = 0
    for w in mass_fractions:
        m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS'].fix(w)
        h = []
        for temp in T:
            m.fs.feed.properties[0].temperature.fix(temp)
            results = solver.solve(m, tee=False)
            if results.solver.termination_condition == 'infeasible':
                print('Infeasible for w = ', str(w), 'kg/kg, T = ', str(temp), ' K')
            #h.append(m.fs.feed.properties[0].enth_mass_phase['Liq'].value*1e-3-h_ref)
            h.append()
        plt.plot(t,h,'.',color = colors[i], label='Seawater w = ' + str(round(w,3)))
        i += 1

    return

def cryst_properties(t, mass_fractions,colors):
    m = ConcreteModel()
    T = t + 273.15

    # Build state block
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties_cryst = props_cryst.NaClParameterBlock()
    m.fs.feed = Feed(default={"property_package": m.fs.properties_cryst})
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(mass_fractions[0])
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    # Scaling
    m.fs.properties_cryst.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_cryst.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    # touch properties
    m.fs.feed.properties[0].enth_mass_phase['Liq']

    # initialize
    solver = get_solver()
    optarg = solver.options

    print('DOF: ', degrees_of_freedom(m.fs.feed))
    results = solver.solve(m, tee=False)
    print(results.solver.termination_condition)
    h_ref = m.fs.feed.properties[0].enth_mass_phase['Liq'].value
    print(h_ref)
    i = 0
    for w in mass_fractions:
        m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(w)
        h = []
        for temp in T:
            m.fs.feed.properties[0].temperature.fix(temp)
            results = solver.solve(m, tee=False)
            if results.solver.termination_condition == 'infeasible':
                print('Infeasible for w = ', str(w), 'kg/kg, T = ', str(temp), ' K')
            h.append(m.fs.feed.properties[0].enth_mass_phase['Liq'].value-h_ref)
        plt.plot(t,h,'x',color = colors[i], label='Cryst NaCl w = ' + str(round(w,3)))
        i += 1

    return

def compare_pitzer():
    pitzer_h_file = "C:/Users/carso/Documents/MVC/watertap_results/Property comparisons/Pitzer Specific Enthalpy.csv"
    df_pitzer_h = pd.read_csv(pitzer_h_file)
    n = df_pitzer_h.shape[0]
    # remove rows with temperature > 100
    df_pitzer_h.drop(labels=range(n - 12, n), axis=0, inplace=True)
    t = np.array(df_pitzer_h[:]['t'])
    P = np.array(df_pitzer_h[:]['P'])

    df_pitzer_h.drop('t', axis=1, inplace=True)
    df_pitzer_h.drop('P', axis=1, inplace=True)
    df_pitzer_h.drop('0.25', axis=1, inplace=True)
    df_pitzer_h.drop('0.5', axis=1, inplace=True)
    df_pitzer_h.drop('0.75', axis=1, inplace=True)
    df_pitzer_h.drop('4', axis=1, inplace=True)
    df_pitzer_h.drop('5', axis=1, inplace=True)

    df_pitzer_h = (df_pitzer_h - df_pitzer_h['0.1'][2]) * 1000
    print(df_pitzer_h)
    mass_fractions = []
    plt.figure()
    # plt.plot([120,120],[-100,900],'r--') # temperature limit
    colors = plt.cm.gnuplot2(np.linspace(0, 1, 6))
    i = 0
    for col_name in df_pitzer_h.columns:
        w = molality_to_mass_fraction_NaCl(float(col_name))
        mass_fractions.append(w)
        label = 'NaCl w = ' + str(round(w, 3))
        # plt.plot(t,df_pitzer_h[col_name],color=colors[i],label=label)
        i += 1
    plt.xlabel('Temperature')
    plt.ylabel('Specific enthalpy [J/kg]')
    seawater_properties(t, mass_fractions, colors)
    cryst_properties(t, mass_fractions, colors)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()