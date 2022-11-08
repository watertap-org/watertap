import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import seaborn as sns
import numpy as np
import math

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    value,
    Expression,
    Param,
    Objective
)

from idaes.core.solvers import get_solver
from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
from watertap.examples.flowsheets.mvc import mvc_single_stage as mvc_full
from watertap.examples.flowsheets.mvc import evap_comp_cond as mvc_unit
from watertap.examples.flowsheets.mvc import mvc_single_stage_distillate_preheater as mvc_distillate_hx_only
import mvc_plotting as mvc_plot


def main():

    # dir = "C:/Users/carso/Documents/MVC/watertap_results/dual_c_evap_U_evap_sensitivity"
    # map_dir = dir + "/plus_20_minus_20"
    # analysis = dir + "/analysis_dual_c_evap_U_evap.csv"
    # output_file = map_dir + "/optimize_sweep.csv"
    # global_results, sweep_params, m = run_multi_param_case(analysis, system='mvc_full_opt', output_filename=output_file,f_evap=3000, f_hx=2000, T_cv_max=450, C_elec=0.15)
    # save_results_for_plotting(output_file, map_dir, 9, 9)

    analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_full_optimize_cases.csv"
    dir = "C:/Users/carso/Documents/MVC/watertap_results/tornado_sensitivity/"

    n_param = 2

    # # Evaporator cost
    # param = 'evaporator_cost'
    # param_min = 2000
    # param_max = 4000
    # run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')
    #
    # # electicity cost
    param = 'electricity_cost'
    param_min = 0.075
    param_max = 0.3
    run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')
    assert False

    # evaporator overall heat transfer coefficient
    # param = 'U_evap'
    # param_min = 800
    # param_max = 1200
    # run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')
    # assert False
    # heat exchanger cost
    param = 'preheater_cost'
    param_min = 1500
    param_max = 2500
    run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')

    # Compressor cost
    # param = 'compressor_cost'
    # param_min = 0.5*7364
    # param_max = 1.5*7364
    # run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')

    # Compressor efficiency
    # param = 'compressor_efficiency'
    # param_min = 0.7
    # param_max = 0.9
    # run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')

    # Distillate heat exchanger overall heat transfer coefficient
    param = 'U_hx_distillate'
    param_min = 800
    param_max = 1200
    run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')

    # Brine heat exchanger overall heat transfer coefficient
    param = 'U_hx_brine'
    run_case(n_param, param, param_min, param_max, system='mvc_full_opt', output_filename=dir+param+'.csv')

    assert False
    cases = {}
    cases['evap_hx_cost'] = [(3000,2000)]
    cases['elec_cost'] = [0.15]
    cases['cv_temp_max'] = [450]
    cases['comp_cost_factor'] = [1]

    map_dir = "C:/Users/carso/Documents/MVC/watertap_results/full_parameter_sweeps_P_out_unfixed"
    map_dir = map_dir + "/evap_3000_hx_2000/elec_0.15/cv_temp_max_450/comp_cost_1"
    # run_full_parameter_sweeps(analysis,cases,map_dir, system='mvc_distillate_hx_only')
    convert_units_results(map_dir)
    assert False
    save_dir = map_dir + '/figures'
    mvc_plot.make_maps(map_dir, save_dir)

    # dir = "C:/Users/carso/Documents/MVC/watertap_results/dual_c_evap_c_comp_sensitivity"
    # map_dir = dir
    # analysis = dir + "/analysis_dual_c_evap_c_comp.csv"
    # output_file = map_dir + "/optimize_sweep.csv"
    # global_results, sweep_params, m = run_multi_param_case(analysis, system='mvc_full_opt', output_filename=output_file,f_evap=3000, f_hx=2000, T_cv_max=450, C_elec=0.15)
    # save_results_for_plotting(output_file, map_dir, 9, 9)

def mvc_unit_presweep():
    m = mvc_unit.build()
    mvc_unit.set_operating_conditions(m)
    mvc_unit.initialize_system(m)
    mvc_unit.scale_costs(m)
    solver = get_solver()
    results = solver.solve(m, tee=False)
    return m

def mvc_full_presweep(f_evap=3000, f_hx=2000, T_cv_max=450, T_b=None, C_elec=0.15, C_evap_hx_ratio=None, C_comp_factor=1):
    m = mvc_full.build()
    mvc_full.add_Q_ext(m, time_point=m.fs.config.time)
    mvc_full.set_operating_conditions(m)

    # Fix values to desired values
    m.fs.costing.heat_exchanger_unit_cost.fix(f_hx)
    if C_evap_hx_ratio is not None:
        print('Fixing preheater cost based on given evaporator cost and C_evap_hx_ratio')
        m.fs.costing.heat_exchanger_unit_cost.fix(f_evap/2)
    m.fs.costing.evaporator_unit_cost.fix(f_evap)
    m.fs.costing.electricity_base_cost = C_elec
    compressor_cost = m.fs.costing.compressor_unit_cost.value
    m.fs.costing.compressor_unit_cost.fix(compressor_cost*C_comp_factor)

    # Initialize
    mvc_full.initialize_system(m)
    mvc_full.scale_costs(m)
    mvc_full.fix_outlet_pressures(m)
    # solve
    # solver = get_solver()
    # m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    # results = solver.solve(m, tee=False)
    results = mvc_full.sweep_solve(m)
    mvc_full.display_results(m)

    m.fs.compressor.control_volume.properties_out[0].temperature.setub(T_cv_max)

    if T_b is not None:
        m.fs.evaporator.properties_brine[0].temperature.fix(T_b)
    return m

def mvc_distillate_hx_only_presweep(f_evap=1000, f_hx=1000, T_cv_max=450, T_b=None, C_elec=0.07, C_evap_hx_ratio=None, C_comp_factor=1):
    m = mvc_distillate_hx_only.build()
    mvc_distillate_hx_only.add_Q_ext(m, time_point=m.fs.config.time)
    mvc_distillate_hx_only.set_operating_conditions(m)

    # Fix values to desired values
    m.fs.costing.heat_exchanger_unit_cost.fix(f_hx)
    if C_evap_hx_ratio is not None:
        print('Fixing preheater cost based on given evaporator cost and C_evap_hx_ratio')
        m.fs.costing.heat_exchanger_unit_cost.fix(f_evap/2)
    m.fs.costing.evaporator_unit_cost.fix(f_evap)
    m.fs.costing.electricity_base_cost = C_elec
    compressor_cost = m.fs.costing.compressor_unit_cost.value
    m.fs.costing.compressor_unit_cost.fix(compressor_cost*C_comp_factor)

    # Initialize
    mvc_distillate_hx_only.initialize_system(m)
    mvc_distillate_hx_only.scale_costs(m)
    mvc_distillate_hx_only.fix_outlet_pressures(m)
    # solve
    solver = get_solver()
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    results = solver.solve(m, tee=False)
    results = mvc_distillate_hx_only.sweep_solve(m)
    mvc_distillate_hx_only.display_results(m)
    m.fs.compressor.control_volume.properties_out[0].temperature.setub(T_cv_max)
    if T_b is not None:
        m.fs.evaporator.properties_brine[0].temperature.fix(T_b)
    return m

def run_full_parameter_sweeps(analysis_file, cases, dir, system='mvc_full_opt'):

    for evap_hx_cost in cases['evap_hx_cost']:
        print('Evaporator, hx cost: ', evap_hx_cost)
        for elec_cost in cases['elec_cost']:
            print('Electricity cost:', elec_cost)
            for cv_temp_max in cases['cv_temp_max']:
                print('Compressed vapor temperature:', cv_temp_max)
                for comp_cost_factor in cases['comp_cost_factor']:
                    print('Compressor cost:', comp_cost_factor)
                    map_dir = dir + '/evap_'+ str(evap_hx_cost[0]) + '_hx_' + str(evap_hx_cost[1]) + '/elec_' + str(elec_cost) + '/cv_temp_max_' + str(cv_temp_max) + '/comp_cost_' + str(comp_cost_factor)
                    output_file = map_dir + '/sweep_results.csv'
                    global_results, sweep_params, m = run_multi_param_case(analysis_file=analysis_file,
                                                                           system=system,
                                                                           output_filename=output_file,
                                                                           f_evap=evap_hx_cost[0],
                                                                           f_hx=evap_hx_cost[1],
                                                                           T_cv_max=cv_temp_max,
                                                                           C_elec=elec_cost,
                                                                           C_comp_factor=comp_cost_factor,
                                                                           )
                    save_dir = map_dir +'/figures'
                    save_results_for_plotting(output_file,map_dir,7,9)
                    convert_units_results(map_dir)
                    mvc_plot.make_maps(map_dir, save_dir)

def run_analysis(analysis_file, fixed_params):
    df = pd.read_csv(analysis_file)
    row = 0
    for param in df["Parameter"]:
        print(param)
        if param in fixed_params:
            filename = "C:/Users/carso/Documents/MVC/watertap_results/type1_" + param + "_results.csv"
            global_results, sweep_params, m = run_case(df['N'][row],
                     param=param,
                     param_min=df['Min'][row],
                     param_max=df['Max'][row],
                     system='mvc_unit',
                     output_filename=filename)
            print(global_results)
            assert False
        row += 1

def run_multi_param_case(analysis_file, system='mvc_unit',output_filename=None,f_evap=1000,f_hx=1000, T_cv_max=500, T_b=None, C_elec=0.07, C_evap_hx_ratio=None, C_comp_factor=1):
    df = pd.read_csv(analysis_file)
    if output_filename is None:
        output_filename = ("C:/Users/carso/Documents/MVC/watertap_results/" + system + "_multi_param_results.csv")

    if system == 'mvc_unit':
        # Get model
        m = mvc_unit_presweep()
        opt_fcn = mvc_unit.solve
        outputs = make_outputs_dict_mvc_unit(m)

    elif system == 'mvc_full':
        # Get model
        m = mvc_full_presweep(f_evap=f_evap,f_hx=f_hx, T_cv_max=T_cv_max, T_b=T_b, C_elec=C_elec, C_evap_hx_ratio=C_evap_hx_ratio)
        opt_fcn = mvc_full.solve
        outputs = make_outputs_dict_mvc_full(m)

    elif system == 'mvc_full_opt':
        m = mvc_full_presweep(f_evap=f_evap,f_hx=f_hx, T_cv_max=T_cv_max, T_b=T_b, C_elec=C_elec,C_evap_hx_ratio=C_evap_hx_ratio,C_comp_factor=C_comp_factor)
        if T_b is not None:
            print(T_b)
            opt_fcn = mvc_full.sweep_solve_fixed_brine_temp
        else:
            opt_fcn = mvc_full.sweep_solve
        outputs= make_outputs_dict_mvc_full(m)

    elif system == 'mvc_distillate_hx_only':
        m = mvc_distillate_hx_only_presweep(f_evap=f_evap,f_hx=f_hx, T_cv_max=T_cv_max, T_b=T_b, C_elec=C_elec,C_evap_hx_ratio=C_evap_hx_ratio,C_comp_factor=C_comp_factor)
        opt_fcn = mvc_distillate_hx_only.sweep_solve
        outputs = make_outputs_dict_mvc_distillate_hx_only(m)

    else:
        print('Flowsheet not correctly specified')
        assert False

    # Sweep parameter
    sweep_params = {}
    param_vars = get_param_var_dict(m)
    row = 0
    for param in df["Parameter"]:
        if df['N'][row] > 0:
            print('adding ', param)
            sweep_params[param] = LinearSample(param_vars[param], df['Min'][row], df['Max'][row], df['N'][row])
        row +=1

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_fcn,
        interpolate_nan_outputs=False,
    )

    return global_results, sweep_params, m

def run_case(n_param, param=None, param_min=None, param_max=None, system='mvc_unit', output_filename=None, T_b=None, T_cv_max=450):
    """
    Run the parameter sweep tool on MVC flowsheet, sweeping over vapor flow rate from 4.5 to 5.5 kg/s

    Arguments
    ---------
    n_param (int): number of points to run
    output_filename (str, optional): the place to write the parameter sweep results csv file

    Returns
    -------

    """
    if output_filename is None:
        output_filename = ("C:/Users/carso/Documents/MVC/watertap_results/results.csv")

    if system == 'mvc_unit':
        # Get model
        m = mvc_unit_presweep()
        opt_fcn = mvc_unit.solve
        outputs = make_outputs_dict_mvc_unit(m)

    elif system == 'mvc_full':
        # Get model
        m = mvc_full_presweep()
        opt_fcn = mvc_full.solve
        outputs = make_outputs_dict_mvc_full(m)

    elif system == 'mvc_full_opt':
        m = mvc_full_presweep(T_cv_max=T_cv_max)
        if T_b is not None:
            print(T_b)
            opt_fcn = mvc_full.sweep_solve_fixed_brine_temp
        else:
            opt_fcn = mvc_full.sweep_solve
        outputs= make_outputs_dict_mvc_full(m)

    else:
        print('Flowsheet not correctly specified')
        assert False

    # Sweep parameter
    sweep_params = {}
    param_vars = get_param_var_dict(m)
    sweep_params[param] = LinearSample(param_vars[param], param_min, param_max, n_param)

    # sweep
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_fcn,
        interpolate_nan_outputs=False,
    )

    return global_results, sweep_params, m

def make_outputs_dict_mvc_unit(m):
    outputs = {}
    # Feed
    outputs['Feed mass flow water'] = m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']
    outputs['Feed mass flow salt'] = m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'TDS']
    outputs['Feed mass fraction'] = m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS']
    outputs['Feed temperature'] = m.fs.feed.properties[0].temperature
    outputs['Feed pressure'] = m.fs.feed.properties[0].pressure

    # Brine from evaporator
    outputs['Brine mass flow water'] = m.fs.brine.properties[0].flow_mass_phase_comp['Liq', 'H2O']
    outputs['Brine mass flow salt'] = m.fs.brine.properties[0].flow_mass_phase_comp['Liq', 'TDS']
    outputs['Brine temperature'] = m.fs.evaporator.properties_brine[0].temperature
    outputs['Brine pressure'] = m.fs.evaporator.properties_brine[0].pressure

    # Vapor
    outputs['Vapor mass flow'] = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap','H2O']
    outputs['Vapor temperature'] = m.fs.evaporator.properties_vapor[0].temperature
    outputs['Vapor pressure'] = m.fs.evaporator.properties_vapor[0].pressure

    # Compressed vapor
    outputs['Compressed vapor temperature'] = m.fs.compressor.control_volume.properties_out[0].temperature
    outputs['Compressed vapor pressure'] = m.fs.compressor.control_volume.properties_out[0].pressure

    # Condensed vapor/distillate
    outputs['Distillate temperature'] = m.fs.condenser.control_volume.properties_out[0].temperature
    outputs['Distillate pressure'] = m.fs.condenser.control_volume.properties_out[0].pressure

    # Exiting distillate
    outputs['Exiting distillate pressure'] = m.fs.distillate.properties[0].pressure

    # Exiting brine
    outputs['Exiting brine pressure'] = m.fs.brine.properties[0].pressure

    # Evaporator performance
    outputs['Evaporator area'] = m.fs.evaporator.area
    outputs['Evaporator LMTD'] = m.fs.evaporator.lmtd
    outputs['Evaporator heat transfer'] = m.fs.evaporator.heat_transfer
    outputs['Evaporator overall heat transfer coefficient'] = m.fs.evaporator.U
    outputs['Evaporator approach temperature in'] = m.fs.evaporator.delta_temperature_in
    outputs['Evaporator approach temperature out'] = m.fs.evaporator.delta_temperature_out

    # Compressor performance
    outputs['Compressor pressure ratio'] = m.fs.compressor.pressure_ratio
    outputs['Compressor work'] = m.fs.compressor.control_volume.work[0]
    outputs['Compressor efficiency'] = m.fs.compressor.efficiency

    # Costing/outcome metrics
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['SEC'] = m.fs.costing.specific_energy_consumption

    return outputs

def make_outputs_dict_mvc_distillate_hx_only(m):
    outputs = make_outputs_dict_mvc_unit(m)

    # Feed exiting distillate heat exchanger
    outputs['Feed exiting distillate hx temperature'] = m.fs.hx_distillate.cold.properties_out[0].temperature
    outputs['Feed exiting distillate hx pressure'] = m.fs.hx_distillate.cold.properties_out[0].pressure

    # Preheated feed
    outputs['Preheated feed temperature'] =  m.fs.evaporator.properties_feed[0].temperature
    outputs['Preheated feed pressure'] = m.fs.evaporator.properties_feed[0].pressure

    # Distillate heat exchanger performance
    outputs['Distillate hx area'] = m.fs.hx_distillate.area
    outputs['Distillate hx delta temp in'] = m.fs.hx_distillate.delta_temperature_in[0]
    outputs['Distillate hx delta temp out'] = m.fs.hx_distillate.delta_temperature_out[0]
    outputs['Distillate hx heat transfer'] =m.fs.hx_distillate.heat_duty[0]
    outputs['Distillate hx overall heat transfer coefficient'] = m.fs.hx_distillate.overall_heat_transfer_coefficient[0]

    # External Q
    outputs['Q external'] = m.fs.Q_ext[0]

    # Capital costs
    outputs['Evaporator cost per area'] = m.fs.costing.evaporator_unit_cost
    outputs['HX cost per area'] = m.fs.costing.heat_exchanger_unit_cost
    outputs['Compressor unit cost'] = m.fs.costing.compressor_unit_cost
    outputs['Feed pump capital cost'] = m.fs.pump_feed.costing.capital_cost
    outputs['Distillate pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Brine pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Distillate hx capital cost'] = m.fs.hx_distillate.costing.capital_cost
    outputs['Evaporator capital cost'] = m.fs.evaporator.costing.capital_cost
    outputs['Compressor capital cost'] = m.fs.compressor.costing.capital_cost
    outputs['Aggregate capital cost'] = m.fs.costing.aggregate_capital_cost
    outputs['Electricity cost'] = m.fs.costing.electricity_base_cost
    outputs['Aggregate electricity flow cost'] = m.fs.costing.aggregate_flow_costs['electricity']
    outputs['Total investment cost'] = m.fs.costing.total_investment_cost
    outputs['MLC cost'] = m.fs.costing.maintenance_labor_chemical_operating_cost
    outputs['Total operating cost'] = m.fs.costing.total_operating_cost

    # Normalized capital costs
    outputs['CC normalized feed pump'] = m.fs.costing.MVC_captial_cost_percentage['feed_pump']
    outputs['CC normalized distillate pump'] = m.fs.costing.MVC_captial_cost_percentage["distillate_pump"]
    outputs['CC normalized brine pump'] = m.fs.costing.MVC_captial_cost_percentage["brine_pump"]
    outputs['CC normalized distiallte hx'] = m.fs.costing.MVC_captial_cost_percentage["hx_distillate"]
    outputs['CC normalized evaportor'] = m.fs.costing.MVC_captial_cost_percentage["evaporator"]
    outputs['CC normalized compressor'] = m.fs.costing.MVC_captial_cost_percentage["compressor"]

    # Normalized LCOW costs
    outputs['LCOW normalized feed pump'] = m.fs.costing.LCOW_percentage["feed_pump"]
    outputs['LCOW normalized distillate pump'] = m.fs.costing.LCOW_percentage["distillate_pump"]
    outputs['LCOW normalized brine pump'] = m.fs.costing.LCOW_percentage["brine_pump"]
    outputs['LCOW normalized distillate hx'] = m.fs.costing.LCOW_percentage["hx_distillate"]
    outputs['LCOW normalized evaporator'] = m.fs.costing.LCOW_percentage["evaporator"]
    outputs['LCOW normalized compressor'] = m.fs.costing.LCOW_percentage["compressor"]
    outputs['LCOW normalized electricity'] = m.fs.costing.LCOW_percentage['electricity']
    outputs['LCOW normalized MLC'] = m.fs.costing.LCOW_percentage['MLC']
    outputs['LCOW normalized capex'] = m.fs.costing.LCOW_percentage["capital_costs"]
    outputs['LCOW normalized opex'] = m.fs.costing.LCOW_percentage["operating_costs"]
    return outputs

def make_outputs_dict_mvc_full(m):
    outputs = make_outputs_dict_mvc_unit(m)
    # Feed
    outputs['Preheater split ratio'] = m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"]

    # Feed exiting distillate heat exchanger
    outputs['Feed exiting distillate hx temperature'] = m.fs.hx_distillate.cold.properties_out[0].temperature
    outputs['Feed exiting distillate hx pressure'] = m.fs.hx_distillate.cold.properties_out[0].pressure

    # Feed exiting brine heat exchanger
    outputs['Feed exiting brine hx temperature'] = m.fs.hx_brine.cold.properties_out[0].temperature
    outputs['Feed exiting brine hx pressure'] = m.fs.hx_brine.cold.properties_out[0].pressure

    # Preheated feed
    outputs['Preheated feed temperature'] =  m.fs.evaporator.properties_feed[0].temperature
    outputs['Preheated feed pressure'] = m.fs.evaporator.properties_feed[0].pressure

    # Distillate heat exchanger performance
    outputs['Distillate hx area'] = m.fs.hx_distillate.area
    outputs['Distillate hx delta temp in'] = m.fs.hx_distillate.delta_temperature_in[0]
    outputs['Distillate hx delta temp out'] = m.fs.hx_distillate.delta_temperature_out[0]
    outputs['Distillate hx heat transfer'] =m.fs.hx_distillate.heat_duty[0]
    outputs['Distillate hx overall heat transfer coefficient'] = m.fs.hx_distillate.overall_heat_transfer_coefficient[0]

    # Brine heat exchanger performance
    outputs['Brine hx area'] = m.fs.hx_brine.area
    outputs['Brine hx delta temp in'] = m.fs.hx_brine.delta_temperature_in[0]
    outputs['Brine hx delta temp out'] = m.fs.hx_brine.delta_temperature_out[0]
    outputs['Brine hx heat transfer'] =m.fs.hx_brine.heat_duty[0]
    outputs['Brine hx overall heat transfer coefficient'] = m.fs.hx_brine.overall_heat_transfer_coefficient[0]

    # External Q
    outputs['Q external'] = m.fs.Q_ext[0]

    # Capital costs
    outputs['Evaporator cost per area'] = m.fs.costing.evaporator_unit_cost
    outputs['HX cost per area'] = m.fs.costing.heat_exchanger_unit_cost
    outputs['Compressor unit cost'] = m.fs.costing.compressor_unit_cost
    outputs['Feed pump capital cost'] = m.fs.pump_feed.costing.capital_cost
    outputs['Distillate pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Brine pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Distillate hx capital cost'] = m.fs.hx_distillate.costing.capital_cost
    outputs['Brine hx capital cost'] = m.fs.hx_brine.costing.capital_cost
    outputs['Mixer capital cost'] = m.fs.mixer_feed.costing.capital_cost
    outputs['Evaporator capital cost'] = m.fs.evaporator.costing.capital_cost
    outputs['Compressor capital cost'] = m.fs.compressor.costing.capital_cost
    outputs['Aggregate capital cost'] = m.fs.costing.aggregate_capital_cost
    outputs['Electricity cost'] = m.fs.costing.electricity_base_cost
    outputs['Aggregate electricity flow cost'] = m.fs.costing.aggregate_flow_costs['electricity']
    outputs['Total investment cost'] = m.fs.costing.total_investment_cost
    outputs['MLC cost'] = m.fs.costing.maintenance_labor_chemical_operating_cost
    outputs['Total operating cost'] = m.fs.costing.total_operating_cost

    # Normalized capital costs
    outputs['CC normalized feed pump'] = m.fs.costing.MVC_captial_cost_percentage['feed_pump']
    outputs['CC normalized distillate pump'] = m.fs.costing.MVC_captial_cost_percentage["distillate_pump"]
    outputs['CC normalized brine pump'] = m.fs.costing.MVC_captial_cost_percentage["brine_pump"]
    outputs['CC normalized distiallte hx'] = m.fs.costing.MVC_captial_cost_percentage["hx_distillate"]
    outputs['CC normalized brine hx'] = m.fs.costing.MVC_captial_cost_percentage["hx_brine"]
    outputs['CC normalized mixer'] = m.fs.costing.MVC_captial_cost_percentage["mixer"]
    outputs['CC normalized evaportor'] = m.fs.costing.MVC_captial_cost_percentage["evaporator"]
    outputs['CC normalized compressor'] = m.fs.costing.MVC_captial_cost_percentage["compressor"]

    # Normalized LCOW costs
    outputs['LCOW normalized feed pump'] = m.fs.costing.LCOW_percentage["feed_pump"]
    outputs['LCOW normalized distillate pump'] = m.fs.costing.LCOW_percentage["distillate_pump"]
    outputs['LCOW normalized brine pump'] = m.fs.costing.LCOW_percentage["brine_pump"]
    outputs['LCOW normalized distillate hx'] = m.fs.costing.LCOW_percentage["hx_distillate"]
    outputs['LCOW normalized brine hx'] = m.fs.costing.LCOW_percentage["hx_brine"]
    outputs['LCOW normalized mixer'] = m.fs.costing.LCOW_percentage["mixer"]
    outputs['LCOW normalized evaporator'] = m.fs.costing.LCOW_percentage["evaporator"]
    outputs['LCOW normalized compressor'] = m.fs.costing.LCOW_percentage["compressor"]
    outputs['LCOW normalized electricity'] = m.fs.costing.LCOW_percentage['electricity']
    outputs['LCOW normalized MLC'] = m.fs.costing.LCOW_percentage['MLC']
    outputs['LCOW normalized capex'] = m.fs.costing.LCOW_percentage["capital_costs"]
    outputs['LCOW normalized opex'] = m.fs.costing.LCOW_percentage["operating_costs"]
    outputs['capex opex ratio'] = m.fs.costing.LCOW_percentage["capex_opex_ratio"]

    return outputs

def get_param_var_dict(m):
    dict = {}
    dict['evaporator_area'] = m.fs.evaporator.area
    dict['vapor_flow_rate'] = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O']
    dict['pressure_ratio'] = m.fs.compressor.pressure_ratio
    # dict['distillate_hx_area'] = m.fs.hx_distillate.area
    # dict['brine_hx_area'] = m.fs.hx_brine.area
    # dict['split_ratio'] = m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"]
    dict['w_f'] = m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS']
    dict['recovery'] = m.fs.recovery[0]
    dict['electricity_cost'] = m.fs.costing.electricity_base_cost
    dict['evaporator_cost'] = m.fs.costing.evaporator_unit_cost
    dict['preheater_cost'] = m.fs.costing.heat_exchanger_unit_cost
    dict['compressed_vapor_temperature'] = m.fs.compressor.control_volume.properties_out[0].temperature
    dict['U_evap'] = m.fs.evaporator.U
    dict['U_hx_distillate'] = m.fs.hx_distillate.overall_heat_transfer_coefficient
    dict['U_hx_brine'] = m.fs.hx_brine.overall_heat_transfer_coefficient
    dict['T_b'] = m.fs.evaporator.properties_brine[0].temperature
    dict['compressor_cost'] = m.fs.costing.compressor_unit_cost
    dict['compressor_efficiency'] = m.fs.compressor.efficiency
    # dict['T_cv_max'] = m.fs.compressor.control_volume[0].
    return dict

def save_results_for_plotting(results_file, save_dir,n_wf,n_rr):
    df = pd.read_csv(results_file)
    data = np.empty((n_wf, n_rr))
    for param in df:
        print(param)
        row = 0
        for i in range(n_wf):
            for j in range(n_rr):
                data[i][j] = df[param][row]
                row += 1
        # save
        data = np.transpose(data)
        pd.DataFrame(data).to_csv(save_dir+'/'+param+'.csv', index=False)
        data = np.transpose(data)

def convert_units_results(map_dir):
    capex_file = map_dir + '/LCOW normalized capex.csv'
    opex_file = map_dir + '/LCOW normalized opex.csv'
    df_capex = pd.read_csv(capex_file)
    df_opex = pd.read_csv(opex_file)
    df_capex_opex_ratio = df_capex/df_opex
    pd.DataFrame(df_capex_opex_ratio).to_csv(map_dir + '/capex opex ratio.csv', index=False)

    rr_file = map_dir + '/recovery.csv'
    df_rr = pd.read_csv(rr_file)
    feed_mass_flow_water = map_dir + '/Feed mass flow water.csv'
    df_feed_mass_flow_water = pd.read_csv(feed_mass_flow_water)
    feed_mass_flow_salt = map_dir + '/Feed mass flow salt.csv'
    df_feed_mass_flow_salt = pd.read_csv(feed_mass_flow_salt)
    df_feed_mass_flow_total = df_feed_mass_flow_water+df_feed_mass_flow_salt
    df_product_mass_flow = df_rr*df_feed_mass_flow_total
    pd.DataFrame(df_product_mass_flow).to_csv(map_dir + '/Product mass flow.csv', index=False)
    # Calcualte mass flux
    df_evaporator_area = pd.read_csv(map_dir + '/Evaporator area.csv')
    rho = 997.05 # kg/m^3
    conversion_factor = 3600*1000/rho
    df_mass_flow_L_hr = df_product_mass_flow*conversion_factor
    df_mass_flux = df_mass_flow_L_hr.div(df_evaporator_area)
    pd.DataFrame(df_mass_flux).to_csv(map_dir + '/Mass flux LMH.csv', index=False)

    brine_temp_file = map_dir + '/Brine temperature.csv'
    df_brine_temp = pd.read_csv(brine_temp_file)
    feed_temp_file = map_dir + '/Preheated feed temperature.csv'
    df_feed_temp = pd.read_csv(feed_temp_file)
    df_delta = df_brine_temp-df_feed_temp
    pd.DataFrame(df_delta).to_csv(map_dir + '/Evaporator-feed temperature difference.csv', index=False)

    # Convert Q_hx to kW
    Q_distillate_hx_file = map_dir +'/Distillate hx heat transfer.csv'
    Q_brine_hx_file = map_dir +'/Brine hx heat transfer.csv'
    df_Q_distillate_hx = pd.read_csv(Q_distillate_hx_file)
    df_Q_brine_hx = pd.read_csv(Q_brine_hx_file)
    df_Q_distillate_hx = df_Q_distillate_hx/1e3
    df_Q_brine_hx = df_Q_brine_hx/1e3
    pd.DataFrame(df_Q_distillate_hx).to_csv(map_dir+'/Distillate hx heat transfer kW.csv', index=False)
    pd.DataFrame(df_Q_brine_hx).to_csv(map_dir+'/Brine hx heat transfer kW.csv', index=False)
    df_distillate_hx_feed_flow = df_feed_mass_flow_total*df_rr
    df_brine_hx_feed_flow = df_feed_mass_flow_total - df_distillate_hx_feed_flow
    df_q_distillate = df_Q_distillate_hx/df_distillate_hx_feed_flow
    df_q_brine = df_Q_brine_hx/df_brine_hx_feed_flow
    pd.DataFrame(df_q_distillate).to_csv(map_dir+'/Normalized distillate hx heat transfer kJ per kg.csv', index=False)
    pd.DataFrame(df_q_brine).to_csv(map_dir+'/Normalized brine hx heat transfer kJ per kg.csv', index=False)

    # convert to kPa
    results_file = map_dir+'/Brine pressure.csv'
    df = pd.read_csv(results_file)
    df = df*1e-3
    pd.DataFrame(df).to_csv(map_dir+'/Brine pressure kPa.csv', index=False)

    # convert to kPa
    results_file = map_dir + '/Compressed vapor pressure.csv'
    df = pd.read_csv(results_file)
    df = df * 1e-3
    pd.DataFrame(df).to_csv(map_dir + '/Compressed vapor pressure kPa.csv', index=False)

    # Convert to Celsius
    results_file = map_dir+'/Brine temperature.csv'
    df = pd.read_csv(results_file)
    df = df -273.15
    pd.DataFrame(df).to_csv(map_dir+'/Brine temperature Celsius.csv', index=False)

    # Convert to Celsius
    results_file = map_dir+'/Compressed vapor temperature.csv'
    df = pd.read_csv(results_file)
    df = df -273.15
    pd.DataFrame(df).to_csv(map_dir+'/Compressed vapor temperature Celsius.csv', index=False)

    # Convert to Celsius
    results_file = map_dir+'/Distillate temperature.csv'
    df = pd.read_csv(results_file)
    df = df -273.15
    pd.DataFrame(df).to_csv(map_dir+'/Distillate temperature Celsius.csv', index=False)

    # Convert to Celsius
    results_file = map_dir + '/Preheated feed temperature.csv'
    df = pd.read_csv(results_file)
    df = df - 273.15
    pd.DataFrame(df).to_csv(map_dir + '/Preheated feed temperature Celsius.csv', index=False)

    # Get LMTD distillate hx
    dt_in_file = map_dir + '/Distillate hx delta temp in.csv'
    df_dt_in = pd.read_csv(dt_in_file)
    dt_out_file = map_dir + '/Distillate hx delta temp out.csv'
    df_dt_out = pd.read_csv(dt_out_file)
    df_lmtd = (df_dt_in*df_dt_out*(df_dt_in+df_dt_out)*0.5)**(1/3)
    pd.DataFrame(df_lmtd).to_csv(map_dir + '/Distillate hx LMTD.csv', index=False)
    #
    # Get LMTD brine hx
    dt_in_file = map_dir + '/Brine hx delta temp in.csv'
    df_dt_in = pd.read_csv(dt_in_file)
    dt_out_file = map_dir + '/Brine hx delta temp out.csv'
    df_dt_out = pd.read_csv(dt_out_file)
    df_lmtd = (df_dt_in * df_dt_out * (df_dt_in + df_dt_out) * 0.5) ** (1 / 3)
    pd.DataFrame(df_lmtd).to_csv(map_dir + '/Brine hx LMTD.csv', index=False)

if __name__ == "__main__":
    main()