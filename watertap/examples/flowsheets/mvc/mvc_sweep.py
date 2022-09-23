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
import mvc_plotting as mvc_plot


def main():
    analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_full_optimize_cases.csv"
    map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/"

    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/C_elec_sensitivity/C_elec_25"
    # output_file = map_dir + "/optimize_full_sweep_rr_wf_40.csv"
    # global_results, sweep_params, m = run_multi_param_case(analysis,system='mvc_full_opt',output_filename=output_file,f_evap=6000, f_hx=3000, C_elec=0.25)
    # save_dir = map_dir + '/figures'
    # save_results_for_plotting(output_file,map_dir,7,9)
    # convert_units_results(map_dir)
    # mvc_plot.make_maps(map_dir,save_dir)
    # print(map_dir)
    #


    # dir = "C:/Users/carso/Documents/MVC/watertap_results/dual_c_evap_c_elec_sensitivity"
    # analysis = dir + "/analysis_dual_c_evap_c_elec.csv"
    # map_dir = dir + "/lower_costs"
    # output_file = map_dir + "/optimize_sweep.csv"
    # global_results, sweep_params, m = run_multi_param_case(analysis, system='mvc_full_opt', output_filename=output_file)

    dir = "C:/Users/carso/Documents/MVC/watertap_results/T_cv_sensitivity"
    analysis = dir + "/analysis_T_cv_sensitivity.csv"
    map_dir = dir
    output_file = map_dir + "/evap_6000_hx_3000/optimize_sweep.csv"
    global_results, sweep_params, m = run_multi_param_case(analysis, system='mvc_full_opt', output_filename=output_file, f_evap=6000, f_hx=3000)

def mvc_unit_presweep():
    m = mvc_unit.build()
    mvc_unit.set_operating_conditions(m)
    mvc_unit.initialize_system(m)
    mvc_unit.scale_costs(m)
    solver = get_solver()
    results = solver.solve(m, tee=False)
    return m

def mvc_full_presweep(f_evap=1000, f_hx=1000, T_cv_max=500, T_b=None, C_elec=0.07, C_evap_hx_ratio=None):
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
    mvc_full.initialize_system(m)
    mvc_full.scale_costs(m)
    # mvc_full.fix_outlet_pressures(m)
    solver = get_solver()
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    results = solver.solve(m, tee=False)
    results = mvc_full.sweep_solve(m)
    mvc_full.display_results(m)
    m.fs.compressor.control_volume.properties_out[0].temperature.setub(T_cv_max)
    if T_b is not None:
        m.fs.evaporator.properties_brine[0].temperature.fix(T_b)
    return m

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

def run_multi_param_case(analysis_file, system='mvc_unit',output_filename=None,f_evap=1000,f_hx=1000, T_cv_max=500, T_b=None, C_elec=0.07, C_evap_hx_ratio=None):
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
        m = mvc_full_presweep(f_evap=f_evap,f_hx=f_hx, T_cv_max=T_cv_max, T_b=T_b, C_elec=C_elec,C_evap_hx_ratio=C_evap_hx_ratio)
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

def run_case(n_param, param=None, param_min=None, param_max=None, system='mvc_unit', output_filename=None):
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
    outputs['Feed pump capital cost'] = m.fs.pump_feed.costing.capital_cost
    outputs['Distillate pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Brine pump captial cost'] = m.fs.pump_distillate.costing.capital_cost
    outputs['Distillate hx capital cost'] = m.fs.hx_distillate.costing.capital_cost
    outputs['Brine hx capital cost'] = m.fs.hx_brine.costing.capital_cost
    outputs['Mixer capital cost'] = m.fs.mixer_feed.costing.capital_cost
    outputs['Evaporator capital cost'] = m.fs.evaporator.costing.capital_cost
    outputs['Compressor capital cost'] = m.fs.compressor.costing.capital_cost
    outputs['Aggregate capital cost'] = m.fs.costing.aggregate_capital_cost
    outputs['Aggregate electricity flow cost'] = m.fs.costing.aggregate_flow_costs['electricity']
    outputs['Total investment cost'] = m.fs.costing.total_investment_cost
    outputs['MLC cost'] = m.fs.costing.maintenance_labor_chemical_operating_cost
    outputs['Total operating cost'] = m.fs.costing.total_operating_cost

    return outputs

def get_param_var_dict(m):
    dict = {}
    dict['evaporator_area'] = m.fs.evaporator.area
    dict['vapor_flow_rate'] = m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O']
    dict['pressure_ratio'] = m.fs.compressor.pressure_ratio
    dict['distillate_hx_area'] = m.fs.hx_distillate.area
    dict['brine_hx_area'] = m.fs.hx_brine.area
    dict['split_ratio'] = m.fs.separator_feed.split_fraction[0, "hx_distillate_cold"]
    dict['w_f'] = m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'TDS']
    dict['recovery'] = m.fs.recovery[0]
    dict['electricity_cost'] = m.fs.costing.electricity_base_cost
    dict['evaporator_cost'] = m.fs.costing.evaporator_unit_cost
    dict['preheater_cost'] = m.fs.costing.heat_exchanger_unit_cost
    dict['compressed_vapor_temperature'] = m.fs.compressor.control_volume.properties_out[0].temperature
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

def get_feasibility_3d_plots():
    filename = "C:/Users/carso/Documents/MVC/watertap_results/mvc_unit_vapor_flow_rate_results.csv"
    # naming of files: type_parameter_being_varied_results.csv
    # Type 1: feed conditions, evaporator area, compressor pressure ratio, vapor flow rate fixed
    # Type 2: feed conditions, evaporator temperature, compressor pressure ratio, vapor flow rate fixed
    # Type 3: feed conditions, evaporator area, evaporator heat duty, compressor work

    # Type 1 - evap-comp-cond
    type_1_parameters = ['evaporator_area', 'pressure_ratio', 'vapor_flow_rate']
    analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_type_1_cases.csv"
    output_file = "C:/Users/carso/Documents/MVC/watertap_results/type1_multi_sweep_results.csv"
    run_multi_param_case(analysis, system='mvc_unit', output_filename=output_file)
    plot_3D_results(output_file)

    # Type 1 - full single stage
    analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_type_1_full_cases.csv"
    output_file = "C:/Users/carso/Documents/MVC/watertap_results/type1_full_multi_sweep_results.csv"
    run_multi_param_case(analysis, system='mvc_full', output_filename=output_file)
    plot_3D_results(output_file)

if __name__ == "__main__":
    main()