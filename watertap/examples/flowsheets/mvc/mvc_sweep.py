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


def main():
    #filename = "C:/Users/carso/Documents/MVC/watertap_results/mvc_unit_vapor_flow_rate_results.csv"
    # naming of files: type_parameter_being_varied_results.csv
    # Type 1: feed conditions, evaporator area, compressor pressure ratio, vapor flow rate fixed
    # Type 2: feed conditions, evaporator temperature, compressor pressure ratio, vapor flow rate fixed
    # Type 3: feed conditions, evaporator area, evaporator heat duty, compressor work

    # Type 1 - evap-comp-cond
    # type_1_parameters = ['evaporator_area', 'pressure_ratio', 'vapor_flow_rate']
    #analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_type_1_cases.csv"
    #output_file = "C:/Users/carso/Documents/MVC/watertap_results/type1_multi_sweep_results.csv"
    #run_multi_param_case(analysis, system='mvc_unit', output_filename=output_file)
    #plot_3D_results(output_file)

    # Type 1 - full single stage
    # analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_type_1_full_cases.csv"
    # output_file = "C:/Users/carso/Documents/MVC/watertap_results/type1_full_multi_sweep_results.csv"
    # run_multi_param_case(analysis, system='mvc_full', output_filename=output_file)
    # plot_3D_results(output_file)

    # Optimize full system
    analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_full_optimize_cases.csv"
    map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/"
    # map_dir_list = [map_dir+'evap_3000_hx_3000',
    #                 map_dir+'evap_4000_hx_4000',
    #                 map_dir+'evap_5000_hx_5000',
    #                map_dir+'evap_6000_hx_6000']
    # title_list = [r'$C_{evap}=C_{hx}=3000 \$/m^2$',
    #               r'$C_{evap}=C_{hx}=4000 \$/m^2$',
    #               r'$C_{evap}=C_{hx}=5000 \$/m^2$',
    #               r'$C_{evap}=C_{hx}=6000 \$/m^2$']
    # save_dir = map_dir+"figures_same_cost_3000_6000"
    # make_maps_comparison(map_dir_list, save_dir, title_list)

    map_dir_list = [map_dir + 'evap_6000_hx_4800',
                    map_dir + 'evap_6000_hx_4000',
                    map_dir + 'evap_6000_hx_3429',
                    map_dir + 'evap_6000_hx_3000']
    title_list = [r'$C_{evap}=6000 \$/m^2, C_{evap}=1.25\times C_{hx}$',
                  r'$C_{evap}=6000 \$/m^2, C_{evap}=1.5\times C_{hx}$',
                  r'$C_{evap}=5000 \$/m^2, C_{evap}=1.75\times C_{hx}$',
                  r'$C_{evap}=6000 \$/m^2, C_{evap}=2\times C_{hx}$']
    save_dir = map_dir+"figures_evap_6000_varying_cost_ratio"
    make_maps_comparison(map_dir_list, save_dir, title_list)

    # map_dir_list = [map_dir + 'evap_1000_hx_1500',
    #                 map_dir + 'evap_1500_hx_2000']
    # title_list = [r'$C_{evap}=1000 \$/m^2, C_{hx}=1500 \$/m^2$',
    #               r'$C_{evap}=1500 \$/m^2, C_{hx}=2000 \$/m^2$',]
    # save_dir = map_dir+"figures_evap_lower_cost"
    # make_maps_comparison(map_dir_list, save_dir, title_list)

    # # evap 800 hx 800
    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_linearized/evap_6000_hx_4800"
    # output_file = map_dir + "/optimize_full_sweep_rr_wf_40.csv"
    # global_results, sweep_params, m = run_multi_param_case(analysis,system='mvc_full_opt',output_filename=output_file,f_evap=6000, f_hx=4800)
    # save_dir = map_dir + '/figures'
    # save_results_for_plotting(output_file,map_dir,7,9)
    # convert_units_results(map_dir)
    # make_maps_F_200(map_dir,save_dir)
    # print(map_dir)

    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/Total investment factor sensitivity/opt_full_40_tif_4"
    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min"
    # save_dir = map_dir + '/figures'
    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/cost_breakdown/Linearized evap 3 hx 3"
    # make_cost_bar_charts(map_dir)

def sr_sensitivity_maps():
    # lcow_file = "C:/Users/carso/Documents/MVC/watertap_results/Total investment factor sensitivity/opt_full_40_tif_4/LCOW.csv"
    # df_lcow = pd.read_csv(lcow_file)
    # lcow_fixed_sr_file = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/LCOW.csv"
    # df_lcow_fixed = pd.read_csv(lcow_fixed_sr_file)
    # df_change = df_lcow_fixed.div(df_lcow) - 1
    # df_change = df_change*100
    # print(df_change)
    # save_file = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/LCOW change.csv"
    # pd.DataFrame(df_change).to_csv(save_file, index=False)

    # # Make change plot
    # map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min"
    # save_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/figures"
    # var = 'LCOW change'
    # label = r'Percentage change in LCOW [%]'
    # vmin = 0  # minimum cost on bar, $/m3
    # vmax = 3  # maximum cost on bar, $/m3
    # ticks = [0, 1, 2, 3]  # tick marks on bar
    # fmt = '.2f'  # format of annotation
    # plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    temp_file = "C:/Users/carso/Documents/MVC/watertap_results/Total investment factor sensitivity/opt_full_40_tif_4/Brine temperature Celsius.csv"
    df_temp = pd.read_csv(temp_file)
    temp_fixed_sr_file = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/Brine temperature Celsius.csv"
    df_temp_fixed = pd.read_csv(temp_fixed_sr_file)
    df_change = df_temp_fixed.div(df_temp) - 1
    df_change = df_change*100
    print(df_change)
    save_file = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/Brine temperature change.csv"
    pd.DataFrame(df_change).to_csv(save_file, index=False)

    # Make change plot
    map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min"
    save_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_4_fixed_split_ratio_Q_ext_min/figures"
    var = 'Brine temperature change'
    label = r'Percentage change in evaporator temperature [%]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 3  # maximum cost on bar, $/m3
    ticks = [0, 1, 2, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

def tif_sensitivity_maps():
    # dir = "C:/Users/carso/Documents/MVC/watertap_results/Total investment factor sensitivity"
    #
    # map_dir = dir + '/opt_full_40_tif_0p5'
    # save_dir = map_dir + '/figures_comparison'
    # #make_maps_comparison(map_dir,save_dir)
    #
    # map_dir = dir + '/opt_full_40_tif_1'
    # save_dir = map_dir + '/figures_comparison'
    # make_maps_comparison(map_dir,save_dir)
    #
    # map_dir = dir + '/opt_full_40_tif_2'
    # save_dir = map_dir + '/figures_comparison'
    # make_maps_comparison(map_dir,save_dir)
    #
    # map_dir = dir + '/opt_full_40_tif_4'
    # save_dir = map_dir + '/figures_comparison'
    # make_maps_comparison(map_dir,save_dir)
    #
    # map_dir = dir + '/opt_full_40_tif_8'
    # save_dir = map_dir + '/figures_comparison'
    # make_maps_comparison(map_dir,save_dir)

    map_dir = "C:/Users/carso/Documents/MVC/watertap_results/opt_full_40_tif_8_fixed_split_ratio"
    save_dir = map_dir + '/figures_comparison'
    make_maps_comparison(map_dir,save_dir)

def mvc_unit_presweep():
    m = mvc_unit.build()
    mvc_unit.set_operating_conditions(m)
    mvc_unit.initialize_system(m)
    mvc_unit.scale_costs(m)
    solver = get_solver()
    results = solver.solve(m, tee=False)
    return m

def mvc_full_presweep(f_evap=1000, f_hx=1000):
    m = mvc_full.build()
    mvc_full.add_Q_ext(m, time_point=m.fs.config.time)
    mvc_full.set_operating_conditions(m)
    m.fs.costing.heat_exchanger_unit_cost.fix(f_evap)
    m.fs.costing.evaporator_unit_cost.fix(f_hx)
    mvc_full.initialize_system(m)
    mvc_full.scale_costs(m)
    # mvc_full.fix_outlet_pressures(m)
    solver = get_solver()
    m.fs.objective = Objective(expr=m.fs.Q_ext[0])
    results = solver.solve(m, tee=False)
    results = mvc_full.sweep_solve(m)
    mvc_full.display_results(m)
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

def run_multi_param_case(analysis_file, system='mvc_unit',output_filename=None,f_evap=1000,f_hx=1000):
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
        m = mvc_full_presweep(f_evap=f_evap,f_hx=f_hx)
        opt_fcn = mvc_full.solve
        outputs = make_outputs_dict_mvc_full(m)

    elif system == 'mvc_full_opt':
        m = mvc_full_presweep(f_evap=f_evap,f_hx=f_hx)
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
    return dict

def plot_3D_results(results_file):
    df = pd.read_csv(results_file)
    outputs = list(df.columns)
    param_x_name = outputs[0]
    param_y_name = outputs[1]
    param_z_name = outputs[2]
    param_check = outputs[9]

    # get feasible results
    df_feasible = df[~df[param_check].isnull()]

    # make scatter plot
    plt.figure()
    axes=plt.axes(projection='3d')
    x = np.array(df_feasible[param_x_name])
    y = np.array(df_feasible[param_y_name])
    z = np.array(df_feasible[param_z_name])
    #color = np.array(df_feasible['Brine temperature']-273.15)
    color = np.array(df_feasible['LCOW'])
    color = np.array(df_feasible['Q external'])

    tick_min = min(color)
    tick_max = max(color)
    fig = axes.scatter3D(x,y,z,c=color) #, cbar_kws={'label': 'Brine temperature'})
    axes.set_xlabel('Evaporator area (m2)')
    axes.set_ylabel('Pressure ratio (-)')
    axes.set_zlabel('Vapor flow rate (kg/s)')
    #plt.colorbar(fig,label='LCOW ($/m3)')
    plt.colorbar(fig, label='Q external')
    plt.show()

def plot_2D_heat_map(map_dir, save_dir, param, label, vmin, vmax, ticks, fmt,make_ticks=True):
    fig = plt.figure()
    ax = plt.axes()
    results_file = map_dir+'/'+param+'.csv'
    df = pd.read_csv(results_file)
    mask = df.isnull()
    if make_ticks:
        decimal=int(fmt[1])
        df_min = float(np.nanmin(df.values))
        vmin = df_min-10**-decimal
        df_max = float(np.nanmax(df.values))
        vmax = df_max + 10**-decimal
        n=5 # number of ticks
        ticks = np.round(np.linspace(df_min, df_max, n),decimal) # round based on formatting decimal places

    xticklabels = ['25', '50', '75', '100', '125', '150', '175']
    yticklabels = ['40', '45', '50', '55', '60', '65', '70', '75', '80']
    # ax = sns.heatmap(df, cmap='Reds', mask=mask,
    #                  vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
    #                  cbar_kws={'label': label, "ticks": ticks}, xticklabels=xticklabels,
    #                  yticklabels=yticklabels
    #                  )  # create heatmap
    ax = sns.heatmap(df, cmap='YlGnBu', mask=mask,
                     vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
                     cbar_kws={'label': label, "ticks": ticks}, xticklabels=xticklabels,
                     yticklabels=yticklabels
                     )  # create heatmap
    ax.invert_yaxis()
    plt.yticks(rotation=0)
    ax.set_xlabel('Feed concentration (g/kg)')
    ax.set_ylabel('Water recovery (%)')
    fig.set_size_inches(3.25, 3.25)
    # plt.show()
    fig.savefig(save_dir + '/' + param + '.png',bbox_inches='tight',dpi=300)

    # fig.savefig(save_dir + '/' + var + '.pdf',bbox_inches='tight')
    # fig.savefig(save_dir + '/' + var + '.svg',bbox_inches='tight')
    #fig.savefig(save_dir + '/' + var + '.png', bbox_inches='tight', dpi=300)

def plot_2D_heat_map_subplots(map_dir, title_list, save_dir, param, param_label, vmin, vmax, ticks, fmt, make_ticks=True,show=False):
    # map dir is a list
    n = len(map_dir)

    if make_ticks:
        decimal=int(fmt[1])
        # first file
        results_file = map_dir[0] + '/' + param + '.csv'
        df = pd.read_csv(results_file)
        df_min = float(np.nanmin(df.values))
        df_max = float(np.nanmax(df.values))

        # search for max and mins
        for i in range(1,n):
            results_file = map_dir[i] + '/' + param + '.csv'
            df = pd.read_csv(results_file)
            min_val = float(np.nanmin(df.values))
            if min_val < df_min:
                df_min = min_val
            max_val = float(np.nanmax(df.values))
            if max_val > df_max:
                df_max = max_val
        # make ticks
        vmin = df_min - 10 ** -decimal
        vmax = df_max + 10 ** -decimal
        n_ticks=5 # number of ticks
        ticks = np.round(np.linspace(df_min, df_max, n_ticks),decimal) # round based on formatting decimal places

    widths = [1 for i in range(n)]
    heights = [1 for i in range(n)]
    widths.append(0.08)
    fig, ax = plt.subplots(1,n+1, gridspec_kw={'width_ratios':widths},figsize=(n*3,3))# figsize=(12,3))
    # ax[0].get_shared_y_axes().join(ax[1],ax[2],ax[3])
    # ax[1].get_shared_y_axes().join(ax[2],ax[3])

    cbar = [False for i in range(n-1)]
    cbar.append(True)
    xticklabels = ['25', '50', '75', '100', '125', '150', '175']
    yticklabels = ['40', '45', '50', '55', '60', '65', '70', '75', '80']
    yticks = [[]for i in range(n)]
    yticks[0] = yticklabels
    g = {}
    for i in range(n):
        results_file = map_dir[i] + '/' + param + '.csv'
        df = pd.read_csv(results_file)
        mask = df.isnull()
        ax[i].axis('equal')
        g[i] = sns.heatmap(df, cmap='YlGnBu', mask=mask, square=True,
                         vmin=vmin, vmax=vmax, annot=True, annot_kws={"fontsize": 8}, fmt=fmt,
                         cbar_kws={'label': param_label, "ticks": ticks}, xticklabels=xticklabels,
                         yticklabels=yticklabels,ax=ax[i],cbar=cbar[i],cbar_ax=ax[i+1]
                         )  # create heatmap

        #g[i].set_xlabel('Feed concentration (g/kg)')
        ax[i].set_yticklabels(yticklabels,rotation=0)
        ax[i].set_xticklabels(xticklabels,rotation=0)
        ax[i].set_title(title_list[i],size=8)
        g[i].invert_yaxis()
        g[i].set_xlabel('Feed concentration (g/kg)')

    g[0].set_ylabel('Water recovery (%)')
    if show:
        plt.show()
    fig.savefig(save_dir + '/' + param + '.png',bbox_inches='tight',dpi=300)

    # fig.savefig(save_dir + '/' + var + '.pdf',bbox_inches='tight')
    # fig.savefig(save_dir + '/' + var + '.svg',bbox_inches='tight')
    #fig.savefig(save_dir + '/' + var + '.png', bbox_inches='tight', dpi=300)

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

def make_maps_F_200(map_dir, save_dir):
    var = 'LCOW'
    label = r'LCOW [\$/$\rmm^3$ of product]'
    vmin = 2.8  # minimum cost on bar, $/m3
    vmax = 6.1 # maximum cost on bar, $/m3
    ticks = [3,4,5,6]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map(map_dir,save_dir,var,label,vmin,vmax,ticks,fmt)

    var = 'SEC'
    label = r'SEC [kWh/$\rmm^3$ of product]'
    vmin = 15  # minimum cost on bar, $/m3
    vmax = 63  # maximum cost on bar, $/m3
    ticks = [20,30,40,50,63]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir,save_dir,var,label,vmin,vmax,ticks,fmt)

    var = 'Brine temperature Celsius'
    label = 'Evaporator temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [25,50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine pressure kPa'
    label = 'Evaporator pressure [kPa]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 403  # maximum cost on bar, $/m3
    ticks = [0, 100, 200, 300, 400]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Preheated feed temperature Celsius'
    label = 'Preheated feed temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor temperature Celsius'
    label = 'Compressed vapor temperature [C]'
    vmin = 87  # minimum cost on bar, $/m3
    vmax = 227  # maximum cost on bar, $/m3
    ticks = [100, 150, 200]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressed vapor pressure kPa'
    label = 'Compressed vapor pressure [kPa]'
    vmin = 3  # minimum cost on bar, $/m3
    vmax = 629  # maximum cost on bar, $/m3
    ticks = [200, 400, 600]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate temperature Celsius'
    label = 'Distillate temperature [C]'
    vmin = 34  # minimum cost on bar, $/m3
    vmax = 161  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Preheater split ratio'
    label = 'Preheater split ratio [C]'
    vmin = 0.4  # minimum cost on bar, $/m3
    vmax = 1  # maximum cost on bar, $/m3
    ticks = [0.4, 0.6, 0.8]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Distillate hx area'
    label = r'Distillate preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1630  # maximum cost on bar, $/m3
    ticks = [500,1000, 2000, 3000, 4000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Brine hx area'
    label = r'Brine preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1001  # maximum cost on bar, $/m3
    ticks = [250,500, 750]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator area'
    label = r'Evaporator area [$\rmm^2$]'
    vmin = 742  # minimum cost on bar, $/m3
    vmax = 3740  # maximum cost on bar, $/m3
    ticks = [1000, 2000, 3000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator LMTD'
    label = r'Evaporator LMTD [K]'
    vmin = 19  # minimum cost on bar, $/m3
    vmax = 61  # maximum cost on bar, $/m3
    ticks = [20, 40, 60]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Compressor pressure ratio'
    label = r'Compressor pressure ratio [-]'
    vmin = 1.5  # minimum cost on bar, $/m3
    vmax = 3.6  # maximum cost on bar, $/m3
    ticks = [2,2.5, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Evaporator-feed temperature difference'
    label = r'$T_{brine}-T_{feed}$ [C]'
    vmin = -7  # minimum cost on bar, $/m3
    vmax = 15  # maximum cost on bar, $/m3
    ticks = [-5,0, 5, 10, 15]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

    var = 'Mass flux LMH'
    label = r'Product flux over evaporator [LMH]'
    vmin = 28  # minimum cost on bar, $/m3
    vmax = 82  # maximum cost on bar, $/m3
    ticks = [30, 40,50,60,70,80]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map(map_dir, save_dir, var, label, vmin, vmax, ticks, fmt)

def make_maps_comparison(map_dir_list, save_dir, title_list):
    var = 'LCOW'
    label = r'LCOW [\$/$\rmm^3$ of product]'
    vmin = 2.8  # minimum cost on bar, $/m3
    vmax = 6.1 # maximum cost on bar, $/m3
    ticks = [3,4,5,6]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'SEC'
    label = r'SEC [kWh/$\rmm^3$ of product]'
    vmin = 15  # minimum cost on bar, $/m3
    vmax = 63  # maximum cost on bar, $/m3
    ticks = [20,30,40,50,63]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Brine temperature Celsius'
    label = 'Evaporator temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [25,50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Brine pressure kPa'
    label = 'Evaporator pressure [kPa]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 403  # maximum cost on bar, $/m3
    ticks = [0, 100, 200, 300, 400]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Preheated feed temperature Celsius'
    label = 'Preheated feed temperature [C]'
    vmin = 25  # minimum cost on bar, $/m3
    vmax = 150  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Compressed vapor temperature Celsius'
    label = 'Compressed vapor temperature [C]'
    vmin = 87  # minimum cost on bar, $/m3
    vmax = 227  # maximum cost on bar, $/m3
    ticks = [100, 150, 200]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Compressed vapor pressure kPa'
    label = 'Compressed vapor pressure [kPa]'
    vmin = 3  # minimum cost on bar, $/m3
    vmax = 629  # maximum cost on bar, $/m3
    ticks = [200, 400, 600]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Distillate temperature Celsius'
    label = 'Distillate temperature [C]'
    vmin = 34  # minimum cost on bar, $/m3
    vmax = 161  # maximum cost on bar, $/m3
    ticks = [50, 100, 150]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Preheater split ratio'
    label = 'Preheater split ratio [C]'
    vmin = 0.4  # minimum cost on bar, $/m3
    vmax = 1  # maximum cost on bar, $/m3
    ticks = [0.4, 0.6, 0.8]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Distillate hx area'
    label = r'Distillate preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1630  # maximum cost on bar, $/m3
    ticks = [500,1000, 2000, 3000, 4000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Brine hx area'
    label = r'Brine preheater area [$\rmm^2$]'
    vmin = 0  # minimum cost on bar, $/m3
    vmax = 1001  # maximum cost on bar, $/m3
    ticks = [250,500, 750]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Evaporator area'
    label = r'Evaporator area [$\rmm^2$]'
    vmin = 742  # minimum cost on bar, $/m3
    vmax = 3740  # maximum cost on bar, $/m3
    ticks = [1000, 2000, 3000]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Evaporator LMTD'
    label = r'Evaporator LMTD [K]'
    vmin = 19  # minimum cost on bar, $/m3
    vmax = 61  # maximum cost on bar, $/m3
    ticks = [20, 40, 60]  # tick marks on bar
    fmt = '.1f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Compressor pressure ratio'
    label = r'Compressor pressure ratio [-]'
    vmin = 1.5  # minimum cost on bar, $/m3
    vmax = 3.6  # maximum cost on bar, $/m3
    ticks = [2,2.5, 3]  # tick marks on bar
    fmt = '.2f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Evaporator-feed temperature difference'
    label = r'$T_{brine}-T_{feed}$ [C]'
    vmin = -7  # minimum cost on bar, $/m3
    vmax = 15  # maximum cost on bar, $/m3
    ticks = [-5,0, 5, 10, 15]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)

    var = 'Mass flux LMH'
    label = r'Product flux over evaporator [LMH]'
    vmin = 28  # minimum cost on bar, $/m3
    vmax = 82  # maximum cost on bar, $/m3
    ticks = [30, 40,50,60,70,80]  # tick marks on bar
    fmt = '.0f'  # format of annotation
    plot_2D_heat_map_subplots(map_dir_list, title_list, save_dir, var, label, vmin, vmax, ticks, fmt, make_ticks=True)


def make_cost_bar_charts(map_dir):
    results_file = map_dir + '/costs_cases.csv'
    df = pd.read_csv(results_file,index_col=[0]) # to make the row names the keys
    n_cases = df.shape[1]
    x = []
    pump_cc = []
    distillate_hx_cc = []
    brine_hx_cc =[]
    mixer_cc = []
    evap_cc = []
    comp_cc = []
    mlc_oc = []
    elec_oc = []
    lcow = []
    for i in range(n_cases):
        case_i = 'Case ' + str(i+1)
        #label = case_i + '\n' + str(round(df[case_i]['Feed concentration'],0)) + ' g/kg\n' + str(round(df[case_i]['Recovery'],0)) + '%'
        label = str(round(df[case_i]['Feed concentration'])) + ' g/kg\n' + str(round(df[case_i]['Recovery'])) + '%'
        x.append(label)
        #normalized LCOW costs
        pump_cc.append(df[case_i]['LCOW normalized feed pump'] + df[case_i]['LCOW normalized distillate pump'] + df[case_i]['LCOW normalized brine pump'])
        distillate_hx_cc.append(df[case_i]['LCOW normalized distillate hx'])
        brine_hx_cc.append(df[case_i]['LCOW normalized brine hx'])
        mixer_cc.append(df[case_i]['LCOW normalized mixer'])
        evap_cc.append(df[case_i]['LCOW normalized evaporator'])
        comp_cc.append(df[case_i]['LCOW normalized compressor'])
        mlc_oc.append(df[case_i]['LCOW normalized MLC'])
        elec_oc.append(df[case_i]['LCOW normalized electricity'])
        lcow.append(df[case_i]['LCOW'])

    pump_cc = tuple(pump_cc)
    distillate_hx_cc = tuple(distillate_hx_cc)
    brine_hx_cc = tuple(brine_hx_cc)
    mixer_cc = tuple(mixer_cc)
    evap_cc = tuple(evap_cc)
    comp_cc = tuple(comp_cc)
    mlc_oc = np.array(mlc_oc)
    elec_oc = np.array(elec_oc)

    # Make plots
    width = 0.8
    ind = np.arange(n_cases)
    fig,ax = plt.subplots(figsize=(5,4))
    p1 = plt.bar(x,elec_oc,width)
    p2 = plt.bar(x,mlc_oc, width, bottom = elec_oc)
    p3 = plt.bar(x,pump_cc, width, bottom = elec_oc + mlc_oc)
    p4 = plt.bar(x,distillate_hx_cc, width, bottom =elec_oc + mlc_oc)
    p5 = plt.bar(x,brine_hx_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc)
    p6 = plt.bar(x,mixer_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc)
    p7 = plt.bar(x, comp_cc, width, bottom= elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc)
    p8 = plt.bar(x, evap_cc, width, bottom=elec_oc + mlc_oc + distillate_hx_cc + brine_hx_cc + mixer_cc + comp_cc)

    for i,data in enumerate(lcow):
        plt.text(x=i,y=1 + 0.01, s=f"{round(data,2)}" + r" $\$/m^3$", ha="center", fontsize=8)

    plt.ylabel ('Normalized LCOW (-)')
    plt.legend((p8[0],p7[0],p6[0],p5[0], p4[0], p3[0], p2[0], p1[0]),
               ('Evaporator', 'Compressor','Mixer', 'Brine HX','Distillate HX','Pumps', 'MLC', 'Electricity'),
               loc='center left',
               bbox_to_anchor=(1, 0.5),
               prop={'size': 8}
               )
    plt.tight_layout()
    ax.xaxis.label.set_size(6)
    plt.show()

if __name__ == "__main__":
    main()