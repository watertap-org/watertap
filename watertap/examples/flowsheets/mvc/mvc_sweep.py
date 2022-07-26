import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    value,
    Expression,
    Param,
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

    # Type 1
    # type_1_parameters = ['evaporator_area', 'pressure_ratio', 'vapor_flow_rate']
    # analysis = "C:/Users/carso/Documents/MVC/watertap_results/analysis_type_1_cases.csv"
    # run_analysis(analysis, type_1_parameters)

    filename="C:/Users/carso/Documents/MVC/watertap_results/type1_multi_sweep_results.csv"
    #global_results, sweep_params, m = run_multi_param_case(5,system='mvc_unit',output_filename=filename)
    # df = pd.read_csv(filename)
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # y = np.array([100, 225 ,350, 475, 600])
    # x = np.array([1.5,1.75,2,2.25,2.5])
    # # evaporator area vs pressure ratio vs compressor work
    # X,Y = np.meshgrid(x,y)
    # R = np.sqrt(X ** 2 + Y ** 2)
    # print(R)
    # Z = np.array([[np.nan,np.nan,np.nan, np.nan,np.nan],
    #      [np.nan,np.nan,np.nan, np.nan,np.nan],
    #      [np.nan,np.nan, 4.04e2, 4.25e2,4.45e2],
    #      [np.nan,3.88e2,4.09e2,4.29e2,4.48e2],
    #      [np.nan,3.91e2,4.11e2,4.30e2,4.49e2]])
    # print(Z)
    # ax.plot_surface(X, Y, Z)
    # ax.set_ylabel(r'Evaporator area [m^2]')
    # ax.set_xlabel(r'Pressure ratio [-]')
    # ax.set_zlabel(r'Compressor work [W]')
    # plt.show()
    #



def mvc_unit_presweep():
    m = mvc_unit.build()
    mvc_unit.set_operating_conditions(m)
    mvc_unit.initialize_system(m)
    mvc_unit.scale_costs(m)
    solver = get_solver()
    results = solver.solve(m, tee=False)
    return m

def mvc_full_presweep():
    m = mvc_full.build()
    mvc_full.set_operating_conditions(m)
    mvc_full.initialize_system(m)
    mvc_full.scale_costs(m)
    solver = get_solver()
    results = solver.solve(m, tee=False)
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

def run_multi_param_case(nx, system='mvc_unit',output_filename=None):
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
    sweep_params["Evaporator area"] = LinearSample(m.fs.evaporator.area, 100, 600, nx)
    sweep_params['Pressure ratio'] = LinearSample(m.fs.compressor.pressure_ratio, 1.5, 2.5, nx)
    #sweep_params['Vapor flow rate'] = LinearSample(m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O'], 3, 7, nx)

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

    if param == 'evaporator_area': # evaporator area
        sweep_params['Evaporator area'] = LinearSample(m.fs.evaporator.area, param_min, param_max, n_param)
    elif param == 'pressure_ratio': # pressure ratio
        sweep_params['Pressure ratio'] = LinearSample(m.fs.compressor.pressure_ratio, param_min, param_max, n_param)
    elif param == 'vapor_flow_rate':
        sweep_params['Vapor flow rate'] = LinearSample(m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O'], param_min, param_max, n_param)
    else:
        print('Provide parameter to be swept.')
        assert False

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

def run_mvc_full_case(n_param, output_filename=None):
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

    # Get model
    m = mvc_full_presweep()

    # Sweep parameters
    sweep_params = {}
    sweep_params['Vapor flow rate'] = LinearSample(m.fs.evaporator.properties_vapor[0].flow_mass_phase_comp['Vap','H2O'],4.5,5.5, n_param)

    # Outputs
    outputs = make_outputs_dict(m)

    # sweep
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=mvc_full.solve,
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
    outputs['Brine temperature'] = m.fs.evaporator.properties_brine[0].temperature
    outputs['Brine pressure'] = m.fs.evaporator.properties_brine[0].pressure

    # Vapor
    outputs['Vapor temperature'] = m.fs.evaporator.properties_vapor[0].temperature
    outputs['Vapor pressure'] = m.fs.evaporator.properties_vapor[0].pressure

    # Compressed vapor
    outputs['Compressed vapor temperature'] = m.fs.compressor.control_volume.properties_out[0].temperature
    outputs['Compressed vapor pressure'] = m.fs.compressor.control_volume.properties_out[0].pressure

    # Condensed vapor/distillate
    outputs['Distillate temperature'] = m.fs.condenser.control_volume.properties_out[0].temperature
    outputs['Distillate pressure'] = m.fs.condenser.control_volume.properties_out[0].pressure

    # Evaporator performance
    outputs['Evaporator area'] = m.fs.evaporator.area
    outputs['Evaporator LMTD'] = m.fs.evaporator.lmtd
    outputs['Evaporator heat transfer'] = m.fs.evaporator.heat_transfer
    outputs['Evaporator overall heat transfer coefficient'] = m.fs.evaporator.U

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

    return outputs



if __name__ == "__main__":
    main()