from parameter_sweep import LinearSample, parameter_sweep
from pyomo.environ import units as pyunits
from watertap.costing import WaterTAPCosting  
from watertap.flowsheets.crystallization.Crystallizer_MVR import (
    build,
    set_operating_conditions,
    initialize_system,
    optimize_set_up,
    solve,
    get_solver,
) 
import pandas as pd

def set_up_sensitivity(m):
    outputs = {}
    optimize_kwargs = {}
    opt_function = solve
    outputs["LCOW"] = m.fs.costing.LCOW
   
    return outputs, optimize_kwargs, opt_function

def run_electricity_price_sweep(nx=20, output_filename="electricity_price_sweep.csv"):
    solver = get_solver()
    m = build()
    set_operating_conditions(m)
    initialize_system(m, solver=solver)
    optimize_set_up(m)
    
    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m)
    
    sweep_params = {}
    sweep_params["electricity_cost"] = LinearSample(
        m.fs.costing.electricity_cost, 0.00, 0.5, nx
    )
    
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=opt_function,
        optimize_kwargs=optimize_kwargs,
    )
    return global_results, sweep_params, m

if __name__ == "__main__":
    results, sweep_params, m = run_electricity_price_sweep()
    df = pd.read_csv("electricity_price_sweep.csv")
    lcow_values = df["LCOW"].tolist()
    print("Extracted LCOW results:", lcow_values)