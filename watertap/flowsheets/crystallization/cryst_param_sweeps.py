from parameter_sweep import LinearSample, parameter_sweep

from watertap.flowsheets.crystallization import (
    crystallizer_live_steam_with_condenser_chiller as steam,
)
import Crystallizer_MVR as mvc
import crystallizer_TVC as tvc

from watertap.core.solvers import get_solver
def set_up_sensitivity(m, cryst_module):
    outputs = {}
    optimize_kwargs = {"tee": False}
    opt_function = cryst_module.solve
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    return outputs, optimize_kwargs, opt_function

def run_cryst_sweep(cryst_module, nx=50, output_filename="live_cryst_sweep.csv"):
    solver = get_solver()
    m = cryst_module.build()
    cryst_module.set_operating_conditions(m)
    cryst_module.initialize_system(m, solver=solver)
    cryst_module.optimize_set_up(m)

    outputs, optimize_kwargs, opt_function = set_up_sensitivity(m, cryst_module)

    sweep_params = {}

    sweep_params["steam_cost"] = LinearSample(
        m.fs.costing.heat_exchanger.steam_cost, 0.00, 0.008, nx
    )
    # sweep_params = {}
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
        number_of_subprocesses=8
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    fc_results, sweep_params, m = run_cryst_sweep(steam, output_filename="FC_cryst_sweep-big.csv")
    # clear_output(wait=False)
    mvc_results, sweep_params, m_mvc = run_cryst_sweep(mvc, output_filename="MVC_cryst_sweep-big.csv")
    # clear_output(wait=False)
    tvc_results, sweep_params, m_tvc = run_cryst_sweep(tvc, output_filename="TVC_cryst_sweep-big.csv")
# clear_output(wait=False)