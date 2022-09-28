###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from idaes.core.solvers import get_solver
from watertap.tools.parameter_sweep import (
    UniformSample,
    NormalSample,
    LatinHypercubeSample,
    parameter_sweep,
)

from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    optimize,
)

from watertap.tools.parameter_sweep import (
    get_sweep_params_from_yaml,
    set_defaults_from_yaml,
)


def get_sweep_params(m, use_LHS=False):
    sweep_params = {}

    # Define the sampling type and ranges for three different variables
    if use_LHS:
        sweep_params["A_comp"] = LatinHypercubeSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
        sweep_params["B_comp"] = LatinHypercubeSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)
        sweep_params["Spacer_porosity"] = LatinHypercubeSample(
            m.fs.RO.spacer_porosity, 0.95, 0.99
        )

    else:
        sweep_params["A_comp"] = NormalSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
        sweep_params["B_comp"] = NormalSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)
        sweep_params["Spacer_porosity"] = UniformSample(
            m.fs.RO.spacer_porosity, 0.95, 0.99
        )

    return sweep_params


def run_parameter_sweep(
    csv_results_file_name=None,
    h5_results_file_name=None,
    seed=None,
    use_LHS=False,
    read_sweep_params_from_file=False,
    sweep_params_fname="mc_sweep_params.yaml",
    read_model_defauls_from_file=False,
    defaults_fname="default_configuration.yaml",
):

    # Set up the solver
    solver = get_solver()

    # Build, set, and initialize the system (these steps will change depending on the underlying model)
    m = build()
    set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
    initialize_system(m, solver=solver)

    # Simulate once outside the parameter sweep to ensure everything is appropriately initialized
    solve(m, solver=solver)

    # Check if we need to read in the default model values from a file
    if read_model_defauls_from_file:
        set_defaults_from_yaml(m, defaults_fname)

    # Define the sampling type and ranges for three different variables
    if read_sweep_params_from_file:
        sweep_params = get_sweep_params_from_yaml(m, sweep_params_fname)
    else:
        sweep_params = get_sweep_params(m, use_LHS=use_LHS)

    # Define the outputs to be saved
    outputs = {}
    outputs["EC"] = m.fs.costing.specific_energy_consumption
    outputs["LCOW"] = m.fs.costing.LCOW

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10

    # Run the parameter sweep
    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        optimize_function=optimize,
        optimize_kwargs={"solver": solver, "check_termination": False},
        num_samples=num_samples,
        seed=seed,
    )

    return global_results


if __name__ == "__main__":
    # For testing this file, a seed needs to be provided as an additional argument, i.e. seed=1
    run_parameter_sweep("output/monte_carlo_results.csv")
