#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

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
    optimize,
)

from watertap.tools.parameter_sweep import (
    get_sweep_params_from_yaml,
    set_defaults_from_yaml,
)


def get_sweep_params(m, num_samples, use_LHS=False):
    sweep_params = {}

    # Define the sampling type and ranges for three different variables
    if use_LHS:
        sweep_params["A_comp"] = LatinHypercubeSample(
            m.fs.RO.A_comp, 4.0e-12, 0.5e-12, num_samples
        )
        sweep_params["B_comp"] = LatinHypercubeSample(
            m.fs.RO.B_comp, 3.5e-8, 0.5e-8, num_samples
        )
        sweep_params["Spacer_porosity"] = LatinHypercubeSample(
            m.fs.RO.feed_side.spacer_porosity, 0.95, 0.99, num_samples
        )

    else:
        sweep_params["A_comp"] = NormalSample(
            m.fs.RO.A_comp, 4.0e-12, 0.5e-12, num_samples
        )
        sweep_params["B_comp"] = NormalSample(
            m.fs.RO.B_comp, 3.5e-8, 0.5e-8, num_samples
        )
        sweep_params["Spacer_porosity"] = UniformSample(
            m.fs.RO.feed_side.spacer_porosity, 0.95, 0.99, num_samples
        )

    return sweep_params


def build_model(
    read_model_defauls_from_file=False,
    defaults_fname="default_configuration.yaml",
):
    # Set up the solver
    solver = get_solver()

    # Build, set, and initialize the system (these steps will change depending on the underlying model)
    m = build()
    set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
    initialize_system(m, solver=solver)

    # Check if we need to read in the default model values from a file
    if read_model_defauls_from_file:
        set_defaults_from_yaml(m, defaults_fname)

    return m


def build_sweep_params(
    m,
    use_LHS=False,
    sweep_params_fname="mc_sweep_params.yaml",
    read_sweep_params_from_file=False,
    num_samples=10,
):
    # Define the sampling type and ranges for three different variables
    if read_sweep_params_from_file:
        sweep_params = get_sweep_params_from_yaml(m, sweep_params_fname)
    else:
        sweep_params = get_sweep_params(m, num_samples, use_LHS=use_LHS)

    return sweep_params


def build_outputs(m):
    outputs = {}
    outputs["EC"] = m.fs.costing.specific_energy_consumption
    outputs["LCOW"] = m.fs.costing.LCOW
    return outputs


def run_parameter_sweep(
    csv_results_file_name=None,
    h5_results_file_name=None,
    seed=None,
    use_LHS=False,
    sweep_params_fname="mc_sweep_params.yaml",
    read_sweep_params_from_file=False,
    read_model_defauls_from_file=False,
    defaults_fname="default_configuration.yaml",
):
    # Set up the solver
    solver = get_solver()

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10

    # Run the parameter sweep
    global_results_arr, _ = parameter_sweep(
        build_model,
        build_sweep_params,
        build_outputs,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        optimize_function=optimize,
        optimize_kwargs={"solver": solver, "check_termination": False},
        num_samples=num_samples,
        seed=seed,
        build_model_kwargs=dict(
            read_model_defauls_from_file=read_model_defauls_from_file,
            defaults_fname=defaults_fname,
        ),
        build_sweep_params_kwargs=dict(
            num_samples=num_samples,
            use_LHS=use_LHS,
            sweep_params_fname=sweep_params_fname,
            read_sweep_params_from_file=read_sweep_params_from_file,
        ),
        number_of_subprocesses=1,
    )

    return global_results_arr


if __name__ == "__main__":
    # For testing this file, a seed needs to be provided as an additional argument, i.e. seed=1
    run_parameter_sweep("output/monte_carlo_results.csv")
