from idaes.core.solvers import get_solver
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    optimize,
)
from watertap.examples.flowsheets.RO_with_energy_recovery.monte_carlo_sampling_RO_ERD import (
    build_model,
    build_outputs,
)
from watertap.tools.parameter_sweep import (
    LinearSample,
    UniformSample,
    ParameterSweep,
    RecursiveParameterSweep,
    DifferentialParameterSweep,
)


def build_sweep_params(m, num_samples=1, scenario="A_comp_vs_LCOW"):
    sweep_params = {}

    # A_comp: Membrane water permeability coefficient
    # B_comp: Membrane salt permeability coefficient
    # NaCl_Loading: Feed water salinity
    # LCOW: Levelized cost of water

    if scenario == "A_comp_vs_LCOW":
        sweep_params["A_comp"] = LinearSample(
            m.fs.RO.A_comp, 1.0e-12, 1e-11, num_samples
        )
    elif scenario == "A_comp_vs_B_comp_vs_LCOW":
        sweep_params["A_comp"] = LinearSample(
            m.fs.RO.A_comp, 1.0e-12, 1e-11, num_samples
        )
        sweep_params["B_comp"] = LinearSample(
            m.fs.RO.B_comp, 8.0e-8, 1.0e-8, num_samples
        )
    elif scenario == "WR_vs_NaCL_loading_vs_LCOW":
        sweep_params["recovery"] = LinearSample(
            m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"], 0.1, 0.65, num_samples
        )
        sweep_params["NaCl_loading"] = LinearSample(
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"],
            0.01,
            0.05,
            num_samples,
        )
    else:
        raise NotImplementedError

    return sweep_params


def run_parameter_sweep(num_samples=100, num_procs=1):
    ps1, kwargs_dict1 = create_parameter_sweep_object(
        num_samples, num_procs, parallel_backend="ConcurrentFutures"
    )

    results_array1, results_dict1 = ps1.parameter_sweep(
        kwargs_dict1["build_model"],
        kwargs_dict1["build_sweep_params"],
        build_outputs=kwargs_dict1["build_outputs"],
        build_outputs_kwargs=kwargs_dict1["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict1["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict1["build_sweep_params_kwargs"],
    )


def run_recursive_parameter_sweep(num_samples=100, num_procs=1):
    ps2, kwargs_dict2 = create_recursive_parameter_sweep_object(num_samples, num_procs)
    results_array2, results_dict2 = ps2.parameter_sweep(
        kwargs_dict2["build_model"],
        kwargs_dict2["build_sweep_params"],
        build_outputs=kwargs_dict2["build_outputs"],
        build_outputs_kwargs=kwargs_dict2["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict2["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict2["build_sweep_params_kwargs"],
    )


def run_differential_parameter_sweep(num_samples=10, num_procs=1):
    model, ps3, kwargs_dict3 = create_differential_parameter_sweep_object(
        num_samples, num_procs
    )
    results_array3, results_dict3 = ps3.parameter_sweep(
        kwargs_dict3["build_model"],
        kwargs_dict3["build_sweep_params"],
        build_outputs=kwargs_dict3["build_outputs"],
        build_outputs_kwargs=kwargs_dict3["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict3["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict3["build_sweep_params_kwargs"],
    )


def create_parameter_sweep_object(
    num_samples,
    num_procs,
    parallel_backend="ConcurrentFutures",
):
    solver = get_solver()
    kwargs_dict = {
        "debugging_data_dir": None,
        "csv_results_file_name": None,
        "h5_results_file_name": None,
        "interpolate_nan_outputs": False,
        "h5_parent_group_name": None,  # Loop Tool
        "build_model": build_model,
        "build_model_kwargs": dict(
            read_model_defauls_from_file=False,
            defaults_fname="default_configuration.yaml",
        ),
        "build_sweep_params": build_sweep_params,
        "build_sweep_params_kwargs": dict(
            num_samples=num_samples, scenario="A_comp_vs_B_comp_vs_LCOW"
        ),
        "build_outputs": build_outputs,
        "build_outputs_kwargs": {},
        "optimize_function": optimize,
        "optimize_kwargs": {"solver": solver, "check_termination": False},
        "initialize_function": None,
        "update_sweep_params_before_init": False,
        "initialize_kwargs": {},
        "initialize_before_sweep": False,
        "reinitialize_function": None,
        "reinitialize_kwargs": {},
        "reinitialize_before_sweep": False,
        "probe_function": None,
        "custom_do_param_sweep": None,
        "custom_do_param_sweep_kwargs": {},
        "publish_progress": False,
        "publish_address": "http://localhost:8888",
        "number_of_subprocesses": num_procs,
        "parallel_back_end": parallel_backend,  # "MultiProcessing",
        "log_model_states": False,
    }
    ps = ParameterSweep(**kwargs_dict)
    return ps, kwargs_dict


def create_recursive_parameter_sweep_object(
    num_samples,
    num_procs,
    parallel_backend="ConcurrentFutures",
):
    solver = get_solver()
    kwargs_dict = {
        "debugging_data_dir": None,
        "csv_results_file_name": None,
        "h5_results_file_name": None,
        "interpolate_nan_outputs": False,
        "h5_parent_group_name": None,  # Loop Tool
        "build_model": build_model,
        "build_model_kwargs": dict(
            read_model_defauls_from_file=False,
            defaults_fname="default_configuration.yaml",
        ),
        "build_sweep_params": build_sweep_params,
        "build_sweep_params_kwargs": dict(
            num_samples=num_samples, scenario="A_comp_vs_B_comp_vs_LCOW"
        ),
        "build_outputs": build_outputs,
        "build_outputs_kwargs": {},
        "optimize_function": optimize,
        "optimize_kwargs": {"solver": solver, "check_termination": False},
        "initialize_function": None,
        "update_sweep_params_before_init": False,
        "initialize_kwargs": {},
        "initialize_before_sweep": False,
        "reinitialize_function": None,
        "reinitialize_kwargs": {},
        "reinitialize_before_sweep": False,
        "probe_function": None,
        "custom_do_param_sweep": None,
        "custom_do_param_sweep_kwargs": {},
        "publish_progress": False,
        "publish_address": "http://localhost:8888",
        "number_of_subprocesses": num_procs,
        "parallel_back_end": parallel_backend,  # "MultiProcessing",
        "log_model_states": False,
    }
    ps = RecursiveParameterSweep(**kwargs_dict)
    return ps, kwargs_dict


def create_differential_parameter_sweep_object(
    num_samples,
    num_procs,
    parallel_backend="ConcurrentFutures",
):
    solver = get_solver()
    m = build_model(read_model_defauls_from_file=False)

    def build_spec(model):
        differential_sweep_specs = {
            "A_comp": {
                "diff_mode": "sum",
                "diff_sample_type": UniformSample,
                "relative_lb": 0.01,
                "relative_ub": 0.01,
                "pyomo_object": model.fs.RO.A_comp,
            }
        }
        return differential_sweep_specs

    kwargs_dict = {
        "debugging_data_dir": None,
        "csv_results_file_name": None,
        "h5_results_file_name": None,
        "interpolate_nan_outputs": False,
        "h5_parent_group_name": None,  # Loop Tool
        "build_model": build_model,
        "build_model_kwargs": dict(
            read_model_defauls_from_file=False,
            defaults_fname="default_configuration.yaml",
        ),
        "build_sweep_params": build_sweep_params,
        "build_sweep_params_kwargs": dict(
            num_samples=num_samples, scenario="A_comp_vs_B_comp_vs_LCOW"
        ),
        "build_outputs": build_outputs,
        "build_outputs_kwargs": {},
        "optimize_function": optimize,
        "optimize_kwargs": {"solver": solver, "check_termination": False},
        "num_diff_samples": 1,
        "build_differential_sweep_specs": build_spec,
        "initialize_function": None,
        "update_sweep_params_before_init": False,
        "initialize_kwargs": {},
        "initialize_before_sweep": False,
        "reinitialize_function": None,
        "reinitialize_kwargs": {},
        "reinitialize_before_sweep": False,
        "probe_function": None,
        "custom_do_param_sweep": None,
        "custom_do_param_sweep_kwargs": {},
        "publish_progress": False,
        "publish_address": "http://localhost:8888",
        "number_of_subprocesses": num_procs,
        "parallel_back_end": parallel_backend,  # "MultiProcessing",
        "log_model_states": False,
    }
    ps = DifferentialParameterSweep(**kwargs_dict)

    return m, ps, kwargs_dict


if __name__ == "__main__":
    import sys
    import time
    import numpy as np
    import pprint

    start_time = time.time()

    if len(sys.argv) == 1:
        num_samples = 2
        num_procs = 1
    else:
        num_samples = int(sys.argv[1])
        num_procs = int(sys.argv[2])

    ps1, kwargs_dict1 = create_parameter_sweep_object(
        num_samples, num_procs, parallel_backend="ConcurrentFutures"
    )

    results_array1, results_dict1 = ps1.parameter_sweep(
        kwargs_dict1["build_model"],
        kwargs_dict1["build_sweep_params"],
        build_outputs=kwargs_dict1["build_outputs"],
        build_outputs_kwargs=kwargs_dict1["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict1["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict1["build_sweep_params_kwargs"],
    )

    ps2, kwargs_dict2 = create_recursive_parameter_sweep_object(num_samples, num_procs)
    results_array2, results_dict2 = ps2.parameter_sweep(
        kwargs_dict2["build_model"],
        kwargs_dict2["build_sweep_params"],
        build_outputs=kwargs_dict2["build_outputs"],
        build_outputs_kwargs=kwargs_dict2["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict2["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict2["build_sweep_params_kwargs"],
    )

    model, ps3, kwargs_dict3 = create_differential_parameter_sweep_object(
        num_samples, num_procs
    )
    results_array3, results_dict3 = ps3.parameter_sweep(
        kwargs_dict3["build_model"],
        kwargs_dict3["build_sweep_params"],
        build_outputs=kwargs_dict3["build_outputs"],
        build_outputs_kwargs=kwargs_dict3["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict3["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict3["build_sweep_params_kwargs"],
    )

    end_time = time.time()
    time_elapsed = end_time - start_time
    print("time_elapsed = ", time_elapsed)
