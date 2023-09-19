from idaes.core.solvers import get_solver
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    build,
    set_operating_conditions,
    initialize_system,
    optimize,
)
from watertap.examples.flowsheets.RO_with_energy_recovery.monte_carlo_sampling_RO_ERD import (
    get_sweep_params,
    build_model,
    # build_sweep_params,
    build_outputs,
    run_parameter_sweep,
)
from watertap.tools.parameter_sweep import (
    LinearSample,
    ParameterSweep,
)


def build_sweep_params(m, num_samples=1, scenario="A_comp_vs_LCOW"):
    sweep_params = {}

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
            0.01, 0.05, num_samples
        )
    else:
        raise NotImplementedError
    
    return sweep_params


def run_parameter_sweep(num_samples=100, num_procs=1):

    # solver = get_solver()
    ps = create_parameter_sweep_object(solver, num_samples, num_procs)
    results_dict, results_array = ps.parameter_sweep(
        build_model,
        build_sweep_params,
        build_outputs=None,
        build_outputs_kwargs=None,
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=None,
        build_sweep_params_kwargs=None,
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
            num_samples=num_samples,
            scenario="A_comp_vs_B_comp_vs_LCOW"
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
        "parallel_back_end": parallel_backend, # "MultiProcessing",
        "log_model_states": False,
    }
    ps = ParameterSweep(**kwargs_dict)
    return ps, kwargs_dict


if __name__ == "__main__":

    import sys
    import time

    start_time = time.time()

    if len(sys.argv) == 1:
        num_samples = 2
        num_procs = 1
    else:
        num_samples = int(sys.argv[1])
        num_procs = int(sys.argv[2])

    ps, kwargs_dict = create_parameter_sweep_object(num_samples, num_procs)

    results_array, results_dict = ps.parameter_sweep(
        kwargs_dict["build_model"],
        kwargs_dict["build_sweep_params"],
        build_outputs=kwargs_dict["build_outputs"],
        build_outputs_kwargs=kwargs_dict["build_outputs_kwargs"],
        num_samples=num_samples,
        seed=None,
        build_model_kwargs=kwargs_dict["build_model_kwargs"],
        build_sweep_params_kwargs=kwargs_dict["build_sweep_params_kwargs"],
    )

    end_time = time.time()
    time_elapsed = end_time - start_time
    print("time_elapsed = ", time_elapsed)

    # import pprint
    # pprint.pprint(results_dict)
