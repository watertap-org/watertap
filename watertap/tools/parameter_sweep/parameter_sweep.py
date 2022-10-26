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
import numpy as np
import pyomo.environ as pyo
import warnings
import copy, pprint

from abc import abstractmethod, ABC
from idaes.core.solvers import get_solver

from idaes.core.surrogate.pysmo import sampling
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.tee import capture_output
from pyomo.common.config import ConfigValue

from watertap.tools.parameter_sweep.parameter_sweep_writer import ParameterSweepWriter
from watertap.tools.parameter_sweep.sampling_types import SamplingType, LinearSample

import watertap.tools.MPI as MPI


def _default_optimize(model, options=None, tee=False):
    """
    Default optimization function used in parameter_sweep.
    Optimizes ``model`` using the IDAES default solver.
    Raises a RuntimeError if the TerminationCondition is not optimal

    Arguments:

        model : A Pyomo ConcreteModel to optimize

        options (optional) : Solver options to pass into idaes.core.utils.get_solver.
                             Default is None
        tee (options) : To display the solver log. Default it False

    """
    solver = get_solver(options=options)
    results = solver.solve(model, tee=tee)
    return results


class _ParameterSweepBase(ABC):

    CONFIG = ParameterSweepWriter.CONFIG()

    CONFIG.declare(
        "optimize_function",
        ConfigValue(
            default=_default_optimize,
            # domain=function,
            description="Optimization function to be used for the parameter sweep.",
        ),
    )

    CONFIG.declare(
        "optimize_kwargs",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Keyword argument for the optimization function for the parameter sweep.",
        ),
    )

    CONFIG.declare(
        "reinitialize_function",
        ConfigValue(
            default=None,
            # domain=function,
            description="Function to reinitialize a flowsheet",
        ),
    )

    CONFIG.declare(
        "reinitialize_kwargs",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Keyword arguments for the reinitialization function.",
        ),
    )

    CONFIG.declare(
        "reinitialize_before_sweep",
        ConfigValue(
            default=False,
            domain=bool,
            description="Reinitializing a model before every iteration.",
        ),
    )

    CONFIG.declare(
        "probe_function",
        ConfigValue(
            default=None,
            # domain=function,
            description="Function to probe if a flowsheet configuration will work",
        ),
    )

    def __init__(
        self,
        **options,
    ):

        self.comm = options.pop("comm", MPI.COMM_WORLD)
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        self.config = self.CONFIG(options)

        # Initialize the writer
        self.writer = ParameterSweepWriter(
            self.comm,
            csv_results_file_name=self.config.csv_results_file_name,
            h5_results_file_name=self.config.h5_results_file_name,
            debugging_data_dir=self.config.debugging_data_dir,
            interpolate_nan_outputs=self.config.interpolate_nan_outputs,
        )

    def _build_combinations(self, d, sampling_type, num_samples):
        num_var_params = len(d)

        if self.rank == 0:
            param_values = []

            for k, v in d.items():
                # Build a vector of discrete values for this parameter
                p = v.sample(num_samples)
                param_values.append(p)

            if sampling_type == SamplingType.FIXED:
                # Form an array with every possible combination of parameter values
                global_combo_array = np.array(np.meshgrid(*param_values, indexing="ij"))
                global_combo_array = global_combo_array.reshape(num_var_params, -1).T

            elif sampling_type == SamplingType.RANDOM:
                sorting = np.argsort(param_values[0])
                global_combo_array = np.vstack(param_values).T
                global_combo_array = global_combo_array[sorting, :]

            elif sampling_type == SamplingType.RANDOM_LHS:
                lb = [val[0] for val in param_values]
                ub = [val[1] for val in param_values]
                lhs = sampling.LatinHypercubeSampling(
                    [lb, ub], number_of_samples=num_samples, sampling_type="creation"
                )
                global_combo_array = lhs.sample_points()
                sorting = np.argsort(global_combo_array[:, 0])
                global_combo_array = global_combo_array[sorting, :]

            else:
                raise ValueError(f"Unknown sampling type: {sampling_type}")

            # Test if the global_combo_array is in row-major order
            if not global_combo_array.flags.c_contiguous:
                # If not, return a copy of this array with row-major memory order
                global_combo_array = np.ascontiguousarray(global_combo_array)

        else:
            if sampling_type == SamplingType.FIXED:
                nx = 1
                for k, v in d.items():
                    nx *= v.num_samples
            elif (
                sampling_type == SamplingType.RANDOM
                or sampling_type == SamplingType.RANDOM_LHS
            ):
                nx = num_samples
            else:
                raise ValueError(f"Unknown sampling type: {sampling_type}")

            if not float(nx).is_integer():
                raise RuntimeError(f"Total number of samples must be integer valued")
            nx = int(nx)

            # Allocate memory to hold the Bcast array
            global_combo_array = np.zeros((nx, num_var_params), dtype=np.float64)

        ### Broadcast the array to all processes
        if self.num_procs > 1:
            self.comm.Bcast(global_combo_array, root=0)

        return global_combo_array

    def _divide_combinations(self, global_combo_array):

        # Split the total list of combinations into NUM_PROCS chunks,
        # one per each of the MPI ranks
        # divided_combo_array = np.array_split(global_combo_array, num_procs, axis=0)
        divided_combo_array = np.array_split(global_combo_array, self.num_procs)

        # Return only this rank's portion of the total workload
        local_combo_array = divided_combo_array[self.rank]

        return local_combo_array

    def _update_model_values(self, m, param_dict, values):

        for k, item in enumerate(param_dict.values()):

            param = item.pyomo_object

            if param.is_variable_type():
                # Fix the single value to values[k]
                param.fix(values[k])

            elif param.is_parameter_type():
                # Fix the single value to values[k]
                param.set_value(values[k])

            else:
                raise RuntimeError(f"Unrecognized Pyomo object {param}")

    def _aggregate_results_arr(self, global_results_dict, num_cases):

        global_results = np.zeros(
            (num_cases, len(global_results_dict["outputs"])), dtype=np.float64
        )

        if self.rank == 0:
            for i, (key, item) in enumerate(global_results_dict["outputs"].items()):
                global_results[:, i] = item["value"][:num_cases]

        if self.num_procs > 1:  # pragma: no cover
            self.comm.Bcast(global_results, root=0)

        return global_results

    def _process_sweep_params(self, sweep_params):

        sampling_type = None

        # Check the list of parameters to make sure they are valid
        for k in sweep_params:

            # Convert to using Sample class
            if isinstance(sweep_params[k], (list, tuple)):
                sweep_params[k] = LinearSample(*sweep_params[k])

            # Get the type of sampling
            current_sampling_type = sweep_params[k].sampling_type

            # Check to make sure only one sampling type is provided
            if sampling_type is None:
                sampling_type = current_sampling_type
            elif current_sampling_type != sampling_type:
                raise ValueError("Cannot mix sampling types")

        return sweep_params, sampling_type

    def _create_local_output_skeleton(self, model, sweep_params, outputs, num_samples):

        output_dict = {}
        output_dict["sweep_params"] = {}
        output_dict["outputs"] = {}

        sweep_param_objs = ComponentSet()

        # Store the inputs
        for sweep_param in sweep_params.values():
            var = sweep_param.pyomo_object
            sweep_param_objs.add(var)
            output_dict["sweep_params"][
                var.name
            ] = self._create_component_output_skeleton(var, num_samples)

        if outputs is None:
            # No outputs are specified, so every Var, Expression, and Objective on the model should be saved
            for pyo_obj in model.component_data_objects(
                (pyo.Var, pyo.Expression, pyo.Objective), active=True
            ):
                # Only need to save this variable if it isn't one of the value in sweep_params
                if pyo_obj not in sweep_param_objs:
                    output_dict["outputs"][
                        pyo_obj.name
                    ] = self._create_component_output_skeleton(pyo_obj, num_samples)

        else:
            # Save only the outputs specified in the outputs dictionary
            for short_name, pyo_obj in outputs.items():
                output_dict["outputs"][
                    short_name
                ] = self._create_component_output_skeleton(pyo_obj, num_samples)

        return output_dict

    def _create_component_output_skeleton(self, component, num_samples):

        comp_dict = {}
        comp_dict["value"] = np.zeros(num_samples, dtype=np.float64)
        if hasattr(component, "lb"):
            comp_dict["lower bound"] = component.lb
        if hasattr(component, "ub"):
            comp_dict["upper bound"] = component.ub
        if hasattr(component, "get_units"):
            unit_obj = component.get_units()
            if unit_obj is not None:
                comp_dict["units"] = component.get_units().name
            else:
                comp_dict["units"] = "None"

        # Add information to this output that WILL NOT be written as part
        # of the file saving step.
        comp_dict["_pyo_obj"] = component

        return comp_dict

    def _update_local_output_dict(
        self, model, sweep_params, case_number, sweep_vals, run_successful, output_dict
    ):

        # Get the inputs
        op_ps_dict = output_dict["sweep_params"]
        for key, item in sweep_params.items():
            var_name = item.pyomo_object.name
            op_ps_dict[var_name]["value"][case_number] = item.pyomo_object.value

        # Get the outputs from model
        if run_successful:
            for label, val in output_dict["outputs"].items():
                output_dict["outputs"][label]["value"][case_number] = pyo.value(
                    val["_pyo_obj"]
                )

        else:
            for label in output_dict["outputs"].keys():
                output_dict["outputs"][label]["value"][case_number] = np.nan

    def _create_global_output(self, local_output_dict, req_num_samples):

        # Before we can create the global dictionary, we need to delete the pyomo
        # object contained within the dictionary
        for key, val in local_output_dict.items():
            if key != "solve_successful":
                for subval in val.values():
                    if "_pyo_obj" in subval:
                        del subval["_pyo_obj"]

        # We make the assumption that the parameter sweep is running the same
        # flowsheet num_samples number of times, i.e., the structure of the
        # local_output_dict remains the same across all mpi_ranks
        local_num_cases = len(local_output_dict["solve_successful"])

        # Gather the size of the value array on each MPI rank
        sample_split_arr = self.comm.allgather(local_num_cases)
        num_total_samples = sum(sample_split_arr)

        # Create the global value array on rank 0
        if self.rank == 0:
            global_output_dict = copy.deepcopy(local_output_dict)
            # Create a global value array of inputs in the dictionary
            for key, item in global_output_dict.items():
                if key != "solve_successful":
                    for subkey, subitem in item.items():
                        subitem["value"] = np.zeros(num_total_samples, dtype=np.float64)

        else:
            global_output_dict = local_output_dict

        # Finally collect the values
        for key, item in local_output_dict.items():
            if key != "solve_successful":
                for subkey, subitem in item.items():
                    self.comm.Gatherv(
                        sendbuf=subitem["value"],
                        recvbuf=(
                            global_output_dict[key][subkey]["value"],
                            sample_split_arr,
                        ),
                        root=0,
                    )

                    # Trim to the exact number
                    global_output_dict[key][subkey]["value"] = global_output_dict[key][
                        subkey
                    ]["value"][0:req_num_samples]

            elif key == "solve_successful":
                local_solve_successful = np.fromiter(
                    item, dtype=np.bool, count=len(item)
                )

                if self.rank == 0:
                    global_solve_successful = np.empty(num_total_samples, dtype=np.bool)
                else:
                    global_solve_successful = None

                self.comm.Gatherv(
                    sendbuf=local_solve_successful,
                    recvbuf=(global_solve_successful, sample_split_arr),
                    root=0,
                )

                if self.rank == 0:
                    global_output_dict[key] = global_solve_successful[0:req_num_samples]

        return global_output_dict

    def _param_sweep_kernel(self, model, reinitialize_values):

        optimize_function = self.config.optimize_function
        optimize_kwargs = self.config.optimize_kwargs
        reinitialize_before_sweep = self.config.reinitialize_before_sweep
        reinitialize_function = self.config.reinitialize_function
        reinitialize_kwargs = self.config.reinitialize_kwargs

        run_successful = False  # until proven otherwise

        # Forced reinitialization of the flowsheet if enabled
        if reinitialize_before_sweep:
            if reinitialize_function is None:
                raise ValueError(
                    "Reinitialization function was not specified. The model will not be reinitialized."
                )
            else:
                for v, val in reinitialize_values.items():
                    if not v.fixed:
                        v.set_value(val, skip_validation=True)
                reinitialize_function(model, **reinitialize_kwargs)

        try:
            # Simulate/optimize with this set of parameter
            with capture_output():
                results = optimize_function(model, **optimize_kwargs)
            pyo.assert_optimal_termination(results)

        except:
            # run_successful remains false. We try to reinitialize and solve again
            if reinitialize_function is not None:
                for v, val in reinitialize_values.items():
                    if not v.fixed:
                        v.set_value(val, skip_validation=True)
                try:
                    reinitialize_function(model, **reinitialize_kwargs)
                    with capture_output():
                        results = optimize_function(model, **optimize_kwargs)
                    pyo.assert_optimal_termination(results)

                except:
                    pass  # run_successful is still False
                else:
                    run_successful = True

        else:
            # If the simulation suceeds, report stats
            run_successful = True

        return run_successful

    def _do_param_sweep(self, model, sweep_params, outputs, local_values):

        # Create easy to read variables for configurations
        probe_function = self.config["probe_function"]

        # Initialize space to hold results
        local_num_cases = np.shape(local_values)[0]

        # Create the output skeleton for storing detailed data
        local_output_dict = self._create_local_output_skeleton(
            model, sweep_params, outputs, local_num_cases
        )

        local_results = np.zeros((local_num_cases, len(local_output_dict["outputs"])))

        local_solve_successful_list = []

        if self.config["reinitialize_function"] is not None:
            reinitialize_values = ComponentMap()
            for v in model.component_data_objects(pyo.Var):
                reinitialize_values[v] = v.value
        else:
            reinitialize_values = None

        # ================================================================
        # Run all optimization cases
        # ================================================================

        for k in range(local_num_cases):
            # Update the model values with a single combination from the parameter space
            self._update_model_values(model, sweep_params, local_values[k, :])

            if probe_function is None or probe_function(model):
                run_successful = self._param_sweep_kernel(
                    model,
                    reinitialize_values,
                )
            else:
                run_successful = False

            # Update the loop based on the reinitialization
            self._update_local_output_dict(
                model,
                sweep_params,
                k,
                local_values[k, :],
                run_successful,
                local_output_dict,
            )

            local_solve_successful_list.append(run_successful)

        local_output_dict["solve_successful"] = local_solve_successful_list

        return local_output_dict

    @abstractmethod
    def parameter_sweep(self, *args, **kwargs):
        pass


class ParameterSweep(_ParameterSweepBase):

    CONFIG = _ParameterSweepBase.CONFIG()

    def _aggregate_local_results(
        self, global_values, local_output_dict, num_samples, local_num_cases
    ):

        # Create the dictionary
        global_results_dict = self._create_global_output(local_output_dict, num_samples)

        # Create the array
        num_global_samples = np.shape(global_values)[0]
        global_results_arr = self._aggregate_results_arr(
            global_results_dict, num_global_samples
        )

        return global_results_dict, global_results_arr

    def parameter_sweep(
        self,
        model,
        sweep_params,
        outputs=None,
        num_samples=None,
        seed=None,
    ):

        # Convert sweep_params to LinearSamples
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)

        # Set the seed before sampling
        np.random.seed(seed)

        # Enumerate/Sample the parameter space
        global_values = self._build_combinations(
            sweep_params, sampling_type, num_samples
        )

        # divide the workload between processors
        local_values = self._divide_combinations(global_values)
        local_num_cases = np.shape(local_values)[0]

        # Do the Loop
        local_results_dict = self._do_param_sweep(
            model,
            sweep_params,
            outputs,
            local_values,
        )

        # Aggregate results on Master
        global_results_dict, global_results_arr = self._aggregate_local_results(
            global_values, local_results_dict, num_samples, local_num_cases
        )

        # Save to file
        global_save_data = self.writer.save_results(
            sweep_params,
            local_values,
            global_values,
            local_results_dict,
            global_results_dict,
            global_results_arr,
        )

        return global_save_data


class RecursiveParameterSweep(_ParameterSweepBase):

    CONFIG = _ParameterSweepBase.CONFIG()

    def _filter_recursive_solves(
        self, model, sweep_params, outputs, recursive_local_dict
    ):

        # Figure out how many filtered solves did this rank actually do
        filter_counter = 0
        for case, content in recursive_local_dict.items():
            filter_counter += sum(content["solve_successful"])

        # Now that we have all of the local output dictionaries, we need to construct
        # a consolidated dictionary of successful solves.
        local_filtered_dict = self._create_local_output_skeleton(
            model, sweep_params, outputs, filter_counter
        )
        local_filtered_dict["solve_successful"] = []

        # Populate local_successful_outputs
        offset = 0
        for case_number, content in recursive_local_dict.items():
            # Filter all of the sucessful solves
            optimal_indices = [
                idx
                for idx, success in enumerate(content["solve_successful"])
                if success
            ]
            n_successful_solves = len(optimal_indices)
            stop = offset + n_successful_solves

            for key, item in content.items():
                if key != "solve_successful":
                    for subkey, subitem in item.items():
                        local_filtered_dict[key][subkey]["value"][
                            offset:stop
                        ] = subitem["value"][optimal_indices]

            # Place the solve status
            local_filtered_dict["solve_successful"].extend(
                [content["solve_successful"][i] for i in optimal_indices]
            )

            offset += n_successful_solves

        return local_filtered_dict, filter_counter

    def _aggregate_filtered_input_arr(self, global_filtered_dict, req_num_samples):

        global_filtered_values = np.zeros(
            (req_num_samples, len(global_filtered_dict["sweep_params"])),
            dtype=np.float64,
        )

        if self.rank == 0:
            for i, (key, item) in enumerate(
                global_filtered_dict["sweep_params"].items()
            ):
                global_filtered_values[:, i] = item["value"][:req_num_samples]

        if self.num_procs > 1:  # pragma: no cover
            self.comm.Bcast(global_filtered_values, root=0)

        return global_filtered_values

    def _aggregate_filtered_results(self, local_filtered_dict, req_num_samples):

        global_filtered_dict = self._create_global_output(
            local_filtered_dict, req_num_samples
        )
        global_filtered_results = self._aggregate_results_arr(
            global_filtered_dict, req_num_samples
        )
        global_filtered_values = self._aggregate_filtered_input_arr(
            global_filtered_dict, req_num_samples
        )

        return global_filtered_dict, global_filtered_results, global_filtered_values

    def parameter_sweep(
        self,
        model,
        sweep_params,
        outputs=None,
        req_num_samples=None,
        seed=None,
    ):

        # Convert sweep_params to LinearSamples
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)

        # Set the seed before sampling
        np.random.seed(seed)

        n_samples_remaining = req_num_samples
        num_total_samples = req_num_samples

        local_output_collection = {}
        loop_ctr = 0
        while n_samples_remaining > 0 and loop_ctr < 10:
            # Enumerate/Sample the parameter space
            global_values = self._build_combinations(
                sweep_params, sampling_type, num_total_samples
            )

            # divide the workload between processors
            local_values = self._divide_combinations(global_values)
            local_num_cases = np.shape(local_values)[0]
            if loop_ctr == 0:
                true_local_num_cases = local_num_cases

            local_output_collection[loop_ctr] = self._do_param_sweep(
                model,
                sweep_params,
                outputs,
                local_values,
            )

            # Get the number of successful solves on this proc (sum of boolean flags)
            success_count = sum(local_output_collection[loop_ctr]["solve_successful"])
            failure_count = local_num_cases - success_count

            # Get the global number of successful solves and update the number of remaining samples
            if self.num_procs > 1:  # pragma: no cover
                global_success_count = np.zeros(1, dtype=np.float64)
                global_failure_count = np.zeros(1, dtype=np.float64)
                self.comm.Allreduce(
                    np.array(success_count, dtype=np.float64), global_success_count
                )
                self.comm.Allreduce(
                    np.array(failure_count, dtype=np.float64), global_failure_count
                )
            else:
                global_success_count = success_count
                global_failure_count = failure_count

            success_prob = global_success_count / (
                global_failure_count + global_success_count
            )

            if success_prob < 0.1:
                warnings.warn(
                    f"Success rate of solves = {100.0*success_prob}%, consider adjusting sweep limits."
                )

            n_samples_remaining -= global_success_count

            # The total number of samples to generate at the next iteration is a multiple of the total remaining samples
            scale_factor = 2.0 / max(success_prob, 0.10)
            num_total_samples = int(np.ceil(scale_factor * n_samples_remaining))
            loop_ctr += 1

        # Now that we have all of the local output dictionaries, we need to construct
        # a consolidated dictionary based on a filter, e.g., optimal solves.
        local_filtered_dict, local_n_successful = self._filter_recursive_solves(
            model, sweep_params, outputs, local_output_collection
        )

        # if we are debugging
        if self.writer.config["debugging_data_dir"] is not None:
            local_filtered_values = np.zeros(
                (local_n_successful, len(local_filtered_dict["sweep_params"])),
                dtype=np.float64,
            )
            for i, (key, item) in enumerate(
                local_filtered_dict["sweep_params"].items()
            ):
                local_filtered_values[:, i] = item["value"][:]
        else:
            local_filtered_values = None

        # Not that we have all of the successful outputs in a consolidated dictionary locally,
        # we can now construct a global dictionary of successful solves.
        (
            global_filtered_dict,
            global_filtered_results,
            global_filtered_values,
        ) = self._aggregate_filtered_results(local_filtered_dict, req_num_samples)

        # Now we can save this
        self.comm.Barrier()

        # Save to file
        global_save_data = self.writer.save_results(
            sweep_params,
            local_filtered_values,
            global_filtered_values,
            local_filtered_dict,
            global_filtered_dict,
            global_filtered_results,
        )

        return global_save_data
