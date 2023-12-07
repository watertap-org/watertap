import copy
import numpy as np


class _ParameterSweepParallelUtils:
    @classmethod
    def remove_unpicklable_state(cls, parameter_sweep_instance):
        """
        Remove and return any state from the ParameterSweep object that cannot be
        pickled, to make the instance picklable. Needed in order to use the
        ConcurrentFuturesParallelManager.
        """
        saved_state = {
            "parallel_manager": parameter_sweep_instance.parallel_manager,
            "writer": parameter_sweep_instance.writer,
        }

        parameter_sweep_instance.parallel_manager = None
        parameter_sweep_instance.writer = None
        return saved_state

    @classmethod
    def restore_unpicklable_state(cls, parameter_sweep_instance, state):
        """
        Restore a collection of saved state that was removed in order to pickle
        the ParameterSweep object.
        """
        parameter_sweep_instance.parallel_manager = state.get("parallel_manager", None)
        parameter_sweep_instance.writer = state.get("writer", None)

    """
    Combine all of the results retrieved from calling gather().
    - all_results is a list of Result objects, each representing the 
    parameters and results from running the optimization on one process. 
    """

    def _combine_gather_results(self, all_results):
        if len(all_results) == 0:
            return None

        # create the output skeleton based on the first set of results
        # we assume the results are in dict format
        initial_results = all_results[0].results

        combined_results = copy.deepcopy(initial_results)

        # remove any lingering pyomo objects, and convert inner results to numpy arrays
        for key, val in combined_results.items():
            if key != "solve_successful":
                if isinstance(val, dict):
                    for subkey, subval in val.items():
                        if "_pyo_obj" in subval:
                            del subval["_pyo_obj"]

        # for each result, concat the "value" array of results into the
        # gathered results to combine them all

        # get length of data in first result for finding missing keys
        total_chunk_length = len(all_results[0].results["solve_successful"])

        for i, result in enumerate(all_results[1:]):
            results = result.results

            for key, val in results.items():
                if key == "solve_successful":
                    combined_results[key] = np.append(
                        combined_results[key], copy.deepcopy(val)
                    )
                    continue
                if key == "nominal_idx" or key == "differential_idx":
                    combined_results[key] = np.append(
                        combined_results[key], copy.deepcopy(val)
                    )
                    continue
                # print("vall all results!", key, val)
                for subkey, subval in val.items():
                    # lets catch any keys that don' exist in result[0] and
                    # create empty array with expected length, after which we will add
                    # additional values, or add nan's instead
                    if subkey not in combined_results[key]:
                        # create empty array, as none of results so far had this key\

                        combined_results[key][subkey] = {}
                        for sub_subkey, value in subval.items():
                            if sub_subkey == "value":
                                combined_results[key][subkey]["value"] = (
                                    np.zeros(total_chunk_length) * np.nan
                                )
                            else:
                                combined_results[key][subkey][sub_subkey] = value
                    combined_results[key][subkey]["value"] = np.append(
                        combined_results[key][subkey]["value"],
                        copy.deepcopy(
                            subval["value"],
                        ),
                    )
                    # keep track of our subchunk_length
                    sub_chunk_length = len(subval["value"])

                # make sure we add any empty value to missing keys

                for subkey in combined_results[key]:
                    if subkey not in val.keys():
                        empty_chunk = np.zeros(sub_chunk_length) * np.nan
                        combined_results[key][subkey]["value"] = np.append(
                            combined_results[key][subkey]["value"], empty_chunk
                        )
            total_chunk_length += sub_chunk_length
        return combined_results

    """
    Build up a list of the outputs for each result of the optimization.
    Returned as a list of lists, where each inner list is the results from
    one process's run.
    """

    def _combine_output_array(self, gathered_results):
        outputs = gathered_results["outputs"]
        if len(outputs) == 0:
            return []

        # assume all output arrays have the same length
        combined_outputs = [
            np.asarray([]) for _ in range(len(list(outputs.values())[0]["value"]))
        ]
        for _, output in outputs.items():
            for i in range(len(output["value"])):
                combined_outputs[i] = np.append(combined_outputs[i], output["value"][i])
        return np.asarray(combined_outputs)

    """
    Build up a list of the sweep_inputs for each result of the optimization.
    Returned as a list of lists, where each inner list is the results from
    one process's run.
    """

    def _combine_input_array(self, gathered_results):
        inputs = gathered_results["sweep_params"]
        if len(inputs) == 0:
            return []

        # assume all output arrays have the same length
        combined_inputs = [
            np.asarray([]) for _ in range(len(list(inputs.values())[0]["value"]))
        ]
        for _, inputv in inputs.items():
            for i in range(len(inputv["value"])):
                combined_inputs[i] = np.append(combined_inputs[i], inputv["value"][i])
        return np.asarray(combined_inputs)

    """
    Use the embedded ParallelManager to fan out and then back in the results.
    Args:
    - build_model: a function for building the flowsheet model
    - build_model_kwargs: any keyword args necessary for the build_model function
    - build_sweep_params: a function for building the sweep parameters
    - build_sweep_params_kwargs: any keyword args necessary for the build_sweep_params
    function
    - build_outputs: a function for building the outputs dictionary
    - all_parameter_combinations: a list where each element represents the parameters
    for a single local run
    Returns:
    - a list of LocalResults representing the results of the simulation runs 
    """

    def run_scatter_gather(self, all_parameter_combinations, class_reference):
        # save a reference to the parallel manager since it will be removed
        # along with the other unpicklable state
        parallel_manager = self.parallel_manager
        saved_state = class_reference.remove_unpicklable_state(self)

        do_build_kwargs = {"param_sweep_instance": self}

        parallel_manager.scatter(
            do_build,
            do_build_kwargs,
            do_execute,
            all_parameter_combinations,
        )

        # gather the results and combine them into the format we want
        all_results = parallel_manager.gather()
        class_reference.restore_unpicklable_state(self, saved_state)

        return all_results


def do_build(
    param_sweep_instance,
):
    """
    Used to pass into the parallel manager to build the parameters necessary
    for the sweep function. Defined at the top level so it's picklable.
    """
    ps_config = param_sweep_instance.config
    model = ps_config.build_model(**ps_config.build_model_kwargs)
    sweep_params = ps_config.build_sweep_params(
        model, **ps_config.build_sweep_params_kwargs
    )
    sweep_params, sampling_type = param_sweep_instance._process_sweep_params(
        sweep_params
    )
    outputs = ps_config.build_outputs(model, **ps_config.build_outputs_kwargs)
    # for when differential parameter tool is used
    if hasattr(param_sweep_instance, "_define_differential_sweep_outputs"):
        param_sweep_instance._define_differential_sweep_outputs(model, sweep_params)

    if outputs is not None:
        param_sweep_instance.assign_variable_names(model, outputs)

    return [param_sweep_instance, model, sweep_params, outputs]


def do_execute(
    local_combo_array,
    param_sweep_instance,
    model,
    sweep_params,
    outputs,
):
    """
    Used to pass into the parallel manager in order to execute the sweep
    for a set of local values. Defined at the top level so it's picklable.
    """

    if param_sweep_instance.config.custom_do_param_sweep is not None:
        return param_sweep_instance.config.custom_do_param_sweep(
            param_sweep_instance, sweep_params, outputs, local_combo_array
        )

    return param_sweep_instance._do_param_sweep(
        sweep_params, outputs, local_combo_array
    )


def return_none(model, outputkeys=None):
    """
    Used so that build_outputs=None is a valid usage of the parameter sweep tool
    without requiring the user to wrap it in a function.
    """
    return None
