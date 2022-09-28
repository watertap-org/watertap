import numpy as np
import pyomo.environ as pyo
import sys
import os
import itertools
import warnings
import copy, pprint
import h5py
import pathlib

from scipy.interpolate import griddata
from enum import Enum, auto
from abc import abstractmethod, ABC
from idaes.core.solvers import get_solver

from idaes.surrogate.pysmo import sampling
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.tee import capture_output
from pyomo.common.config import ConfigValue

from watertap.tools.parameter_sweep.parameter_sweep_writer import ParameterSweepWriter
from watertap.tools.parameter_sweep.sampling_types import *  # (SamplingType, LinearSample)
from watertap.tools.parameter_sweep.parameter_sweep import (
    _ParameterSweepBase,
    ParameterSweep,
)


class DifferentialParameterSweep(_ParameterSweepBase):

    CONFIG = _ParameterSweepBase.CONFIG()

    CONFIG.declare(
        "num_diff_samples",
        ConfigValue(
            default=1,
            domain=int,
            description="Number of differntial sweep samples",
        ),
    )

    CONFIG.declare(
        "guarantee_solves",
        ConfigValue(
            default=False,
            domain=bool,
            description="Guarantee a pre-specified number of solves.",
        ),
    )

    def __init__(
        self,
        # csv_results_file_name=None,
        # h5_results_file_name=None,
        # debugging_data_dir = None,
        # interpolate_nan_outputs = False,
        # guarantee_solves=False,
        # num_diff_samples=1,
        **options,
    ):

        # Initialize the base Class
        super().__init__(**options)

        if self.config.guarantee_solves:
            raise NotImplementedError

        # self.num_diff_samples = num_diff_samples

    def _create_differential_sweep_params(self, local_values, differential_sweep_specs):

        # Assume that we have local values and example differential sweep_info
        # local_values = np.array([4.0e-12, 3.5e-8, 0.95])
        # differential_sweep_info = {
        #     # "n_differential_samples" = same across everything
        #     "A_comp" : { "diff_mode" : sum/product ,
        #                  "diff_sample_type": Anything of the Sampling types,
        #                  "relative_lb" : For linear/geometric sampling,
        #                  "relative_ub" : For linear/geometric sampling,
        #                  "relative_std_dev" : For normal sampling
        #                 },
        #     "B_comp" : {},
        #     "Spacer_porosity"  : {},
        # }

        diff_sweep_param = {}
        ctr = 0
        for param, specs in differential_sweep_specs.items():
            nominal_val = local_values[ctr]
            pyomo_object = specs["pyomo_object"]
            if specs["diff_sample_type"] == NormalSample:
                std_dev = specs["relative_std_dev"]
                diff_sweep_param[param] = NormalSample(
                    pyomo_object, nominal_val, std_dev
                )
            else:
                relative_lb = specs["relative_lb"]
                relative_ub = specs["relative_ub"]
                if specs["diff_mode"] == "sum":
                    lb = nominal_val * (1 - relative_lb)
                    ub = nominal_val * (1 + relative_lb)
                elif specs["diff_mode"] == "product":
                    lb = nominal_val * relative_lb
                    ub = nominal_val * relative_ub
                else:
                    raise NotImplementedError
                diff_sweep_param[param] = specs["diff_sample_type"](
                    pyomo_object, lb, ub
                )

        return diff_sweep_param

    def _append_differential_results(self, local_output_dict, diff_results_dict):

        for idx, diff_sol in diff_results_dict.items():
            for key, item in diff_sol.items():
                # Solve status
                if key == "solve_successful":
                    local_output_dict["solve_successful"].extend(item)
                else:
                    for subkey, subitem in item.items():
                        local_output_dict[key][subkey]["value"] = np.concatenate(
                            (local_output_dict[key][subkey]["value"], subitem["value"])
                        )

    def _create_global_output(self, local_output_dict):
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
                    ]["value"]

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
                    global_output_dict[key] = global_solve_successful

        return global_output_dict

    def _collect_local_inputs(self, local_results_dict):

        num_local_samples = len(local_results_dict["solve_successful"])
        local_inputs = np.zeros(
            (num_local_samples, len(local_results_dict["sweep_params"])),
            dtype=np.float64,
        )

        for i, (key, item) in enumerate(local_results_dict["sweep_params"].items()):
            local_inputs[:, i] = item["value"]

        return local_inputs

    def _aggregate_input_arr(self, global_results_dict, num_global_samples):

        global_values = np.zeros(
            (num_global_samples, len(global_results_dict["sweep_params"])),
            dtype=np.float64,
        )

        if self.rank == 0:
            for i, (key, item) in enumerate(
                global_results_dict["sweep_params"].items()
            ):
                global_values[:, i] = item["value"]

        if self.num_procs > 1:  # pragma: no cover
            self.comm.Bcast(global_values, root=0)

        return global_values

    def _aggregate_results(self, local_output_dict):

        num_local_samples = len(local_output_dict["solve_successful"])

        # Create the global results dictionary
        global_results_dict = self._create_global_output(local_output_dict)

        # Broadcast the number of global samples to all ranks
        num_global_samples = len(global_results_dict["solve_successful"])
        num_global_samples = self.comm.bcast(num_global_samples, root=0)

        global_results_arr = self._aggregate_results_arr(
            global_results_dict, num_global_samples
        )
        global_input_values = self._aggregate_input_arr(
            global_results_dict, num_global_samples
        )

        return (
            global_results_dict,
            global_results_arr,
            global_input_values,
            num_global_samples,
        )

    def _do_param_sweep(
        self,
        model,
        sweep_params,
        differential_sweep_specs,
        outputs,
        local_values,
        # optimize_function,
        # optimize_kwargs,
        # reinitialize_function,
        # reinitialize_kwargs,
        # reinitialize_before_sweep,
        # probe_function,
    ):

        # Initialize space to hold results
        local_num_cases = np.shape(local_values)[0]

        # Create the output skeleton for storing detailed data
        local_output_dict = self._create_local_output_skeleton(
            model, sweep_params, outputs, local_num_cases
        )

        # local_results = np.zeros((local_num_cases, len(local_output_dict["outputs"])))

        local_solve_successful_list = []

        if self.config.reinitialize_function is not None:
            reinitialize_values = ComponentMap()
            for v in model.component_data_objects(pyo.Var):
                reinitialize_values[v] = v.value
        else:
            reinitialize_values = None

        # ================================================================
        # Run all optimization cases
        # ================================================================
        counter = 0
        differential_sweep_output_dict = {}

        for k in range(local_num_cases):

            # Step 1 : Run baseline/nominal case
            # Update the model values with a single combination from the parameter space
            self._update_model_values(model, sweep_params, local_values[k, :])

            if self.config.probe_function is None or probe_function(model):
                run_successful = self._param_sweep_kernel(
                    model,
                    # optimize_function,
                    # optimize_kwargs,
                    # reinitialize_before_sweep,
                    # reinitialize_function,
                    # reinitialize_kwargs,
                    reinitialize_values,
                )
            else:
                run_successful = False

            # Update the loop based on the reinitialization for baseline values
            self._update_local_output_dict(
                model,
                sweep_params,
                k,
                local_values[k, :],
                run_successful,
                local_output_dict,
            )

            local_solve_successful_list.append(run_successful)

            # Step 2: Run differential case
            # self.diff_ps_dict[counter] = ParamweterSweep()
            diff_sweep_param_dict = self._create_differential_sweep_params(
                local_values[k, :], differential_sweep_specs
            )
            diff_ps = ParameterSweep(
                optimize_function=self.config.optimize_function,  # self._default_optimize,
                optimize_kwargs=self.config.optimize_kwargs,
                reinitialize_function=self.config.reinitialize_function,
                reinitialize_kwargs=self.config.reinitialize_kwargs,
                reinitialize_before_sweep=self.config.reinitialize_before_sweep,
            )

            _, differential_sweep_output_dict[k] = diff_ps.parameter_sweep(
                model,
                diff_sweep_param_dict,
                outputs=outputs,
                # optimize_function=optimize_function, # self._default_optimize,
                # optimize_kwargs=optimize_kwargs,
                # reinitialize_function=reinitialize_function,
                # reinitialize_kwargs=None,
                # reinitialize_before_sweep=reinitialize_before_sweep,
                num_samples=self.config.num_diff_samples,
                seed=self.seed,
            )

        local_output_dict["solve_successful"] = local_solve_successful_list

        # Now append the outputs of the differential solves
        self._append_differential_results(
            local_output_dict, differential_sweep_output_dict
        )

        return local_output_dict

    def parameter_sweep(
        self,
        model,
        sweep_params,
        differential_sweep_specs,
        outputs=None,
        # optimize_function=None,
        # optimize_kwargs=None,
        # reinitialize_function=None,
        # reinitialize_kwargs=None,
        # reinitialize_before_sweep=False,
        # probe_function=None,
        num_samples=None,
        seed=None,
    ):

        # Create a base sweep_params
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)

        # Set the seed before sampling
        self.seed = seed
        np.random.seed(self.seed)

        # Enumerate/Sample the parameter space
        global_values = self._build_combinations(
            sweep_params, sampling_type, num_samples
        )

        # divide the workload between processors
        local_values = self._divide_combinations(global_values)
        local_num_cases = np.shape(local_values)[0]

        # Create a dictionary to store all the differential ps_objects
        self.diff_ps_dict = {}

        # Do the Loop
        local_results_dict = self._do_param_sweep(
            model,
            sweep_params,
            differential_sweep_specs,
            outputs,
            local_values,
            # optimize_function,
            # optimize_kwargs,
            # reinitialize_function,
            # reinitialize_kwargs,
            # reinitialize_before_sweep,
            # probe_function,
        )

        # re-writing local_values
        local_values = self._collect_local_inputs(local_results_dict)

        # Aggregate results on Master
        (
            global_results_dict,
            global_results_arr,
            global_input_arr,
            num_global_samples,
        ) = self._aggregate_results(local_results_dict)

        # Save to file
        global_save_data = self.writer.save_results(
            sweep_params,
            local_values,
            global_input_arr,
            local_results_dict,
            global_results_dict,
            global_results_arr,
        )

        return global_results_dict, global_save_data