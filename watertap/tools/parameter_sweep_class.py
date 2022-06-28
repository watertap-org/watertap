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
from pyomo.common.collections import ComponentSet
from pyomo.common.tee import capture_output

from watertap.tools.parameter_sweep_writer import ParameterSweepWriter
from watertap.tools.sampling_types import (SamplingType, LinearSample)

np.set_printoptions(linewidth=200)


class _ParameterSweepBase(ABC):
    def __init__(self, *args, **kwargs):
        self.comm = self._init_mpi()
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

    # ================================================================

    def _init_mpi(self):
        try:
            from mpi4py import MPI
            return MPI.COMM_WORLD
        except:
            dummy_comm = DummyMPI()
            return dummy_comm

    # ================================================================

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

    # ================================================================

    def _divide_combinations(self, global_combo_array):

        # Split the total list of combinations into NUM_PROCS chunks,
        # one per each of the MPI ranks
        # divided_combo_array = np.array_split(global_combo_array, num_procs, axis=0)
        divided_combo_array = np.array_split(global_combo_array, self.num_procs)

        # Return only this rank's portion of the total workload
        local_combo_array = divided_combo_array[self.rank]

        return local_combo_array

    # ================================================================

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

    # ================================================================

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

    # ================================================================

    def _default_optimize(self, model, options=None, tee=False):
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

    # ================================================================

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

    # ================================================================

    def _create_local_output_skeleton(self, model, sweep_params, outputs, num_samples):

        output_dict = {}
        output_dict["sweep_params"] = {}
        output_dict["outputs"] = {}

        sweep_param_objs = ComponentSet()

        # Store the inputs
        for sweep_param in sweep_params.values():
            var = sweep_param.pyomo_object
            sweep_param_objs.add(var)
            output_dict["sweep_params"][var.name] = self._create_component_output_skeleton(
                var, num_samples
            )

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
                output_dict["outputs"][short_name] = self._create_component_output_skeleton(
                    pyo_obj, num_samples
                )

        return output_dict


    # ================================================================


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

    # ================================================================

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

    # ================================================================

    def _create_global_output(self, local_output_dict, req_num_samples):

        # Before we can create the global dictionary, we need to delete the pyomo
        # object contained within the dictionary
        for key, val in local_output_dict.items():
            if key != "solve_successful":
                for subval in val.values():
                    if "_pyo_obj" in subval:
                        del subval["_pyo_obj"]

        if self.num_procs == 1:
            global_output_dict = local_output_dict
        else:  # pragma: no cover
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
            for key, item in local_output_dict.items():  # This probably doesnt work
                if key != "solve_successful":
                    for subkey, subitem in item.items():
                        comm.Gatherv(
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

    # ================================================================

    def _do_param_sweep(
        self,
        model,
        sweep_params,
        outputs,
        local_values,
        optimize_function,
        optimize_kwargs,
        reinitialize_function,
        reinitialize_kwargs,
        reinitialize_before_sweep,
    ):

        # Initialize space to hold results
        local_num_cases = np.shape(local_values)[0]

        # Create the output skeleton for storing detailed data
        local_output_dict = self._create_local_output_skeleton(
            model, sweep_params, outputs, local_num_cases
        )

        local_results = np.zeros((local_num_cases, len(local_output_dict["outputs"])))

        local_solve_successful_list = []

        # ================================================================
        # Run all optimization cases
        # ================================================================

        for k in range(local_num_cases):
            # Update the model values with a single combination from the parameter space
            self._update_model_values(model, sweep_params, local_values[k, :])

            run_successful = False  # until proven otherwise

            # Forced reinitialization of the flowsheet if enabled
            if reinitialize_before_sweep:
                try:
                    assert reinitialize_function is not None
                except:
                    raise ValueError(
                        "Reinitialization function was not specified. The model will not be reinitialized."
                    )
                else:
                    reinitialize_function(model, **reinitialize_kwargs)

            try:
                # Simulate/optimize with this set of parameter
                with capture_output():
                    results = optimize_function(model, **optimize_kwargs)
                pyo.assert_optimal_termination(results)

            except:
                # run_successful remains false. We try to reinitialize and solve again
                if reinitialize_function is not None:
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

    # ================================================================

    @abstractmethod
    def parameter_sweep(self, *args, **kwargs):
        raise NotImplementedError

    # ================================================================

class ParameterSweep(_ParameterSweepBase):

    def __init__(self,
        csv_results_file_name=None,
        h5_results_file_name=None,
        debugging_data_dir = None,
        interpolate_nan_outputs = False,
        *args,
        **kwargs
    ):

        # Initialize the base Class
        _ParameterSweepBase.__init__(self)

        # Initialize the writer
        if self.rank == 0:
            if h5_results_file_name is None and csv_results_file_name is None:
                warnings.warn(
                    "No results will be writen to disk as h5_results_file_name and csv_results_file_name are both None"
                )
                self.write_outputs = False
            else:
                self.write_outputs = True
                self.writer = ParameterSweepWriter(
                    self.comm,
                    csv_results_file_name=csv_results_file_name,
                    h5_results_file_name=h5_results_file_name,
                    debugging_data_dir=debugging_data_dir,
                    interpolate_nan_outputs=interpolate_nan_outputs
                )



    def _aggregate_local_results(
        self,
        global_values,
        local_output_dict,
        num_samples,
        local_num_cases
    ):

        # Create the dictionary
        global_results_dict = self._create_global_output(local_output_dict, num_samples)

        # Create the array
        num_global_samples = np.shape(global_values)[0]
        global_results_arr = self._aggregate_results_arr(global_results_dict, num_global_samples)

        return global_results_dict, global_results_arr

    def parameter_sweep(
        self,
        model,
        sweep_params,
        outputs=None,
        csv_results_file_name=None,
        h5_results_file_name=None,
        optimize_function=None, # self._default_optimize,
        optimize_kwargs=None,
        reinitialize_function=None,
        reinitialize_kwargs=None,
        reinitialize_before_sweep=False,
        debugging_data_dir=None,
        interpolate_nan_outputs=False,
        num_samples=None,
        seed=None,
    ):

        """
        This function offers a general way to perform repeated optimizations
        of a model for the purposes of exploring a parameter space while
        monitoring multiple outputs.
        If provided, writes single CSV file to ``results_file`` with all inputs and resulting outputs.

        Arguments:

            model : A Pyomo ConcreteModel containing a watertap flowsheet, for best
                    results it should be initialized before being passed to this
                    function.

            sweep_params: A dictionary containing the values to vary with the format
                          ``sweep_params['Short/Pretty-print Name'] =
                          (model.fs.variable_or_param[index], lower_limit, upper_limit, num_samples)``.
                          A uniform number of samples ``num_samples`` will be take between
                          the ``lower_limit`` and ``upper_limit``.

            outputs : An optional dictionary containing "short names" as keys and and Pyomo objects
                      on ``model`` whose values to report as values. E.g.,
                      ``outputs['Short/Pretty-print Name'] = model.fs.variable_or_expression_to_report``.
                      If not provided, i.e., outputs = None, the default behavior is to save all model
                      variables, parameters, and expressions which provides very thorough results
                      at the cost of large file sizes.

            csv_results_file_name (optional) : The path and file name to write a csv file. The default `None`
                                               does not write a csv file.

            h5_results_file_name (optional) : The path and file name to write a h5 file. The default `None`
                                              does not write a file.
                                              Writing an h5 file will also create a companion text file `{h5_results_file_name}.txt`
                                              which contains the variable names contained within the H5 file.

            optimize_function (optional) : A user-defined function to perform the optimization of flowsheet
                                           ``model`` and loads the results back into ``model``. The first
                                           argument of this function is ``model``. The default uses the
                                           default IDAES solver, raising an exception if the termination
                                           condition is not optimal.

            optimize_kwargs (optional) : Dictionary of kwargs to pass into every call to
                                         ``optimize_function``. The first arg will always be ``model``,
                                         e.g., ``optimize_function(model, **optimize_kwargs)``. The default
                                         uses no kwargs.

            reinitialize_function (optional) : A user-defined function to perform the re-initialize the
                                               flowsheet ``model`` if the first call to ``optimize_function``
                                               fails for any reason. After ``reinitialize_function``, the
                                               parameter sweep tool will immediately call
                                               ``optimize_function`` again.

            reinitialize_kwargs (optional) : Dictionary or kwargs to pass into every call to
                                             ``reinitialize_function``. The first arg will always be
                                             ``model``, e.g.,
                                             ``reinitialize_function(model, **reinitialize_kwargs)``.
                                             The default uses no kwargs.

            reinitialize_before_sweep (optional): Boolean option to reinitialize the flow sheet model before
                                                  every parameter sweep realization. The default is False.
                                                  Note the parameter sweep model will try to reinitialize the
                                                  solve regardless of the option if the run fails.

            mpi_comm (optional) : User-provided MPI communicator for parallel parameter sweeps.
                                  If None COMM_WORLD will be used. The default is sufficient for most
                                  users.

            debugging_data_dir (optional) : Save results on a per-process basis for parallel debugging
                                            purposes. If None no `debugging` data will be saved.

            interpolate_nan_outputs (optional) : When the parameter sweep has finished, interior values
                                                 of np.nan will be replaced with a value obtained via
                                                 a linear interpolation of their surrounding valid neighbors.
                                                 If true, a second output file with the extension "_clean"
                                                 will be saved alongside the raw (un-interpolated) values.

            num_samples (optional) : If the user is using sampling techniques rather than a linear grid
                                     of values, they need to set the number of samples

            seed (optional) : If the user is using a random sampling technique, this sets the seed

        Returns:

            save_data : A list were the first N columns are the values of the parameters passed
                        by ``sweep_params`` and the remaining columns are the values of the
                        simulation identified by the ``outputs`` argument.
        """

        # # Get an MPI communicator
        # comm, rank, num_procs = self._init_mpi(mpi_comm)

        # if self.rank == 0:
        #     if h5_results_file_name is None and csv_results_file_name is None:
        #         warnings.warn(
        #             "No results will be writen to disk as h5_results_file_name and csv_results_file_name are both None"
        #         )

        # Convert sweep_params to LinearSamples
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)

        # Set the seed before sampling
        np.random.seed(seed)

        # Enumerate/Sample the parameter space
        global_values = self._build_combinations(sweep_params, sampling_type, num_samples)

        # divide the workload between processors
        local_values = self._divide_combinations(global_values)
        local_num_cases = np.shape(local_values)[0]

        # Set up optimize_kwargs
        if optimize_kwargs is None:
            optimize_kwargs = dict()
        # Set up reinitialize_kwargs
        if reinitialize_kwargs is None:
            reinitialize_kwargs = dict()

        # Do the Loop
        local_results_dict = self._do_param_sweep(
            model,
            sweep_params,
            outputs,
            local_values,
            optimize_function,
            optimize_kwargs,
            reinitialize_function,
            reinitialize_kwargs,
            reinitialize_before_sweep,
        )

        # Aggregate results on Master
        global_results_dict, global_results_arr = self._aggregate_local_results(
            global_values,
            local_results_dict,
            num_samples,
            local_num_cases
        )

        # Save to file
        global_save_data = self.writer.save_results(sweep_params, local_values, global_values, local_results_dict,
            global_results_dict, global_results_arr)

        return global_save_data
