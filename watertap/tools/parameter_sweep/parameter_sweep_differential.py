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
import numpy as np
from pyomo.common.config import ConfigValue

from watertap.tools.parameter_sweep.sampling_types import NormalSample
from watertap.tools.parameter_sweep.parameter_sweep import (
    _ParameterSweepBase,
    ParameterSweep,
)
from watertap.tools.MPI.dummy_mpi import DummyCOMM
from watertap.tools.parallel.parallel_manager import SingleProcessParallelManager


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

    CONFIG.declare(
        "differential_sweep_specs",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Dictionary containing the specifications for the differential sweep",
            doc="""
            A specification dictionary that contains details for how to construct the parameter sweep dictionary for differential sweep.
            This is a nested dictionary where the first level denotes the variable names for which the differential sweep needs to be carried out.
            The second level denotes various options to be used for wach variable.
            The number of samples for each differential sweep is specified while initializing the DifferentialParameterSweep object wsing the keyword `num_diff_samples`
            e.g.
            
            {
                "fs.a": {
                    "diff_mode": "sum",
                    "diff_sample_type": NormalSample,
                    "std_dev": 0.01,
                    "pyomo_object": m.fs.input["a"],
                },
                "fs.b": {
                    "diff_mode": "product",
                    "diff_sample_type": UniformSample,
                    "relative_lb": 0.01,
                    "relative_ub": 0.01,
                    "pyomo_object": m.fs.input["b"],
                },
                "fs.c": {
                    "diff_mode": "sum",
                    "diff_sample_type": GeomSample,
                    "relative_lb": 0.01,
                    "relative_ub": 10.0,
                    "pyomo_object": m.fs.input["c"],
                },
            }

            """,
        ),
    )

    def __init__(
        self,
        **options,
    ):

        # Initialize the base Class
        super().__init__(**options)

        if self.config.guarantee_solves:
            raise NotImplementedError

    def _create_differential_sweep_params(self, local_values):

        differential_sweep_specs = self.config.differential_sweep_specs

        diff_sweep_param = {}
        for ctr, (param, specs) in enumerate(differential_sweep_specs.items()):
            nominal_val = local_values[self.diff_spec_index[ctr]]
            pyomo_object = specs["pyomo_object"]
            if specs["diff_sample_type"] == NormalSample:
                std_dev = specs["std_dev"]
                diff_sweep_param[param] = NormalSample(
                    pyomo_object, nominal_val, std_dev, self.config.num_diff_samples
                )
            else:
                relative_lb = specs["relative_lb"]
                relative_ub = specs["relative_ub"]
                if specs["diff_mode"] == "sum":
                    lb = nominal_val * (1 - relative_lb)
                    ub = nominal_val * (1 + relative_ub)
                elif specs["diff_mode"] == "product":
                    lb = nominal_val * relative_lb
                    ub = nominal_val * relative_ub
                elif specs["diff_mode"] == "percentile":
                    lower_nominal = specs["nominal_lb"]
                    upper_nominal = specs["nominal_ub"]
                    delta_nominal = abs(upper_nominal - lower_nominal)
                    lb = nominal_val + delta_nominal * relative_lb
                    ub = nominal_val + delta_nominal * relative_ub
                else:
                    raise NotImplementedError
                diff_sweep_param[param] = specs["diff_sample_type"](
                    pyomo_object, lb, ub, self.config.num_diff_samples
                )

        return diff_sweep_param

    def _check_differential_sweep_key_validity(self, sweep_params):
        diff_specs_keys = list(self.config.differential_sweep_specs.keys())
        sweep_param_keys = list(sweep_params.keys())

        if all(key in sweep_param_keys for key in diff_specs_keys):
            self.diff_spec_index = [
                sweep_param_keys.index(key) for key in diff_specs_keys
            ]
        else:
            raise ValueError(
                "differential_sweep_specs keys don't match with sweep_param keys"
            )

    def _define_differential_sweep_outputs(self, sweep_params):
        self.differential_outputs = self.outputs
        if self.outputs is not None:
            for key in sweep_params.keys():
                if key not in self.config.differential_sweep_specs.keys():
                    self.differential_outputs[key] = sweep_params[key].pyomo_object

    def _create_local_output_skeleton(self, model, sweep_params, outputs, num_samples):
        output_dict = super()._create_local_output_skeleton(
            model, sweep_params, outputs, num_samples
        )
        output_dict["nominal_idx"] = np.arange(
            num_samples, dtype=float
        )  # [*range(num_samples)]
        output_dict["differential_idx"] = np.array([np.nan] * num_samples)
        return output_dict

    def _append_differential_results(self, local_output_dict, diff_results_dict):

        for idx, diff_sol in diff_results_dict.items():
            for key, item in diff_sol.items():
                # Solve status
                if key == "solve_successful":
                    n_diff_samples = len(item)
                    local_output_dict["solve_successful"].extend(item)
                    local_output_dict["nominal_idx"] = np.concatenate(
                        (
                            local_output_dict["nominal_idx"],
                            np.array([np.nan] * n_diff_samples),
                        ),
                        axis=0,
                    )
                    local_output_dict["differential_idx"] = np.concatenate(
                        (
                            local_output_dict["differential_idx"],
                            np.array([idx] * n_diff_samples, dtype=float),
                        ),
                        axis=0,
                    )

                else:
                    for subkey, subitem in item.items():
                        local_output_dict[key][subkey]["value"] = np.concatenate(
                            (
                                local_output_dict[key][subkey]["value"],
                                subitem["value"],
                            )
                        )

                    # We also need to capture sweep_params variables that are not a part of differential_sweep_specs
                    if key == "sweep_params":
                        missing_sub_keys = set(item.keys()) ^ set(
                            local_output_dict[key].keys()
                        )
                        # This loop shouldn't run if the above set is empty
                        for subkey in missing_sub_keys:
                            # We are picking the unchanged sweep_params from the outputs. In the ideal world, they would be the same.
                            local_output_dict["sweep_params"][subkey][
                                "value"
                            ] = np.concatenate(
                                (
                                    local_output_dict["sweep_params"][subkey]["value"],
                                    diff_sol["outputs"][subkey]["value"],
                                )
                            )

    def _collect_local_inputs(self, local_results_dict):

        num_local_samples = len(local_results_dict["solve_successful"])
        local_inputs = np.zeros(
            (num_local_samples, len(local_results_dict["sweep_params"])),
            dtype=float,
        )

        for i, (key, item) in enumerate(local_results_dict["sweep_params"].items()):
            local_inputs[:, i] = item["value"]

        return local_inputs

    def _aggregate_input_arr(self, global_results_dict, num_global_samples):

        global_values = np.zeros(
            (num_global_samples, len(global_results_dict["sweep_params"])),
            dtype=float,
        )

        if self.rank == 0:
            for i, (key, item) in enumerate(
                global_results_dict["sweep_params"].items()
            ):
                global_values[:, i] = item["value"]

        self.comm.Bcast(global_values, root=0)

        return global_values

    def _aggregate_results(self, local_output_dict):

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

    def _create_global_output(self, local_output_dict):  # , req_num_samples=None):
        global_output_dict = super()._create_global_output(local_output_dict)

        # We now need to get the mapping array. This only needs to happen on root
        local_num_cases_all = len(local_output_dict["solve_successful"])
        # AllGather the total size of the value array on each MPI rank
        sample_split_arr = self.comm.allgather(local_num_cases_all)
        num_total_samples = sum(sample_split_arr)
        # AllGather nominal values for creating the parallel offset
        nominal_sample_split_arr = self.comm.allgather(self.n_nominal_local)

        # We need to create a global index and offset items accordingly. This
        # needs to happen on all ranks/workers.
        my_rank = self.parallel_manager.get_rank()
        offset = 0
        if my_rank > 0:
            offset = sum(nominal_sample_split_arr[:my_rank])
        local_output_dict["nominal_idx"] = local_output_dict["nominal_idx"] + offset
        local_output_dict["differential_idx"] = (
            local_output_dict["differential_idx"] + offset
        )

        # Resize global index array
        if self.parallel_manager.is_root_process():
            global_output_dict["nominal_idx"] = np.zeros(num_total_samples, dtype=float)
            global_output_dict["differential_idx"] = np.zeros(
                num_total_samples, dtype=float
            )

        # Now we need to collect it on global_output_dict
        self.comm.Gatherv(
            sendbuf=local_output_dict["nominal_idx"],
            recvbuf=(
                global_output_dict["nominal_idx"],
                sample_split_arr,
            ),
            root=0,
        )
        self.comm.Gatherv(
            sendbuf=local_output_dict["differential_idx"],
            recvbuf=(
                global_output_dict["differential_idx"],
                sample_split_arr,
            ),
            root=0,
        )

        return global_output_dict

    def _run_differential_sweep(self, model, local_value):

        diff_sweep_param_dict = self._create_differential_sweep_params(local_value)

        # We want this instance of the parameter sweep to run in serial
        diff_ps = ParameterSweep(
            optimize_function=self.config.optimize_function,
            optimize_kwargs=self.config.optimize_kwargs,
            reinitialize_function=self.config.reinitialize_function,
            reinitialize_kwargs=self.config.reinitialize_kwargs,
            reinitialize_before_sweep=self.config.reinitialize_before_sweep,
            parallel_manager_class=SingleProcessParallelManager,
            comm=DummyCOMM,
        )

        _, differential_sweep_output_dict = diff_ps.parameter_sweep(
            model,
            diff_sweep_param_dict,
            combined_outputs=self.differential_outputs,
            num_samples=self.config.num_diff_samples,
            seed=self.seed,
        )

        return differential_sweep_output_dict

    def _run_sample(
        self,
        model,
        reinitialize_values,
        local_value_k,
        k,
        sweep_params,
        local_output_dict,
    ):

        run_successful = super()._run_sample(
            model,
            reinitialize_values,
            local_value_k,
            k,
            sweep_params,
            local_output_dict,
        )
        self.differential_sweep_output_dict[k] = self._run_differential_sweep(
            model, local_value_k
        )

        return run_successful

    def _do_param_sweep(self, model, sweep_params, outputs, local_values):
        self.differential_sweep_output_dict = {}

        local_output_dict = super()._do_param_sweep(
            model, sweep_params, outputs, local_values
        )

        # Now append the outputs of the differential solves
        self._append_differential_results(
            local_output_dict, self.differential_sweep_output_dict
        )

        return local_output_dict

    def parameter_sweep(
        self,
        model,
        sweep_params,
        outputs=None,
        num_samples=None,
        seed=None,
    ):

        # Create a base sweep_params
        sweep_params, sampling_type = self._process_sweep_params(sweep_params)

        # Check if the keys in the differential sweep specs exist in sweep params
        self._check_differential_sweep_key_validity(sweep_params)

        # Define differential sweep outputs
        self.outputs = outputs
        self._define_differential_sweep_outputs(sweep_params)

        # Set the seed before sampling
        self.seed = seed
        np.random.seed(self.seed)

        # Enumerate/Sample the parameter space
        global_values = self._build_combinations(
            sweep_params, sampling_type, num_samples
        )

        # divide the workload between processors
        local_values = self._divide_combinations(global_values)
        self.n_nominal_local = np.shape(local_values)[0]

        # Check if the outputs have the name attribute. If not, assign one.
        if outputs is not None:
            self._assign_variable_names(model, outputs)

        # Create a dictionary to store all the differential ps_objects
        self.diff_ps_dict = {}

        # Do the Loop
        if self.config.custom_do_param_sweep is None:
            local_results_dict = self._do_param_sweep(
                model,
                sweep_params,
                outputs,
                local_values,
            )
        else:
            local_results_dict = self.config.custom_do_param_sweep(
                model,
                sweep_params,
                outputs,
                local_values,
                **self.config.custom_do_param_sweep_kwargs,
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
