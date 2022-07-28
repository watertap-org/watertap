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
from watertap.tools.parameter_sweep_class import (_ParameterSweepBase, ParameterSweep)

class DifferentialParameterSweep(_ParameterSweepBase):

    def __init__(self,
        csv_results_file_name=None,
        h5_results_file_name=None,
        debugging_data_dir = None,
        interpolate_nan_outputs = False,
        guarantee_solves=False,
        *args,
        **kwargs,
    ):

    # Initialize the base Class
    _ParameterSweepBase.__init__(self)

    if guarantee_solves:
        raise NotImplementedError

    self.writer = ParameterSweepWriter(
        self.comm,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        debugging_data_dir=debugging_data_dir,
        interpolate_nan_outputs=interpolate_nan_outputs
    )

    def _process_sweep_params(self, sweep_params):
        pass

    def create_differential_sweep_params(self, local_values, differential_sweep_info):
        local_values = np.array([4.0e-12, 3.5e-8, 0.95])


        """ Kinshuk
        Differential sampling type -> how do we modify
        Attributes:
            Diff_mode = Add to the base value, multiply to the base value
            differential_sweep_info['value']['sample']= NormalSample(0.99,1.01,3) #product
            differential_sweep_info['value']['sample']= NormalSample(-1e-12,1e-12,3) #sum
        """

        """ Alex
        if Diff_mode=='sum':
            local_dif_values[i]=local_values[i]+differential_sweep_info['value']['sample']
        if Diff_mode=='product':
            local_dif_values[i]=local_values[i]*differential_sweep_info['value']['sample']

        """

        """ Kinshuk
        if Diff_mode=='sum':
            local_dif_values[i]=local_values[i]+differential_sweep_info['value']['sample']

        """
        differential_sweep_info = {
            "n_differential_samples" = same across everything
            "A_comp" : { "diff_mode" : sum/product ,
                         "diff_sample_type": Anything of the Sampling types,
                         "relative_lb" : For linear/geometric sampling,
                         "relative_ub" : For linear/geometric sampling,
                         "relative_std_dev" : For normal sampling
                        },
            "B_comp" : {},
            "Spacer_porosity"  : {},
        }


        return diff_sweep_param_dict

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
        probe_function,
    ):

        # Initialize space to hold results
        local_num_cases = np.shape(local_values)[0]

        # Create the output skeleton for storing detailed data
        local_output_dict = _create_local_output_skeleton(
            model, sweep_params, outputs, local_num_cases
        )

        local_results = np.zeros((local_num_cases, len(local_output_dict["outputs"])))

        local_solve_successful_list = []

        if reinitialize_function is not None:
            reinitialize_values = ComponentMap()
            for v in model.component_data_objects(pyo.Var):
                reinitialize_values[v] = v.value
        else:
            reinitialize_values = None

        # ================================================================
        # Run all optimization cases
        # ================================================================
        counter = 0
        for k in range(local_num_cases):

            # Step 1 : Run baseline case
            # Update the model values with a single combination from the parameter space
            _update_model_values(model, sweep_params, local_values[k, :])

            if probe_function is None or probe_function(model):
                run_successful = self._param_sweep_kernel(
                    model,
                    optimize_function,
                    optimize_kwargs,
                    reinitialize_before_sweep,
                    reinitialize_function,
                    reinitialize_kwargs,
                    reinitialize_values,
                )
            else:
                run_successful = False

            # Step 2: Run differential case
            self.diff_ps_dict[counter] = ParamweterSweep()
            self.create_differential_sweep_params(local_values[k,:])

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


    def parameter_sweep(
        self,
        model,
        sweep_params,
        outputs=None,
        optimize_function=None, # self._default_optimize,
        optimize_kwargs=None,
        reinitialize_function=None,
        reinitialize_kwargs=None,
        reinitialize_before_sweep=False,
        num_samples=None,
        seed=None,
    ):

        # Create a base sweep_params
        sweep_params_base, sampling_type = self._process_sweep_params(sweep_params)

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


        # Create a dictionary to store all the differential ps_objects
        self.diff_ps_dict = {}

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



if __name__ == '__main__':

    sweep_params = {}
    sweep_params["A_comp"] = LatinHypercubeSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
    sweep_params["B_comp"] = LatinHypercubeSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)
    sweep_params["Spacer_porosity"] = LatinHypercubeSample(
        m.fs.RO.spacer_porosity, 0.95, 0.99
    )

    # # TODO: We need to define what this dictionary looks like
    # differntial_params = {}
    # differntial_params["A_comp"] = LatinHypercubeSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
    # differntial_params["B_comp"] = LatinHypercubeSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)
    # differntial_params["Spacer_porosity"] = LatinHypercubeSample(
    #     m.fs.RO.spacer_porosity, 0.95, 0.99
    # )



    sweep_params_new = {}
    sweep_params_new['base'] = sweep_params
    sweep_params_new['differential'] = differntial_params
