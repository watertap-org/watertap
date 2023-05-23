import pyomo.environ as pyo
import copy
import numpy as np
import os
import idaes.logger as idaeslog
import platform
import time
from pyomo.common.collections import ComponentSet, ComponentMap


from idaes.core.util import to_json, from_json

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

from idaes.core.util import to_json, from_json

from watertap.tools.parameter_sweep.async_param_sweep.async_utils import (
    _create_component_output_skeleton,
    _create_local_output_skeleton,
    _merge_copied_dicts,
    _convert_sweep_params_to_dict,
    _convert_outputs_to_dict,
    _update_model_values_from_dict,
)


def do_async_sweep(self, model, sweep_params, outputs, local_values):
    # check if analysis tools are available in watertap
    try:
        from watertap.tools.analysis_tools.model_state_tool import modelStateStorage
        import watertap.tools.analysis_tools.step_optimize_tool as stepTool

        if self.config.custom_do_param_sweep_kwargs.get("use_analysis_tools") == None:
            self.config.custom_do_param_sweep_kwargs["use_analysis_tools"] = True

    except ImportError:
        self.config.custom_do_param_sweep_kwargs["use_analysis_tools"] = False

    # serialize model
    if self.config.custom_do_param_sweep_kwargs.get("load_form_json"):
        if self.config.custom_do_param_sweep_kwargs["use_analysis_tools"]:
            model_state = modelStateStorage(model)
            json_string = model_state.get_dict_state()
        else:
            json_string = to_json(model, return_json_string=True)
    else:
        json_string = None
    # create ouptut dict for paramActor
    outputs = _convert_outputs_to_dict(outputs)
    # create paramActor option dict (Must be pickle safe)
    parall_kwargs = {}

    # standard paramsweep info
    parall_kwargs["optimize_function"] = self.config.optimize_function
    parall_kwargs["optimize_kwargs"] = self.config.optimize_kwargs
    parall_kwargs["reinitialize_before_sweep"] = self.config.reinitialize_before_sweep
    parall_kwargs["reinitialize_function"] = self.config.reinitialize_function
    parall_kwargs["reinitialize_kwargs"] = self.config.reinitialize_kwargs
    parall_kwargs["probe_function"] = self.config.probe_function
    # custom param sweep info, required!
    parall_kwargs["build_function"] = self.config.custom_do_param_sweep_kwargs.get(
        "build_function"
    )
    parall_kwargs["build_kwargs"] = self.config.custom_do_param_sweep_kwargs.get(
        "build_kwargs"
    )
    parall_kwargs["sweep_kwargs"] = self.config.custom_do_param_sweep_kwargs.get(
        "build_kwargs"
    )

    parall_kwargs["html_notice"] = None  # not implemented
    parall_kwargs["use_analysis_tools"] = self.config.custom_do_param_sweep_kwargs[
        "use_analysis_tools"
    ]
    if self.config.custom_do_param_sweep_kwargs.get("step_tool_options") != None:
        parall_kwargs[
            "step_tool_options"
        ] = self.config.custom_do_param_sweep_kwargs.get("step_tool_options")
    parall_kwargs["json_string"] = json_string
    parall_kwargs["outputs"] = outputs
    # check if ray is avaiallbe for paraell solving

    # create loal run dict
    run_dict = _convert_sweep_params_to_dict(sweep_params, local_values)
    # if self.config.custom_do_param_sweep_kwargs.get("num_workers")!= None:
    num_workers = self.config.custom_do_param_sweep_kwargs.get("num_workers")

    try:
        import ray

        ray_unavailable = True

    except ImportError:
        ray_unavailable = False

    if (
        self.config.custom_do_param_sweep_kwargs.get("use_mp") == True
        or ray_unavailable
    ):
        from watertap.tools.parameter_sweep.async_param_sweep import (
            multiprocessing_sweep,
        )

        result = multiprocessing_sweep.do_mp_sweep(parall_kwargs, run_dict, num_workers)
    else:
        from watertap.tools.parameter_sweep.async_param_sweep import rayio_sweep

        result = rayio_sweep.do_rayio_sweep(parall_kwargs, run_dict, num_workers)
    return _merge_copied_dicts(result)


class paramActor:
    def __init__(self, options):
        self.build_function = options["build_function"]
        self.build_kwargs = options["build_kwargs"]
        self.reinitialize_before_sweep = options["reinitialize_before_sweep"]
        self.reinitialize_function = options["reinitialize_function"]
        self.reinitialize_kwargs = options["reinitialize_kwargs"]
        self.optimize_function = options["optimize_function"]
        self.optimize_kwargs = options["optimize_kwargs"]
        self.outputs = options["outputs"]
        self.probe_function = options["probe_function"]
        self.json_string = options["json_string"]
        self.sweep_kwargs = options["sweep_kwargs"]
        self.use_analysis_tools = options["use_analysis_tools"]
        if options.get("step_tool_options") != None:
            self.try_final_first = options.get("try_final_first")
            self.num_steps = options.get("num_steps")
            self.re_steps = options.get("re_steps")
        else:
            self.try_final_first = True
            self.num_steps = 5
            self.re_steps = 2
        self.model_init = False
        self.model = None
        self.build_model()

    def build_model(self):
        del self.model
        self.model = self.build_function(**self.build_kwargs)
        if self.json_string is not None:
            if self.use_analysis_tools:
                self.model.state_after_box_solve = modelStateStorage(
                    self.model, restore_dict=self.json_string
                )
            else:
                from_json(self.model, sd=None, fname=None, s=self.json_string)

            self.model_init = True

    def init_model(self):
        self.build_model()
        if self.reinitialize_function is not None:
            self.reinitialize_function(self.model, **self.reinitialize_kwargs)
            self.model_init = True

    def get_probe_result(self, m, sweep_params):
        if self.probe_function is None:
            return True
        else:
            try:
                probe_result = self.probe_function(m, sweep_params)
            except:
                print("----------------------------------")
                print("error in probe function!!!!!!!!!!!")
                probe_result = True
            return probe_result

    def step_optimize(
        self,
        m,
        solver=None,
    ):
        solver = get_solver()
        results = stepTool.step_optimize(
            m,
            solver,
            self.optimize_function,
            self.model.state_after_box_solve,
            steps=self.num_steps,
            try_final_first=self.try_final_first,
            re_steps=self.re_steps,
        )
        return results

    def opt(self):
        solver = get_solver()
        if self.use_analysis_tools:
            return self.step_optimize(self.model)
        else:
            return self.optimize_function(self.model, solver)

    def solve_problem(self, sweep_params):
        order_index = sweep_params["order_index"]
        sweep_params.pop("order_index", None)

        run_successful = False

        _update_model_values_from_dict(self.model, sweep_params)
        if self.get_probe_result(self.model, self.reinitialize_kwargs):
            for i in ["Try #0", "Try #1"]:
                if self.reinitialize_function is not None and self.model_init == False:
                    try:
                        self.init_model()
                        self.model_init = True
                    except:
                        self.model_init = False
                        run_successful = False
                        init_okay = False

                _update_model_values_from_dict(self.model, sweep_params)

                if self.model_init and self.get_probe_result(
                    self.model, self.reinitialize_kwargs
                ):
                    try:
                        results = self.step_optimize(
                            self.model
                        )  # , **self.optimize_kwargs)
                        pyo.assert_optimal_termination(results)
                        run_successful = True
                        break
                    except:
                        run_successful = False
                        self.model_init = False
                        if (
                            self.reinitialize_before_sweep == False
                            and self.reinitialize_function is None
                        ):
                            break

                else:
                    print("Model failed to init, skipping condition!")
                    break
        output_dict = _create_local_output_skeleton(
            self.model, sweep_params, self.outputs, run_successful
        )
        output_dict["order_index"] = order_index
        print("RUN Complete, solution optimal: {}".format(run_successful))
        print(
            "Param {} \n Build option {} \n Init option {}".format(
                sweep_params, self.build_kwargs, self.reinitialize_kwargs
            )
        )
        return output_dict
