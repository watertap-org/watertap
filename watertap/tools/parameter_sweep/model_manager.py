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
import pyomo.environ as pyo


__author__ = "Alexander V. Dudchenko (SLAC)"


class ModelManager:
    def __init__(self, ps_instance):
        """This class manages the model state
        It uses the ps_instance info so it can manage
        model build, init and solve logic, as well as get any other config paramters
        """
        self.ps_conf = ps_instance.config
        self.ps = ps_instance

        self.is_initialized = False
        self.is_solved = False
        self.is_prior_parameter_solved = False
        # this is isused for loggin states if enabled
        self.solved_states = {"state": [], "local_value_k": []}
        self.initialized_states = {"state": [], "local_value_k": []}
        self.current_k = None

    def build_and_init(self, sweep_params=None, local_value_k=None):
        """build and init model, if required by user, update paramaters before init"""
        self.model = self.ps_conf.build_model(**self.ps_conf.build_model_kwargs)
        # intilized model if init function is passed in
        if self.ps_conf.initialize_function is not None:
            # update paramters before init if enabled by user
            if (
                self.ps_conf.update_sweep_params_before_init
                and sweep_params != None
                and local_value_k != None
            ):
                self.update_model_params(sweep_params, local_value_k)
            # init
            self.init_model()
        # raise error if user sets to init before sweep, but does not provide
        # initilize function
        elif self.ps_conf.update_sweep_params_before_init:
            raise ValueError(
                "Initialization function was not specified. The model will not be reinitialized with specified paramters."
            )
        self.update_initialized_state(True)
        self.update_solved_state(False)

    def update_initialized_state(self, state):
        """used to update and if requested log intilized state"""
        if self.ps_conf.log_model_states:
            self.initialized_states["state"].append(state)
            self.initialized_states["local_value_k"].append(self.current_k)

        self.is_initialized = state

    def update_solved_state(self, state):
        """used to update and if requested log solved state"""
        if self.ps_conf.log_model_states:
            self.solved_states["state"].append(state)
            self.solved_states["local_value_k"].append(self.current_k)
        self.is_solved = state

    def init_model(self):
        """attempt to init model"""
        try:
            self.ps_conf.initialize_function(
                self.model, **self.ps_conf.initialize_kwargs
            )
            self.update_initialized_state(True)
        except TypeError:
            # this happens if the optimize_kwargs are misspecified,
            # which is an error we want to raise
            self.update_solved_state(False)
            self.update_initialized_state(False)
            raise
        except:
            self.update_solved_state(False)
            self.update_initialized_state(False)

    def update_model_params(self, sweep_params, local_value_k):
        """Update local params, and store current local k state"""
        # print(self.model, sweep_params, local_value_k)
        self.ps._update_model_values(self.model, sweep_params, local_value_k)
        self.current_k = local_value_k

    def solve_model(self):
        """Attempt to solve the unit model"""
        # update to determine if we are solving from initilized or pre-solved state
        self.is_prior_parameter_solved = self.is_solved
        try:
            results = self.ps_conf.optimize_function(
                self.model, **self.ps_conf.optimize_kwargs
            )
            pyo.assert_optimal_termination(results)
            self.update_solved_state(True)
            self.update_initialized_state(True)
            return results
        except TypeError:
            # this happens if the optimize_kwargs are misspecified,
            # which is an error we want to raise
            self.update_solved_state(False)
            self.update_initialized_state(False)
            raise
        except:
            self.update_solved_state(False)
            self.update_initialized_state(False)
        return None
