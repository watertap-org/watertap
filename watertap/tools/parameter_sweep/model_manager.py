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


class ModelManager:
    def __init__(self, ps_instance):
        """this gets ps_instance info so it can manage
        model build, init and solve logic
        """
        self.ps_conf = ps_instance.config
        self.ps = ps_instance
        self.is_intilized = False
        self.is_solved = False
        self.is_params_updated = False
        self.is_prior_parameter_solved = False

    def build_and_init(self, params=None, local_value_k=None):
        """build and init model, if required by user, update paramaters before init"""
        print(self.ps_conf.build_model)
        self.model = self.ps_conf.build_model(**self.ps_conf.build_model_kwargs)
        self.is_params_updated = False
        if self.ps_conf.initialize_function is not None:
            if (
                self.ps_conf.update_sweep_params_before_init
                and params != None
                and local_value_k != None
            ):
                self.update_model_params(sweep_params, local_value_k)
            self.ps_conf.initialize_function(
                self.model, **self.ps_conf.initialize_kwargs
            )
        # raise v
        elif self.ps_conf.initialize_before_sweep:
            raise ValueError(
                "Initialization function was not specified. The model will not be reinitialized."
            )
        self.is_intilized = True
        self.is_solved = False

    def add_initialized_model(self, model):
        self.model = model
        self.is_intilized = True
        self.is_prior_parameter_solved = True

    def init_model(self):
        try:
            self.ps_conf.initialize_function(
                self.model, **self.ps_conf.initialize_kwargs
            )
            self.is_intilized = True
        except TypeError:
            # this happens if the optimize_kwargs are misspecified,
            # which is an error we want to raise
            self.is_intilized = False
            self.is_solved = False
            raise
        except:
            self.is_intilized = False
            self.is_solved = False

    def update_model_params(self, sweep_params, local_value_k):
        self.ps._update_model_values(self.model, sweep_params, local_value_k)
        self.is_params_updated = True

    def solve_model(self):
        # update to determine if we are solving from initilized or pre-solved state
        self.is_prior_parameter_solved = self.is_solved
        try:
            results = self.ps_conf.optimize_function(
                self.model, **self.ps_conf.optimize_kwargs
            )
            pyo.assert_optimal_termination(results)
            self.is_solved = True
        except TypeError:
            # this happens if the optimize_kwargs are misspecified,
            # which is an error we want to raise
            self.is_intilized = False
            self.is_solved = False
            raise
        except:
            self.is_intilized = False
            self.is_solved = False

        return results
