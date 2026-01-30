import sys
from io import StringIO
import pandas as pd
import numpy as np
from scipy.stats import gamma
import os
import math

from pyomo.environ import (
    Var,
    Param,
    Constraint,
    Expression,
    value,
    check_optimal_termination,
    units as pyunits,
    exp,
)

from idaes.core.util.exceptions import ConfigurationError
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)

from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from scipy.optimize import rosen, differential_evolution

import idaes.core.util.math as idaesMath

__author__ = "Mukta Hardikar, Kurban Sitterley, Alexander V. Dudchenko"


@declare_process_block_class("Flushing")
class FlushingData(UnitModelBlockData):
    """
    This is a surrogate unit model for flushing.
    1. If filepath is passed to configuration variable 'dataset_filename', the number of tanks in series and mean residence time are estimated using the experimental data provided in the file.
    2. If no experimental data is provided, the default values for number_of_tanks_in_series will be used to create a profile for the given mean residence time.

    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Solar energy surrogate models are steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Solar energy models do not include holdup""",
        ),
    )

    CONFIG.declare(
        "surrogate_model_file",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to surrogate model file",
            doc="""User provided surrogate model .json file.""",
        ),
    )
    CONFIG.declare(
        "surrogate_filename_save",
        ConfigValue(
            default=None,
            domain=str,
            description="Filename used to save surrogate model to .json",
            doc="""Filename used to save surrogate model file to .json""",
        ),
    )
    CONFIG.declare(
        "dataset_filename",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to data file",
            doc="""Path to data file. Must be a .csv""",
        ),
    )
    CONFIG.declare(
        "number_tanks_in_series",
        ConfigValue(
            default=2,
            domain=PositiveInt,
            description="Number of tanks in series to represent the a plug flow reactor with mixing",
            doc="Number of tanks in series to represent the a plug flow reactor with mixing",
        ),
    )

    def get_data(self, file_path):
        """Load data from a CSV file and return as a pandas DataFrame."""
        data = pd.read_csv(file_path)
        t = data["t"].to_numpy()
        f_t = data["F_t"].to_numpy()

        return t, f_t

    def ft_profile(self, params, t):
        """Calculate the profile function F(t) based on the given parameters."""
        t_m, N = params
        tau = t_m / int(N)
        N = int(N)
        sum_terms = [np.zeros(len(t))]
        for n in range(N - 1):
            term = ((t / tau) ** n) / math.factorial(n)
            sum_terms.append(term)
        # print(np.sum(sum_terms, axis=0))
        if N == 1:
            F_t = 1 - np.exp(-t / tau)
        else:
            F_t = 1 - np.exp(-t / tau) * np.sum(sum_terms, axis=0)
        # print(F_t)
        return F_t

    def objective_function(self, params, t, f_t):
        """Calculate the sum of squared errors between the model and data."""
        f_t_model = self.ft_profile(params, t)
        return np.sum((f_t - f_t_model) ** 2)

    def fit_rtd_profile(self):
        """Fit the RTD profile using the provided dataset filename in the configuration.

        The function uses the profile_fitting module to optimize the number of tanks in series
        and mean residence time based on the experimental data provided in the CSV file.
        """
        data_file = self.config.dataset_filename
        if not os.path.isfile(data_file):
            raise ConfigurationError(
                f"The specified dataset file '{data_file}' does not exist."
            )

        t_data, f_t_data = self.get_data(data_file)

        bounds = [(1e-5, 500), (1, 50)]  # Bounds for tau and N
        result = differential_evolution(
            self.objective_function,
            bounds,
            # workers=8,
            args=(t_data, f_t_data),
            tol=0.0001,
        )

        tau_opt, N_opt = result.x
        N_opt = int(N_opt)

        self.mean_residence_time.fix(tau_opt * N_opt)
        self.number_tanks_in_series.set_value(N_opt)

        idaeslog.getLogger(self.name).info(
            f"Fitted RTD profile: mean residence time = {tau_opt * N_opt}, number of tanks in series = {N_opt}"
        )

    def build(self):

        super().build()
        if self.config.dataset_filename is not None:
            self.fit_rtd_profile()
        # Parameters
        self.number_tanks_in_series = Param(
            initialize=self.config.number_tanks_in_series,
            units=pyunits.dimensionless,
            mutable=False,
            doc="Number of tanks in series to represent the a plug flow reactor with mixing",
        )

        self.mean_residence_time = Var(
            initialize=30.0,
            units=pyunits.s,
            doc="Mean residence time in the system",
        )

        # Variables
        self.flushing_feed_concentration = Var(
            initialize=0.0,
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the raw feed used for flushing",
        )

        self.pre_flushing_concentration = Var(
            initialize=0.0,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the accumulation volume prior to flushing at the end of the concentration cycle",
        )

        self.post_flushing_concentration = Var(
            initialize=0.0,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the accumulation volume after flushing at the start of the concentration cycle",
        )

        self.flushing_efficiency = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Flushing efficiency of the system",
        )
        self.minimum_flushing_efficiency = Var(
            initialize=0.05,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Flushing efficiency of the system",
        )
        self.minimum_flushing_efficiency.fix()
        self.flushing_time = Var(
            initialize=20.0,
            bounds=(1, None),
            units=pyunits.s,
            doc="Duration of flushing",
        )

        # Add the constraint to calculate flushing efficiency using parameters in the unit model
        @self.Constraint(
            doc="Constraint to calculate flushing efficiency using number of tanks in series and mean residence time"
        )
        def flushing_efficiency_constraint(b):
            N = b.number_tanks_in_series
            t = b.flushing_time
            tau = b.mean_residence_time / N  # mean residence time PER TANK
            theta = t / tau
            series_sum = sum(theta**m / math.factorial(m) for m in range(N()))
            return b.flushing_efficiency == 1.0 - exp(-theta) * series_sum

        # Constraint to calculate the concentration of the accumulation volume at the end of flushing
        self.minimum_flushing_efficiency_constraint = Constraint(
            expr=self.minimum_flushing_efficiency <= self.flushing_efficiency
        )
        self.flushing_concentration_constraint = Constraint(
            # expr=idaesMath.smooth_min(  # self.post_flushing_concentration
            #     self.post_flushing_concentration,
            #     self.pre_flushing_concentration
            #     * (1 - self.minimum_flushing_efficiency),
            # )
            expr=self.post_flushing_concentration
            == (1 - self.flushing_efficiency) * self.pre_flushing_concentration
            + self.flushing_efficiency * self.flushing_feed_concentration
        )

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Check if pre-flushing and post-flushing concentrations are fixed before intialization
        assert (
            degrees_of_freedom(self) == 0
        ), "Flushing unit model has non-zero degrees of freedom at the time of initialization. DOFs are {}".format(
            degrees_of_freedom(self)
        )

        # Solve unit
        opt = get_solver(solver, optarg)
        # with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        res = opt.solve(self, tee=False)
        self.display()

        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        self.display()

    def calculate_scaling_factors(self):

        iscale.set_scaling_factor(self.flushing_feed_concentration, 1)
        iscale.set_scaling_factor(self.mean_residence_time, 1)
        iscale.set_scaling_factor(self.flushing_efficiency, 1)

        iscale.set_scaling_factor(self.minimum_flushing_efficiency, 1)
        iscale.constraint_scaling_transform(self.flushing_efficiency_constraint, 10)
        iscale.constraint_scaling_transform(self.flushing_concentration_constraint, 1)

        iscale.constraint_scaling_transform(
            self.minimum_flushing_efficiency_constraint, 1
        )
