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

import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock

# TODO:
# Create surrogate model during initialization if experimental dataset is provided

@declare_process_block_class("FlushingSurrogate")
class FlushingSurrogateData(UnitModelBlockData):
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


    # def create_surrogate_model(self):

    #     # Create the surrogate model for the flushing process
    #     # Check is a dataset file was passed

    #     if self.config.dataset_filename is not None:
    #         # Load dataset
    #         data = pd.read_csv(self.config.dataset_filename)
            
    #         pass



    def build(self):

        super().build()

        # Parameters
        self.number_tanks_in_series = Param(
            initialize=2,
            units=pyunits.dimensionless,
            doc="Number of tanks in series to represent the a plug flow reactor with mixing",
        )

        self.mean_residence_time = Var(
            initialize=30.0,
            units=pyunits.s,
            doc="Mean residence time in the system",
        )

        # Variables
        self.raw_feed_concentration = Var(
            initialize=0.0,
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the raw feed used for flushing",
        )

        self.pre_flushing_concentration = Var(
            initialize=0.0,
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the accumulation volume prior to flushing at the end of the concentration cycle",
        )

        self.post_flushing_concentration = Var(
            initialize=0.0,
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of the accumulation volume after flushing at the start of the concentration cycle",
        )

        self.flushing_efficiency = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Flushing efficiency of the system",
        )

        self.flushing_time = Var(
            initialize=20.0, units=pyunits.s, doc="Duration of flushing"
        )

        if self.config.dataset_filename is not None:
            self.create_surrogate_model()
        else:
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
        self.flushing_concentration_constraint = Constraint(
            expr=self.post_flushing_concentration
            == (1 - self.flushing_efficiency) * self.pre_flushing_concentration
            + self.flushing_efficiency * self.raw_feed_concentration
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

        pre_con_fixed = self.pre_flushing_concentration.fixed
        post_con_fixed = self.post_flushing_concentration.fixed
        self.pre_flushing_concentration.fix()
        self.post_flushing_concentration.unfix()

        # Check if surrogate model exists
        if hasattr(self, "surrogate"):
            # Initialize surrogate model
            self.init_data = pd.DataFrame(
                {
                    "time": [value(self.flushing_time)],
                    "mean_residence_time": [value(self.mean_residence_time)],
                }
            )
            self.init_output = self.surrogate.evaluate_surrogate(self.init_data)

            # Set initial values for model variables
            self.flushing_efficiency.set_value(self.init_output["F_t"].values[0])

        else:
            init_log.warning(
                "No surrogate model found. Using default flushing efficiency constraint."
            )
            self.flushing_time.fix()
            self.mean_residence_time.fix()


        # Solve unit
        opt = get_solver(solver, optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        if pre_con_fixed == False:
            self.pre_flushing_concentration.unfix()
        if post_con_fixed == False:
            self.post_flushing_concentration.unfix()

        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")
        
        self.flushing_time.unfix()
        self.mean_residence_time.unfix()

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        pass
