import pandas as pd
import numpy as np
from scipy.stats import gamma
import os

from pyomo.environ import (
    Var,
    Param,
    Constraint,
    Expression,
    value,
    check_optimal_termination,
    units as pyunits,
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
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock


@declare_process_block_class("FlushingSurrogate")
class FlushingSurrogateData(UnitModelBlockData):
    """
    This is a surrogate unit model for flushing.
    1. A surrogate model can be created as a function of experimental data passed to configuration variable 'dataset_filename'
    2. If not data is passed the default values for number_of_tanks_in_series will be used to create a profile for the given flow_rate and accumulation volume.
    3. If a surrogate file is passed to the configuration variable "surrogate_model_file" it is used.
    4. If no surrogate file is provided, a surrogate model is created during initialization using the experimental data provided in "dataset_filename" or default values number_of_tanks_in_series.
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


    def generate_rtd_profile(self):
        # Generate the RTD profile for the flushing process

        mean_residence_time = self.accumulation_volume / self.flushing_flow_rate  # Mean residence time in s
        time = np.linspace(0, 3 * mean_residence_time, 1000)
        scale = mean_residence_time / self.number_tanks_in_series
        F_t = gamma.cdf(time, a=self.number_tanks_in_series, scale=scale)
        data = pd.DataFrame({"time": time, "F_t": F_t})

        self.rtd_profile = data


    def create_surrogate_model(self):

        # Create the surrogate model for the flushing process
        # Check is a dataset file was passed

        if self.config.dataset_filename is not None:
            # TODO: Create surrogate using experimental data
            pass

        else:
            # Create a default surrogate model with generated rtd profile
            self.generate_rtd_profile()

            # Create a surrogate using the default rtd profile
            # TODO: Surrogate creation and fit check - Update the 



    def load_surrogate(self):
        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.config.surrogate_model_file)
        self.surrogate_blk.build_model(
            self.surrogate, 
            input_vars=self.flushing_time, 
            output_vars=self.flushing_efficiency
            )
        

    def build(self):
        
        super().build()

        # Parameters
        self.number_tanks_in_series = Param(
            initialize=1, 
            units=pyunits.dimensionless,
            doc = "Number of tanks in series to represent the a plug flow reactor with mixing"
            )

        self.accumulation_volume = Param(
            initialize=0.0, 
            units=pyunits.m**3,
            doc="Accumulation volume through which is being flushing"
            )
        
        self.flushing_flow_rate = Param(
            initialize=0.0, 
            units=pyunits.m**3/pyunits.s,
            doc = "Flow rate of the flushing water"
            )

        # Variables
        self.raw_feed_concentration = Var(
            initialize=0.0, 
            units=pyunits.kg/pyunits.m**3,
            doc = "Concentration of the raw feed used for flushing"
            )
        
        self.pre_flushing_concentration = Var(
            initialize=0.0, 
            units=pyunits.kg/pyunits.m**3,
            doc = "Concentration of the accumulation volume prior to flushing"
            )
        
        self.post_flushing_concentration = Var(
            initialize=0.0, 
            units=pyunits.kg/pyunits.m**3,
            doc = "Concentration of the accumulation volume after flushing"
            )


        # Variables - Not sure if these are needed or if their created with the surrogate model
        self.flushing_efficiency = Var(
            initialize=0.1, 
            units=pyunits.dimensionless,
            bounds = (0,1),
            doc="Flushing efficiency of the system"
            )
           
        self.flushing_time = Var(
            initialize=20.0, 
            units=pyunits.s,
            doc = "Duration of flushing"
            )


        # If a surrogate model is passed, it is used to calculate the concentration 
        if self.config.surrogate_model_file is not None:
            try:
                self.load_surrogate()
            except Exception as e:
                err_msg = f"Error loading surrogate model: {e}"
                raise ConfigurationError(err_msg)
        else:
            self.create_surrogate_model()


        # Constraint to calculate the concentration of the accumulation volume at the end of flushing
        self.flushing_concentration_constraint = Constraint(
                expr= self.post_flushing_concentration
                      == (1-self.flushing_efficiency) * self.pre_flushing_concentration + self.flushing_efficiency*self.raw_feed_concentration
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

        # Initialize surrogate
        self.init_data = pd.DataFrame(
            {
                "time" : [value(self.flushing_time)],
            }
        )

        self.init_output = self.surrogate.evaluate_surrogate(self.init_data)

        # Set initial values for model variables
        self.flushing_efficiency.set_value(self.init_output["F_t"].values[0])

        # Solve unit
        opt = get_solver(solver, optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        pass


        


