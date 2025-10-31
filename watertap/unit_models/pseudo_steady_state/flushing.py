import sys
from io import StringIO
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
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock

# TODO:
# Surrogate model as a function of residence time and number of tanks in series
# Create surrogate model during initialization if not provided

@declare_process_block_class("FlushingSurrogate")
class FlushingSurrogateData(UnitModelBlockData):
    """
    This is a surrogate unit model for flushing.
    1. A surrogate model can be created as a function of experimental data passed to configuration variable 'dataset_filename'
    2. If no experimental data is provided, the default values for number_of_tanks_in_series will be used to create a profile for the given mean residence time.
    3. If a surrogate file is passed to the configuration variable "surrogate_model_file", the previously created surrogate file will be used.
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
        rtd_profile = pd.DataFrame()
        time = np.linspace(0, 3 * 35, 100)
        N= self.number_tanks_in_series

        for t_m in np.linspace(0, 2*35, 5):
            scale = t_m / N
            F_t = gamma.cdf(time, a=N, scale=scale)
            df = pd.DataFrame({"time": time, "F_t": F_t})
            df["mean_residence_time"] = t_m
            df["N"] = N
            rtd_profile = pd.concat([rtd_profile, df], ignore_index=True)

        self.rtd_profile = rtd_profile

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
            training_fraction = 0.8
            input_labels = ["time", "mean_residence_time"]
            output_labels = ["F_t"]

            # Create training and validation sets
            data_training, val_df = split_training_validation(
                self.rtd_profile, training_fraction
            )

            PATH = os.path.dirname(os.path.abspath(__file__))
            # If no filename to save surrogate model is provided, use default with number of tanks in series
            if self.config.surrogate_filename_save is None:
                self.config.surrogate_filename_save = os.path.join(
                    PATH, f"flushing_surrogate_n_{self.number_tanks_in_series()}_tau_{self.mean_residence_time()}.json"
                )
            

            # Create surrogate model
            self.surrogate = self.create_rbf_surrogate(
                data_training=data_training,
                input_labels=input_labels,
                output_labels=output_labels,
                output_filename=self.config.surrogate_filename_save,
                tee=True,
            )

            self.surrogate_blk = SurrogateBlock(concrete=True)

            self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.flushing_time,
            output_vars=self.flushing_efficiency,
        )



    def create_rbf_surrogate(
        self, data_training, input_labels, output_labels, output_filename=None, tee=False
    ):
        # Capture long output
        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        # Create PySMO trainer object
        trainer = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=data_training,
        )

        # Set PySMO options
        trainer.config.basis_function = "gaussian"  # default = gaussian
        trainer.config.solution_method = "algebraic"  # default = algebraic
        trainer.config.regularization = True  # default = True

        # Train surrogate
        rbf_train = trainer.train_surrogate()

        # Remove autogenerated 'solution.pickle' file
        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass
        except Exception as e:
            raise e

        # Create callable surrogate object
        xmin, xmax = 0,1000
        input_bounds = {
            input_labels[i]: (xmin, xmax) for i in range(len(input_labels))
        }
        rbf_surr = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds)

        # Save model to JSON
        if output_filename is not None:
            model = rbf_surr.save_to_file(output_filename, overwrite=True)

        # Revert back to standard output
        sys.stdout = oldstdout

        if tee:
            # Display first 50 lines and last 50 lines of output
            celloutput = stream.getvalue().split("\n")
            for line in celloutput[:50]:
                print(line)
            print(".")
            print(".")
            print(".")
            for line in celloutput[-50:]:
                print(line)

        return rbf_surr


    def load_surrogate(self):
        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.config.surrogate_model_file)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=[self.flushing_time, self.mean_residence_time],  
            output_vars=self.flushing_efficiency,
        )

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

        # Variables - Not sure if these are needed or if their created with the surrogate model
        self.flushing_efficiency = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Flushing efficiency of the system",
        )

        self.flushing_time = Var(
            initialize=20.0, units=pyunits.s, doc="Duration of flushing"
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

        # Initialize surrogate
        pre_con_fixed = self.pre_flushing_concentration.fixed
        post_con_fixed = self.post_flushing_concentration.fixed
        self.pre_flushing_concentration.fix()
        self.post_flushing_concentration.unfix()
        self.init_data = pd.DataFrame(
            {
                "time": [value(self.flushing_time)],
                "mean_residence_time": [value(self.mean_residence_time)],
            }
        )

        self.init_output = self.surrogate.evaluate_surrogate(self.init_data)

        # Set initial values for model variables
        self.flushing_efficiency.set_value(self.init_output["F_t"].values[0])

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

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        pass
