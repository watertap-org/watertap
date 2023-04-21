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
"""

"""

from enum import Enum
from pandas import DataFrame

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
)
from idaes.models.unit_models.separator import SeparatorData


from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from pyomo.environ import (
    Constraint,
    Param,
    Block,
    value,
    units as pyunits,
    check_optimal_termination,
    Set,
)

from idaes.core.util.exceptions import (
    ConfigurationError,
    PropertyNotSupportedError,
    InitializationError,
)

__author__ = "Alejandro Garciadiego"


# Set up logger
_log = idaeslog.getLogger(__name__)

# Enumerate options for balances
class SplittingType(Enum):
    """
    Enum of supported material split types.
    """

    totalFlow = 1
    phaseFlow = 2
    componentFlow = 3
    phaseComponentFlow = 4


class EnergySplittingType(Enum):
    """
    Enum of support energy split types.
    """

    equal_temperature = 1
    equal_molar_enthalpy = 2
    enthalpy_split = 3


@declare_process_block_class("Dewatering_Unit")
class DewateringUnit(SeparatorData):
    """
    Dewatering unit block for BSM2
    """

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(DewateringUnit, self).build()

        self.config.outlet_list = ["underflow", "overflow"]

        self.config.split_basis == SplittingType.componentFlow

        self.p_dewat = Param(
            initialize=0.28,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Percentage of suspended solids in the underflow",
        )

        self.TSS_rem = Param(
            initialize=0.98,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Percentage of suspended solids removed",
        )

        @self.Expression(self.flowsheet().time, doc="Suspended solid concentration")
        def TSS(blk, t):
            return 0.75 * (
                blk.inlet.conc_mass_comp[0, "X_I"]
                + blk.inlet.conc_mass_comp[0, "X_P"]
                + blk.inlet.conc_mass_comp[0, "X_BH"]
                + blk.inlet.conc_mass_comp[0, "X_BA"]
                + blk.inlet.conc_mass_comp[0, "X_S"]
            )

        @self.Expression(self.flowsheet().time, doc="Dewatering factor")
        def f_dewat(blk, t):
            return blk.p_dewat * (10 / (blk.TSS[t]))

        @self.Expression(self.flowsheet().time, doc="Remove factor")
        def f_q_du(blk, t):
            return blk.TSS_rem / 100 / blk.f_dewat[t]

        self.non_particulate_components = Set(
            initialize=[
                "S_I",
                "S_S",
                "S_O",
                "S_NO",
                "S_NH",
                "S_ND",
                "H2O",
                "S_ALK",
            ]
        )

        self.particulate_components = Set(
            initialize=["X_I", "X_S", "X_P", "X_BH", "X_BA", "X_ND"]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.particulate_components,
            doc="particulate fraction",
        )
        def overflow_particulate_fraction(blk, t, i):
            return blk.split_fraction[t, "overflow", i] == 1 - blk.TSS_rem

        self.display()

        @self.Constraint(
            self.flowsheet().time,
            self.non_particulate_components,
            doc="soluble fraction",
        )
        def non_particulate_components(blk, t, i):
            return blk.split_fraction[t, "overflow", i] == 1 - blk.f_q_du[t]

    def initialize_build(
        self, outlvl=idaeslog.NOTSET, optarg=None, solver=None, hold_state=False
    ):
        """
        Initialization routine for separator

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - False. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # Initialize mixed state block
        if self.config.mixed_state_block is not None:
            mblock = self.config.mixed_state_block
        else:
            mblock = self.mixed_state
        flags = mblock.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Solve for split fractions only
        component_status = {}
        for c in self.component_objects((Block, Constraint)):
            for i in c:
                if not c[i].local_name == "sum_split_frac":
                    # Record current status of components to restore later
                    component_status[c[i]] = c[i].active
                    c[i].deactivate()

        if degrees_of_freedom(self) != 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
                init_log.info(
                    "Initialization Step 1 Complete: {}".format(idaeslog.condition(res))
                )

        for c, s in component_status.items():
            if s:
                c.activate()

        if self.config.ideal_separation:
            # If using ideal splitting, initialization should be complete
            return flags

        # Initialize outlet StateBlocks
        outlet_list = self.create_outlet_list()

        # Premises for initializing outlet states:
        # 1. Intensive states remain unchanged - this is either a valid premise
        # or the actual state is impossible to calculate without solving the
        # full separator model.
        # 2. Extensive states are use split fractions if index matches, or
        # average of split fractions for outlet otherwise
        for o in outlet_list:
            # Get corresponding outlet StateBlock
            o_block = getattr(self, o + "_state")

            # Create dict to store fixed status of state variables
            o_flags = {}
            for t in self.flowsheet().time:

                # Calculate values for state variables
                s_vars = o_block[t].define_state_vars()
                for v in s_vars:
                    for k in s_vars[v]:
                        # Record whether variable was fixed or not
                        o_flags[t, v, k] = s_vars[v][k].fixed

                        # If fixed, use current value
                        # otherwise calculate guess from mixed state and fix
                        if not s_vars[v][k].fixed:
                            m_var = getattr(mblock[t], s_vars[v].local_name)
                            if "flow" in v:
                                # Need average split fraction
                                avg_split = value(
                                    sum(
                                        self.split_fraction[t, o, j]
                                        for (p, j) in mblock.phase_component_set
                                    )
                                    / len(mblock.phase_component_set)
                                )

            # Call initialization routine for outlet StateBlock
            o_block.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                hold_state=False,
            )

            # Revert fixed status of variables to what they were before
            for t in self.flowsheet().time:
                s_vars = o_block[t].define_state_vars()
                for v in s_vars:
                    for k in s_vars[v]:
                        s_vars[v][k].fixed = o_flags[t, v, k]

        if self.config.mixed_state_block is None:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)

            if not check_optimal_termination(res):
                raise InitializationError(
                    f"{self.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
                )

            init_log.info(
                "Initialization Step 2 Complete: {}".format(idaeslog.condition(res))
            )
        else:
            init_log.info("Initialization Complete.")

        if hold_state is True:
            return flags
        else:
            self.release_state(flags, outlvl=outlvl)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        """
        if self.config.mixed_state_block is None:
            mblock = self.mixed_state
        else:
            mblock = self.config.mixed_state_block

        mblock.release_state(flags, outlvl=outlvl)

    def calculate_scaling_factors(self):
        mb_type = self.config.material_balance_type
        mixed_state = self.get_mixed_state_block()
        t_ref = self.flowsheet().time.first()
        mb_type = mixed_state[t_ref].default_material_balance_type()
        super().calculate_scaling_factors()

        if hasattr(self, "temperature_equality_eqn"):
            for (t, i), c in self.temperature_equality_eqn.items():
                s = iscale.get_scaling_factor(
                    mixed_state[t].temperature, default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s)

        if hasattr(self, "pressure_equality_eqn"):
            for (t, i), c in self.pressure_equality_eqn.items():
                s = iscale.get_scaling_factor(
                    mixed_state[t].pressure, default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s)

        if hasattr(self, "material_splitting_eqn"):
            if mb_type == MaterialBalanceType.componentTotal:
                for (t, _, j), c in self.material_splitting_eqn.items():
                    for i, p in enumerate(mixed_state.phase_list):
                        if i == 0:
                            s = iscale.get_scaling_factor(ft, default=1)
                        else:
                            _s = iscale.get_scaling_factor(ft, default=1)
                            s = _s if _s < s else s
                    iscale.constraint_scaling_transform(c, s)

    def _get_performance_contents(self, time_point=0):
        if hasattr(self, "split_fraction"):
            var_dict = {}
            for k in self.split_fraction.keys():
                if k[0] == time_point:
                    var_dict[f"Split Fraction [{str(k[1:])}]"] = self.split_fraction[k]
            return {"vars": var_dict}
        else:
            return None

    def _get_stream_table_contents(self, time_point=0):
        outlet_list = self.create_outlet_list()

        if not self.config.ideal_separation:
            io_dict = {}
            if self.config.mixed_state_block is None:
                io_dict["Inlet"] = self.mixed_state
            else:
                io_dict["Inlet"] = self.config.mixed_state_block

            for o in outlet_list:
                io_dict[o] = getattr(self, o + "_state")

            return create_stream_table_dataframe(io_dict, time_point=time_point)

        else:
            stream_attributes = {}
            stream_attributes["Units"] = {}

            for n in ["inlet"] + outlet_list:
                port_obj = getattr(self, n)

                stream_attributes[n] = {}

                for k in port_obj.vars:
                    for i in port_obj.vars[k]:
                        if isinstance(i, float):
                            quant = report_quantity(port_obj.vars[k][time_point])
                            stream_attributes[n][k] = quant.m
                            stream_attributes["Units"][k] = quant.u
                        else:
                            if len(i) == 2:
                                kname = str(i[1])
                            else:
                                kname = str(i[1:])
                            quant = report_quantity(port_obj.vars[k][time_point, i[1:]])
                            stream_attributes[n][k + " " + kname] = quant.m
                            stream_attributes["Units"][k + " " + kname] = quant.u

            return DataFrame.from_dict(stream_attributes, orient="columns")
