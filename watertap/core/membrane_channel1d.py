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

from pyomo.common.config import Bool, ConfigValue, In
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    NegativeReals,
    Param,
    Set,
    Var,
    value,
    units as pyunits,
)
from idaes.core import (
    declare_process_block_class,
    DistributedVars,
    FlowDirection,
    useDefault,
)
from idaes.core.base.control_volume1d import ControlVolume1DBlockData
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import (
    MembraneChannelMixin,
    PressureChangeType,
    CONFIG_Template as Base_CONFIG_Template,
)

CONFIG_Template = Base_CONFIG_Template()

CONFIG_Template.declare(
    "area_definition",
    ConfigValue(
        default=DistributedVars.uniform,
        domain=In(DistributedVars),
        description="Argument for defining form of area variable",
        doc="""Argument defining whether area variable should be spatially
variant or not. **default** - DistributedVars.uniform.
**Valid values:** {
DistributedVars.uniform - area does not vary across spatial domain,
DistributedVars.variant - area can vary over the domain and is indexed
by time and space.}""",
    ),
)

CONFIG_Template.declare(
    "transformation_method",
    ConfigValue(
        default=useDefault,
        description="Discretization method to use for DAE transformation",
        doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations.""",
    ),
)

CONFIG_Template.declare(
    "transformation_scheme",
    ConfigValue(
        default=useDefault,
        description="Discretization scheme to use for DAE transformation",
        doc="""Discretization scheme to use when transforming domain. See
Pyomo documentation for supported schemes.""",
    ),
)

CONFIG_Template.declare(
    "finite_elements",
    ConfigValue(
        default=10,
        domain=int,
        description="Number of finite elements in length domain",
        doc="""Number of finite elements to use when discretizing length 
        domain (default=10)""",
    ),
)

CONFIG_Template.declare(
    "collocation_points",
    ConfigValue(
        default=5,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
        discretizing length domain (default=5)""",
    ),
)


@declare_process_block_class("MembraneChannel1DBlock")
class MembraneChannel1DBlockData(MembraneChannelMixin, ControlVolume1DBlockData):
    def _skip_element(self, x):
        if self.config.transformation_scheme != "FORWARD":
            return x == self.length_domain.first()
        else:
            return x == self.length_domain.last()

    def apply_transformation(self, *args, **kwargs):
        super().apply_transformation(*args, **kwargs)
        self.difference_elements = Set(
            ordered=True,
            initialize=(x for x in self.length_domain if not self._skip_element(x)),
        )
        self._set_nfe()

    def add_geometry(
        self, length_var, width_var, flow_direction=FlowDirection.forward, **kwargs
    ):
        """
        Method to create spatial domain and volume Var in ControlVolume.

        Args:
            length_var - An external variable to use for the length of
                         the channel. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.

            width_var - An external variable to use for the width of
                         the channel. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.
            flow_direction - argument indicating direction of material flow
                            relative to length domain. Valid values:
                                - FlowDirection.forward (default), flow goes
                                  from 0 to 1.
                                - FlowDirection.backward, flow goes from 1 to 0
            length_domain - (optional) a ContinuousSet to use as the length
                            domain for the ControlVolume. If not provided, a
                            new ContinuousSet will be created (default=None).
                            ContinuousSet should be normalized to run between
                            0 and 1.
            length_domain_set - (optional) list of point to use to initialize
                            a new ContinuousSet if length_domain is not
                            provided (default = [0.0, 1.0]).

        Returns:
            None
        """
        super().add_geometry(
            length_var=length_var, flow_direction=flow_direction, **kwargs
        )
        self._add_var_reference(width_var, "width", "width_var")

    def add_state_blocks(self, has_phase_equilibrium=None):
        """
        This method constructs the state blocks for the
        control volume.

        Args:
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
        Returns:
            None
        """
        super().add_state_blocks(has_phase_equilibrium=has_phase_equilibrium)
        self._add_interface_stateblock(has_phase_equilibrium)

    def add_total_enthalpy_balances(self, **kwrags):
        # make this a no-op for MC1D
        return None

    def add_isothermal_conditions(self, **kwargs):

        super().add_isothermal_conditions(**kwargs)

        ## ==========================================================================
        # Feed-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isothermal assumption for feed channel",
        )
        def eq_feed_isothermal(b, t, x):
            if self._skip_element(x):
                return Constraint.Skip
            return (
                b.properties[t, b.length_domain.first()].temperature
                == b.properties[t, x].temperature
            )

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        add_object_reference(self, "dP_dx", self.deltaP)

    def initialize(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
        hold_state=True,
        initialize_guess=None,
    ):
        """
        Initialization routine for the membrane channel control volume

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output log level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None)
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - True. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.
            initialize_guess : a dict of guesses for solvent_recovery, solute_recovery,
                     and cp_modulus. These guesses offset the initial values
                     for the retentate, permeate, and membrane interface
                     state blocks from the inlet feed
                     (default =
                     {'deltaP': -1e4,
                     'solvent_recovery': 0.5,
                     'solute_recovery': 0.01,
                     'cp_modulus': 1.1})

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """

        # Get inlet state if not provided
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="control_volume")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="control_volume")

        state_args = self._get_state_args(initialize_guess, state_args)

        # intialize self.properties
        source_flags = super().initialize(
            state_args=state_args["feed_side"],
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        self.properties_interface.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["interface"],
        )

        if hold_state:
            return source_flags
        else:
            self.release_state(source_flags, outlvl)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables

        # will not override if the user provides the scaling factor
        ## default of 1 set by ControlVolume1D
        if iscale.get_scaling_factor(self.area) == 1:
            iscale.set_scaling_factor(self.area, 100)

        if hasattr(self, "pressure_change_total"):
            for v in self.pressure_change_total.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, "dP_dx"):
            for v in self.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e-5)
