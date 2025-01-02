#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import (
    NegativeReals,
    Set,
    Var,
)
from idaes.core import (
    declare_process_block_class,
    FlowDirection,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
from idaes.core.base.control_volume1d import ControlVolume1DBlockData
import idaes.logger as idaeslog
from .MD_channel_base import (
    MDChannelMixin,
    PressureChangeType,
)

__author__ = "Elmira Shamlou"


@declare_process_block_class("MDChannel1DBlock")
class MDChannel1DBlockData(MDChannelMixin, ControlVolume1DBlockData):
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
        add_object_reference(self, "width", width_var)

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        add_object_reference(self, "dP_dx", self.deltaP)

        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.deltaP_channel = Var(
            self.flowsheet().config.time,
            initialize=-1e5,
            bounds=(-1e6, 0),
            domain=NegativeReals,
            units=units_meta("pressure"),
            doc="total prossure drop across the channel",
        )

    def _add_deltaP(self, pressure_change_type=PressureChangeType.calculated):

        if pressure_change_type == PressureChangeType.fixed_per_stage:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="pressure change due to friction",
            )
            def eq_pressure_change(b, t, x):
                return b.deltaP_channel[t] == b.dP_dx[t, x] * b.length

        else:

            @self.Constraint(
                self.flowsheet().config.time, doc="Total Pressure drop across channel"
            )
            def eq_pressure_change(b, t):
                return b.deltaP_channel[t] == sum(
                    b.dP_dx[t, x] * b.length / b.nfe for x in b.difference_elements
                )

    def initialize(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
        hold_state=True,
        initialize_guess=None,
        type=None,
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
            initialize_guess : a dict of guesses
        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """

        # Get inlet state if not provided
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="control_volume")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="control_volume")

        state_args = self._get_state_args(initialize_guess, state_args)
        state_args_properties_in = state_args["inlet"]

        source_flags = super().initialize(
            state_args=state_args_properties_in,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Differentiate between hot and cold channels for properties_out
        if type == "hot_ch":
            state_args_properties_out = state_args["hot_outlet"]

        elif type == "cold_ch":
            state_args_properties_out = state_args["cold_outlet"]
        else:
            raise ConfigurationError(
                "Either hot_ch or cold_ch must be set in the configuration."
            )

        if hasattr(self, "properties_interface"):
            state_args_interface = self._get_state_args_interface(
                state_args_properties_in, state_args_properties_out
            )

            state_args_interface = self._get_state_args_interface(
                state_args_properties_in, state_args_properties_out
            )
            self.properties_interface.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_interface,
            )

        if hasattr(self, "properties_vapor"):
            state_args_vapor = self._get_state_args_vapor(
                state_args_properties_in, state_args_properties_out
            )

            state_args_vapor = self._get_state_args_vapor(
                state_args_properties_in, state_args_properties_out
            )
            self.properties_vapor.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_vapor,
            )

        init_log.info("Initialization Complete")

        if hold_state:
            return source_flags
        else:
            self.release_state(source_flags, outlvl)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "area"):
            if iscale.get_scaling_factor(self.area) is None:
                iscale.set_scaling_factor(self.area, 100)

        if hasattr(self, "deltaP_channel"):
            for v in self.deltaP_channel.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, "dP_dx"):
            for v in self.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e-5)
