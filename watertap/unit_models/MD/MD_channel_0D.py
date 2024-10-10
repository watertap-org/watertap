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
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
import idaes.logger as idaeslog
from .MD_channel_base import (
    MDChannelMixin,
    PressureChangeType,
)

__author__ = "Elmira Shamlou"


@declare_process_block_class("MDChannel0DBlock")
class MDChannel0DBlockData(MDChannelMixin, ControlVolume0DBlockData):
    def _skip_element(self, x):
        return False

    def add_geometry(
        self, length_var=None, width_var=None, flow_direction=FlowDirection.forward
    ):

        self._flow_direction = flow_direction

        # If the length_var and width_var are provided, create references to them
        if length_var is not None:
            add_object_reference(self, "length", length_var)

        if width_var is not None:
            add_object_reference(self, "width", width_var)

    def add_state_blocks(
        self,
        has_phase_equilibrium=None,
    ):
        super().add_state_blocks(has_phase_equilibrium=has_phase_equilibrium)

        # quack like a 1D model
        self.length_domain = Set(ordered=True, initialize=(0.0, 1.0))
        add_object_reference(self, "difference_elements", self.length_domain)

        self._set_nfe()

        # Determine flow direction from the argument or from the configuration

        if self._flow_direction == FlowDirection.forward:
            properties_dict = {
                **{
                    (t, 0.0): self.properties_in[t]
                    for t in self.flowsheet().config.time
                },
                **{
                    (t, 1.0): self.properties_out[t]
                    for t in self.flowsheet().config.time
                },
            }
        elif self._flow_direction == FlowDirection.backward:
            properties_dict = {
                **{
                    (t, 0): self.properties_out[t] for t in self.flowsheet().config.time
                },
                **{(t, 1): self.properties_in[t] for t in self.flowsheet().config.time},
            }
        else:
            raise ConfigurationError(
                "FlowDirection must be set to FlowDirection.forward or FlowDirection.backward."
            )

        add_object_reference(self, "properties", properties_dict)

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        if pressure_change_type == PressureChangeType.fixed_per_stage:
            return

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if pressure_change_type == PressureChangeType.fixed_per_unit_length:
            # Pressure change equation when dP/dx = user-specified constant,
            self.dP_dx = Var(
                self.flowsheet().config.time,
                initialize=-5e4,
                bounds=(-2e5, None),
                domain=NegativeReals,
                units=units_meta("pressure") * units_meta("length") ** -1,
                doc="pressure drop per unit length across channel",
            )

        elif pressure_change_type == PressureChangeType.calculated:

            self.dP_dx = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta("pressure") * units_meta("length") ** -1,
                doc="Pressure drop per unit length of channel at inlet and outlet",
            )

    def _add_deltaP(self, pressure_change_type=PressureChangeType.calculated):
        if pressure_change_type == PressureChangeType.fixed_per_stage:
            return

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if pressure_change_type == PressureChangeType.fixed_per_unit_length:

            @self.Constraint(
                self.flowsheet().config.time, doc="pressure change due to friction"
            )
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx[t] * b.length

        elif pressure_change_type == PressureChangeType.calculated:

            @self.Constraint(
                self.flowsheet().config.time, doc="Total Pressure drop across channel"
            )
            def eq_pressure_change(b, t):
                return b.deltaP[t] == sum(
                    b.dP_dx[t, x] * b.length / b.nfe for x in b.length_domain
                )

        else:
            raise ConfigurationError(
                f"Unrecognized pressure_change_type {pressure_change_type}"
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

        source_flags = self.properties_in.initialize(
            state_args=state_args_properties_in,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Differentiate between hot and cold channels for properties_out
        if type == "hot_ch":
            state_args_properties_out = state_args["hot_outlet"]

            self.properties_out.initialize(
                state_args=state_args_properties_out,
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
            )
        elif type == "cold_ch":
            state_args_properties_out = state_args["cold_outlet"]

            self.properties_out.initialize(
                state_args=state_args_properties_out,
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
            )
        else:
            raise ConfigurationError(
                "Either hot_ch or cold_ch must be set in the configuration."
            )
        if hasattr(self, "properties_interface"):
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

        if hasattr(self, "dP_dx"):
            for v in self.dP_dx.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
