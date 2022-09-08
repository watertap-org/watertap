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

from pyomo.environ import (
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
    EnergyBalanceType,
    FlowDirection,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import (
    MembraneChannelMixin,
    PressureChangeType,
    CONFIG_Template,
)


@declare_process_block_class("MembraneChannel0DBlock")
class MembraneChannel0DBlockData(MembraneChannelMixin, ControlVolume0DBlockData):
    def _skip_element(self, x):
        return False

    # overwrite CV0D `add_geometry`
    def add_geometry(
        self, length_var=None, width_var=None, flow_direction=FlowDirection.forward
    ):
        """
        Method to create spatial domain and volume Var in ControlVolume.

        Args:
            length_var - (optional) external variable to use for the length of
                         the channel. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.

            width_var - (optional) external variable to use for the width of
                         the channel. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.
            flow_direction - argument indicating direction of material flow
                            relative to length domain. Valid values:
                                - FlowDirection.forward (default), flow goes
                                  from 0 to 1.
                                - FlowDirection.backward, flow goes from 1 to 0
        Returns:
            None
        """
        # Validate and create flow direction attribute, like 1D
        if flow_direction in (flwd for flwd in FlowDirection):
            self._flow_direction = flow_direction
        else:
            raise ConfigurationError(
                "{} invalid value for flow_direction "
                "argument. Must be a FlowDirection Enum.".format(self.name)
            )

        self._add_var_reference(length_var, "length", "length_var")
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

        # quack like a 1D model
        self.length_domain = Set(ordered=True, initialize=(0.0, 1.0))
        add_object_reference(self, "difference_elements", self.length_domain)

        self._set_nfe()

        if self._flow_direction == FlowDirection.forward:
            add_object_reference(
                self,
                "properties",
                {
                    **{
                        (t, 0.0): self.properties_in[t]
                        for t in self.flowsheet().config.time
                    },
                    **{
                        (t, 1.0): self.properties_out[t]
                        for t in self.flowsheet().config.time
                    },
                },
            )
        else:
            add_object_reference(
                self,
                "properties",
                {
                    **{
                        (t, 0.0): self.properties_out[t]
                        for t in self.flowsheet().config.time
                    },
                    **{
                        (t, 1.0): self.properties_in[t]
                        for t in self.flowsheet().config.time
                    },
                },
            )

        self._add_interface_stateblock(has_phase_equilibrium)

    def apply_transformation(self):
        pass

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        if pressure_change_type == PressureChangeType.fixed_per_stage:
            return

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if pressure_change_type == PressureChangeType.fixed_per_unit_length:
            # Pressure change equation when dP/dx = user-specified constant,
            self.dP_dx = Var(
                self.flowsheet().config.time,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta("pressure") * units_meta("length") ** -1,
                doc="pressure drop per unit length across channel",
            )

            @self.Constraint(
                self.flowsheet().config.time, doc="pressure change due to friction"
            )
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx[t] * b.length

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

        # TODO: this function needs to be changed for use on the permeate side
        state_args = self._get_state_args(initialize_guess, state_args)

        # intialize self.properties
        source_flags = self.properties_in.initialize(
            state_args=state_args["feed_side"],
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        self.properties_out.initialize(
            state_args=state_args["retentate"],
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )

        self.properties_interface.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["interface"],
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
