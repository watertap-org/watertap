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

import collections
import weakref

from pyomo.environ import (
    NegativeReals,
    NonNegativeReals,
    Set,
    Var,
)
from pyomo.core.base.indexed_component import slicer_types
from idaes.core import (
    declare_process_block_class,
    FlowDirection,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import (
    MembraneChannelMixin,
    PressureChangeType,
    CONFIG_Template as Base_CONFIG_Template,
)

CONFIG_Template = Base_CONFIG_Template()


class Not0Or1Error(KeyError):
    pass


class _0DPropertyHelper(collections.abc.Mapping):
    """
    Class to make blk.properties_in and blk.properties_out
    like blk.properties from a 1D model.
    """

    def __init__(self, blk, reverse=False):
        self._blk_ref = weakref.ref(blk)
        if reverse:
            self._get_property_x = self._get_property_x_backward
        else:
            self._get_property_x = self._get_property_x_forward

    def __len__(self):
        return len(self.blk.properties_in) + len(self.blk.properties_out)

    def __iter__(self):
        for t in self._get_property_x(0):
            yield (t, 0)
        for t in self._get_property_x(1):
            yield (t, 1)

    def __contains__(self, index):
        idx0 = index[0]
        idx1 = index[1]
        try:
            return idx0 in self._get_property_x(idx1)
        except Not0Or1Error:
            return False

    def __getitem__(self, index):
        if index is Ellipsis:
            return (_ for _ in self._props_in_out(index))
        try:
            index_len = len(index)
        except:
            raise KeyError(index)
        if index_len != 2:
            raise KeyError(index)
        idx0 = index[0]
        idx1 = index[1]
        if type(idx1) in slicer_types:
            if isinstance(idx1, slice):
                if (
                    idx1.start is not None
                    or idx1.stop is not None
                    or idx1.step is not None
                ):
                    raise IndexError("Indexed components only support simple slices")
            return (_ for _ in self._props_in_out(idx0))
        try:
            props = self._get_property_x(idx1)
            return props[idx0]
        except Not0Or1Error:
            raise KeyError(index)

    def keys(self):
        for t in self._get_property_x(0).keys():
            yield (t, 0)
        for t in self._get_property_x(1).keys():
            yield (t, 1)

    def items(self):
        for t, v in self._get_property_x(0).items():
            yield (t, 0, v)
        for t, v in self._get_property_x(1).items():
            yield (t, 1, v)

    def values(self):
        yield from self._get_property_x(0).values()
        yield from self._get_property_x(1).values()

    @property
    def blk(self):
        return self._blk_ref()

    def _get_property_x_forward(self, x):
        if x == 0:
            return self.blk.properties_in
        elif x == 1:
            return self.blk.properties_out
        else:
            raise Not0Or1Error

    def _get_property_x_backward(self, x):
        if x == 1:
            return self.blk.properties_in
        elif x == 0:
            return self.blk.properties_out
        else:
            raise Not0Or1Error

    def _props_in_out(self, idx0):
        if type(idx0) in slicer_types:
            yield from self._get_property_x(0)[idx0]
            yield from self._get_property_x(1)[idx0]
        else:
            yield self._get_property_x(0)[idx0]
            yield self._get_property_x(1)[idx0]


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
        # volume will be added be calling add_geometry from ControlVolume0D
        if not hasattr(self, "length") and not hasattr(self, "width"):
            self._add_var_reference(length_var, "length", "length_var")
            self._add_var_reference(width_var, "width", "width_var")

        if hasattr(self, "length") and hasattr(self, "width"):
            super().add_geometry()
            units_meta = self.config.property_package.get_metadata().get_derived_units

            if not hasattr(self, "channel_height"):
                self.channel_height = Var(
                    initialize=1e-3,
                    bounds=(1e-4, 5e-3),
                    domain=NonNegativeReals,
                    units=units_meta("length"),
                    doc="membrane-channel height",
                )

            # # TODO: negating spacer volume for now as it can be assumed negligible (Park et al., 2020: https://doi.org/10.1016/j.desal.2020.114625). Revisit.
            # # Note: considering that this volume relationship holds for both flat_sheet and spiral_wound types
            @self.Constraint(
                self.flowsheet().config.time, doc="Membrane-channel volume"
            )
            def eq_volume(b, t):
                return b.volume[t] == b.length * b.width * b.channel_height

        # Validate and create flow direction attribute, like 1D
        if flow_direction in (flwd for flwd in FlowDirection):
            self._flow_direction = flow_direction
        else:
            raise ConfigurationError(
                "{} invalid value for flow_direction "
                "argument. Must be a FlowDirection Enum.".format(self.name)
            )

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
            self.properties = _0DPropertyHelper(self, reverse=False)
        elif self._flow_direction == FlowDirection.backward:
            self.properties = _0DPropertyHelper(self, reverse=True)
        else:
            raise ConfigurationError(
                "FlowDirection must be set to FlowDirection.forward or FlowDirection.backward."
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
                bounds=(-2e5, None),
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
        state_args_properties_in = state_args["feed_side"]
        if self._flow_direction == FlowDirection.forward:
            state_args_properties_out = state_args["retentate"]
        else:
            state_args_properties_out = state_args["permeate"]

        source_flags = self.properties_in.initialize(
            state_args=state_args_properties_in,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        self.properties_out.initialize(
            state_args=state_args_properties_out,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )

        state_args_interface = self._get_state_args_interface(
            initialize_guess, state_args_properties_in, state_args_properties_out
        )
        self.properties_interface.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_interface,
        )

        init_log.info("Initialization Complete")

        if hold_state:
            return source_flags
        else:
            self.release_state(source_flags, outlvl)

    def calculate_scaling_factors(self):
        # set volume scale factor first otherwise ControlVolume0D default will be set
        if hasattr(self, "volume"):
            for v in self.volume.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e3)

        super().calculate_scaling_factors()

        if hasattr(self, "area"):
            if iscale.get_scaling_factor(self.area) is None:
                iscale.set_scaling_factor(self.area, 100)

        if hasattr(self, "dP_dx"):
            for v in self.dP_dx.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
