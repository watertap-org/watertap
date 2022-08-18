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
from idaes.core import (declare_process_block_class,
    EnergyBalanceType,
    FlowDirection,
    )
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
from watertap.core.membrane_channel_base import (
    MembraneChannelMixin, PressureChangeType, CONFIG_Template)


@declare_process_block_class("MembraneChannel0DBlock")
class MembraneChannel0DBlockData(MembraneChannelMixin, ControlVolume0DBlockData):

    # overwrite CV0D `add_geometry`
    def add_geometry(self, include_length_and_width=True):
        """
        Method to create spatial domain and volume Var in ControlVolume.

        Args:
            include_length_and_width - (optional) add a length and width
                variables to the membrane channel. Default: `True`

        Returns:
            None
        """

        if include_length_and_width:
            units_meta = self.config.property_package.get_metadata().get_derived_units
            self.length = Var(
                 initialize=10,
                 bounds=(0.1, 5e2),
                 domain=NonNegativeReals,
                 units=units_meta("length"),
                 doc="Effective membrane length",
                )
            self.width = Var(
                initialize=1,
                bounds=(1e-1, 1e3),
                domain=NonNegativeReals,
                units=units_meta("length"),
                doc="Membrane width",
            )

    def add_state_blocks(
        self, information_flow=FlowDirection.forward, has_phase_equilibrium=None
    ):
        """
        This method constructs the state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
        Returns:
            None
        """
        super().add_state_blocks(information_flow, has_phase_equilibrium)
        # quack like a 1D model
        self.length_domain = Set(ordered=True, initialize=(0.0, 1.0))
        add_object_reference(self, "difference_elements", self.length_domain)
        self.first_element = self.length_domain.first()

        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

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
            self.dP_dx= Var(
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
            raise ConfigurationError(f"Unrecognized pressure_change_type {pressure_change_type}")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "area"):
            if iscale.get_scaling_factor(self.area) is None:
                iscale.set_scaling_factor(self.area, 100)

        if hasattr(self, "dP_dx"):
            for v in self.dP_dx.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
