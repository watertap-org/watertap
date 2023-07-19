#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from enum import Enum, auto

from pyomo.common.config import ConfigValue, In, ConfigDict
from pyomo.environ import (
    NonNegativeReals,
    Var,
    units as pyunits,
)

from idaes.core import (
    FlowDirection,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    useDefault,
)

from idaes.core.util import scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference

from watertap.core.membrane_channel_base import (
    MembraneChannelMixin,
    CONFIG_Template as Base_CONFIG_Template,
)

CONFIG_Template = Base_CONFIG_Template()


class TemperaturePolarizationType(Enum):
    """
    none: no temperature polarization
    fixed: temperature polarization modulus is a user specified value
    calculated: calculate temperature polarization (temperature at membrane interface)
    """

    none = auto()
    fixed = auto()
    calculated = auto()


CONFIG_Template.declare(
    "temperature_polarization_type",
    ConfigValue(
        default=TemperaturePolarizationType.calculated,
        domain=In(TemperaturePolarizationType),
        description="External temperature polarization effect",
        doc="""
        Options to account for temperature polarization.

        **default** - ``TemperaturePolarizationType.calculated``

    .. csv-table::
        :header: "Configuration Options", "Description"

        "``TemperaturePolarizationType.none``", "Simplifying assumption to ignore temperature polarization"
        "``TemperaturePolarizationType.fixed``", "Specify an estimated value for the temperature polarization modulus"
        "``TemperaturePolarizationType.calculated``", "Allow model to perform calculation of membrane-interface temperature"
    """,
    ),
)

CONFIG_Template.declare(
    "property_package_vapor",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
    ),
)

CONFIG_Template.declare(
    "property_package_args_vapor",
    ConfigDict(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigDict with arguments to be passed to a property block(s)
and used when constructing these.
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)


class TemperaturePolarizationMixin(MembraneChannelMixin):
    def add_interface_isothermal_conditions(self):
        pass

    def add_control_volume_isothermal_conditions(self):
        pass

    def add_temperature_polarization(
        self,
        temperature_polarization_type=TemperaturePolarizationType.none,
    ):

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if temperature_polarization_type == TemperaturePolarizationType.none:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="No temperature polarization",
            )
            def eq_no_temp_pol(b, t, x):
                return (
                    b.properties_interface[t, x].temperature
                    == b.properties[t, x].temperature
                )

            return self.eq_no_temp_pol

        elif temperature_polarization_type not in (
            TemperaturePolarizationType.fixed,
            TemperaturePolarizationType.calculated,
        ):
            raise ConfigurationError(
                f"Unrecognized temperature_polarization_type {temperature_polarization_type}"
            )

        self.h_conv = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=5e-5,
            bounds=(1e-6, 1e-3),
            domain=NonNegativeReals,
            units=units_meta("power")
            * units_meta("length") ** -2
            * units_meta("temperature_difference") ** -1,
            doc="Convective heat transfer coefficient",
        )

        if temperature_polarization_type == TemperaturePolarizationType.calculated:
            self._add_calculated_convective_heat_transfer_coefficient()

        return self.h_conv

    def _add_calculated_convective_heat_transfer_coefficient(self):
        self._add_calculated_pressure_change_mass_transfer_components()

        self.N_Pr = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=5e2,
            bounds=(1, 50),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Prandtl number in membrane channel",
        )
        self.N_Nu = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=1e2,
            bounds=(1, 3e2),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Nusselt number in membrane channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Convective heat transfer coefficient",
        )
        def eq_h_conv(b, t, x):
            return (
                b.h_conv[t, x]
                == b.properties[t, x].therm_cond_phase["Liq"] * b.N_Nu[t, x] / b.dh
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Nusselt number",
        )
        def eq_N_Nu(b, t, x):
            return b.N_Nu[t, x] == 0.162 * (b.N_Re[t, x] * b.N_Pr[t, x]) ** 0.656

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Prandtl number",
        )
        def eq_N_Pr(b, t, x):
            return (
                b.N_Pr[t, x]
                * b.properties[t, x].visc_d_phase["Liq"]
                * b.properties[t, x].cp_mass_phase["Liq"]
                == b.properties[t, x].therm_cond_phase["Liq"]
            )

        return self.eq_h_conv

    def _add_vapor_stateblock(self, has_phase_equilibrium=None):
        """
        This method constructs the vapor state blocks for the
        control volume.

        Args:
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
        Returns:
            None
        """
        tmp_dict = dict(**self.config.property_package_args_vapor)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets

        self.properties_vapor = self.config.property_package_vapor.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of vapor phase",
            **tmp_dict,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "N_Pr"):
            for v in self.N_Pr.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "N_Nu"):
            for v in self.N_Nu.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "h_conv"):
            for v in self.h_conv.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

    def validate_membrane_config_args(self):

        super().validate_membrane_config_args()

        if self.config.temperature_polarization_type not in (
            TemperaturePolarizationType.none,
            TemperaturePolarizationType.fixed,
            TemperaturePolarizationType.calculated,
        ):
            raise ConfigurationError(
                f"Unrecognized temperature_polarization_type {self.config.temperature_polarization_type}"
                "\nValid options are TemperaturePolarizationType.none, TemperaturePolarizationType.fixed, "
                "or TemperaturePolarizationType.calculated"
            )
