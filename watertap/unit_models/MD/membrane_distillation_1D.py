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

# Import Pyomo libraries
from pyomo.environ import Constraint
from pyomo.common.config import Bool, ConfigDict, ConfigValue, ConfigBlock, In
from idaes.core import FlowDirection

from .MD_channel_base import (
    ConcentrationPolarizationType,
    TemperaturePolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
    FrictionFactor,
)
from idaes.core import (
    FlowDirection,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    useDefault,
    DistributedVars,
)

from idaes.core import declare_process_block_class
from .MD_channel_base import (
    TemperaturePolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from .MD_channel_1D import MDChannel1DBlock
from .membrane_distillation_base import (
    MembraneDistillationBaseData,
    MDconfigurationType,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

__author__ = "Elmira Shamlou"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("MembraneDistillation1D")
class MembraneDistillationData(MembraneDistillationBaseData):
    """
    Standard DCMD Unit Model Class:
    - one dimensional model
    - steady state only
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not.
    **default** - False. Membrane units do not yet support dynamic
    behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Membrane units do not have defined volume, thus
    this must be False.""",
        ),
    )
    _CONFIG_Template = ConfigBlock()

    _CONFIG_Template.declare(
        "dynamic",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not.
    **default** - False. Membrane units do not yet support dynamic
    behavior.""",
        ),
    )

    _CONFIG_Template.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Membrane units do not have defined volume, thus
    this must be False.""",
        ),
    )

    _CONFIG_Template.declare(
        "property_package",
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
    _CONFIG_Template.declare(
        "property_package_args",
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

    _CONFIG_Template.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault** - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )

    _CONFIG_Template.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed.
    **default** - useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault** - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )

    _CONFIG_Template.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
    **default** - MomentumBalanceType.pressureTotal.
    **Valid values:** {
    **MomentumBalanceType.none** - exclude momentum balances,
    **MomentumBalanceType.pressureTotal** - single pressure balance for material,
    **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
    **MomentumBalanceType.momentumTotal** - single momentum balance for material,
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    _CONFIG_Template.declare(
        "flow_direction",
        ConfigValue(
            default=FlowDirection.forward,
            domain=In(FlowDirection),
            description="Direction of flow",
            doc="""
            Options for the direction of flow:

            **default** - ``FlowDirection.forward``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``FlowDirection.forward``", "Flow is in the forward direction"
            "``FlowDirection.backward``", "Flow is in the backward direction"
            """,
        ),
    )

    _CONFIG_Template.declare(
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

    _CONFIG_Template.declare(
        "concentration_polarization_type",
        ConfigValue(
            default=ConcentrationPolarizationType.calculated,
            domain=In(ConcentrationPolarizationType),
            description="External concentration polarization effect in RO",
            doc="""
            Options to account for concentration polarization.

            **default** - ``ConcentrationPolarizationType.calculated``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``ConcentrationPolarizationType.none``", "Simplifying assumption to ignore concentration polarization"
            "``ConcentrationPolarizationType.fixed``", "Specify an estimated value for the concentration polarization modulus"
            "``ConcentrationPolarizationType.calculated``", "Allow model to perform calculation of membrane-interface concentration"
        """,
        ),
    )

    _CONFIG_Template.declare(
        "mass_transfer_coefficient",
        ConfigValue(
            default=MassTransferCoefficient.calculated,
            domain=In(MassTransferCoefficient),
            description="Mass transfer coefficient in RO feed channel",
            doc="""
            Options to account for mass transfer coefficient.

            **default** - ``MassTransferCoefficient.calculated``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``MassTransferCoefficient.none``", "Mass transfer coefficient not used in calculations"
            "``MassTransferCoefficient.fixed``", "Specify an estimated value for the mass transfer coefficient in the feed channel"
            "``MassTransferCoefficient.calculated``", "Allow model to perform calculation of mass transfer coefficient"
        """,
        ),
    )

    _CONFIG_Template.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
        ),
    )

    _CONFIG_Template.declare(
        "pressure_change_type",
        ConfigValue(
            default=PressureChangeType.fixed_per_stage,
            domain=In(PressureChangeType),
            description="Pressure change term construction flag",
            doc="""
        Indicates what type of pressure change calculation will be made. To use any of the
        ``pressure_change_type`` options to account for pressure drop, the configuration keyword
        ``has_pressure_change`` must also be set to ``True``. Also, if a value is specified for pressure
        change, it should be negative to represent pressure drop.

        **default** - ``PressureChangeType.fixed_per_stage`` 

        .. csv-table::
            :header: "Configuration Options", "Description"
        
            "``PressureChangeType.fixed_per_stage``", "Specify an estimated value for pressure drop across the membrane feed channel"
            "``PressureChangeType.fixed_per_unit_length``", "Specify an estimated value for pressure drop per unit length across the membrane feed channel"
            "``PressureChangeType.calculated``", "Allow model to perform calculation of pressure drop across the membrane feed channel"
        """,
        ),
    )

    _CONFIG_Template.declare(
        "friction_factor",
        ConfigValue(
            default=FrictionFactor.flat_sheet,
            domain=In(FrictionFactor),
            description="Darcy friction factor correlation",
            doc="""
            Options to account for friction factor correlations.

            **default** - ``FrictionFactor.flat_sheet`` 

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``FrictionFactor.flat_sheet``", "Friction factor correlation for flat-sheet membrane modules"
            "``FrictionFactor.spiral_wound``", "Friction factor correlation for spiral-wound membranes"
        """,
        ),
    )

    _CONFIG_Template.declare(
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
    _CONFIG_Template.declare(
        "property_package_args_vapor",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}""",
        ),
    )

    CONFIG.declare("hot_ch", _CONFIG_Template(doc="hot channel config arguments"))
    CONFIG.declare("cold_ch", _CONFIG_Template(doc="cold channel config arguments"))
    CONFIG.declare("gap_ch", _CONFIG_Template(doc="cold channel config arguments"))

    # Common config args for both sides

    CONFIG.declare(
        "MD_configuration_Type",
        ConfigValue(
            default=MDconfigurationType.DCMD,
            domain=In(MDconfigurationType),
            description="Options for selecting membrane distillation process configurations",
            doc="""
            Options for selecting membrane distillation process configurations
         **default** - ``TemperaturePolarizationType.calculated``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``MDconfigurationType.DCMD``", "Direct Contact Membrane Distillation"
            "``MDconfigurationType.VMD``", "Vacuum Membrane Distillation"
            "``MDconfigurationType.PGMD_CGMD``", "Permeate Gap or Coductive Gap Membrane Distillation"
            "``MDconfigurationType.AGMD``", "Air Gap Membrane Distillation"
        """,
        ),
    )

    CONFIG.declare(
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
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See
        Pyomo documentation for supported transformations.""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transformating domain. See
        Pyomo documentation for supported schemes.""",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
        domain (default=20)""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
        discretizing length domain (default=3)""",
        ),
    )

    def _make_MD_channel_control_volume(
        self, name_ch, common_config, individual_config
    ):

        if self.config.transformation_method is useDefault:
            _log.warning(
                "Discretization method was "
                "not specified for the "
                "membrane distillation module. "
                "Defaulting to finite "
                "difference method."
            )
            self.config.transformation_method = "dae.finite_difference"

        if self.config.transformation_scheme is useDefault:
            _log.warning(
                "Discretization scheme was "
                "not specified for the "
                "membrane distillation module."
                "Defaulting to backward finite "
                "difference."
            )
            self.config.transformation_scheme = "BACKWARD"

        if not isinstance(name_ch, str):
            raise TypeError(
                f"{name_ch} is not a string. Please provide a string for the name_ch argument."
            )

        # Build membrane channel control volume
        self.add_component(
            name_ch,
            MDChannel1DBlock(
                dynamic=False,
                has_holdup=False,
                area_definition=common_config.area_definition,
                property_package=individual_config.property_package,
                property_package_args=individual_config.property_package_args,
                transformation_method=common_config.transformation_method,
                transformation_scheme=common_config.transformation_scheme,
                finite_elements=common_config.finite_elements,
                collocation_points=common_config.collocation_points,
            ),
        )

        channel = getattr(self, name_ch)
        if not hasattr(self, "length") and not hasattr(self, "width"):
            self._add_length_and_width()
        channel.add_geometry(
            flow_direction=individual_config.flow_direction,
            length_var=self.length,
            width_var=self.width,
        )

        if not hasattr(self, "eq_area"):
            add_eq_area = True
        else:
            add_eq_area = False
        self._add_area(include_constraint=add_eq_area)

    def _add_mass_transfer(self):
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.cold_ch.property_package.phase_list,
            self.config.cold_ch.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, x, p, j):

            if p == "Liq":
                if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                    return (
                        b.cold_ch.mass_transfer_term[t, x, p, j]
                        == -b.hot_ch.mass_transfer_term[t, x, p, j]
                    )
                elif self.config.MD_configuration_Type == MDconfigurationType.PGMD_CGMD:
                    return (
                        b.gap_ch.mass_transfer_term[t, x, "Liq", "H2O"]
                        == -b.hot_ch.mass_transfer_term[t, x, "Liq", "H2O"]
                    )
                elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                    # b.cold_ch.mass_transfer_term[t, x, "Liq", j].fix(0)
                    return Constraint.Skip
            else:
                if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                    b.cold_ch.mass_transfer_term[t, x, p, j].fix(0)
                    return Constraint.Skip
                elif self.config.MD_configuration_Type == MDconfigurationType.PGMD_CGMD:
                    b.gap_ch.mass_transfer_term[t, x, p, j].fix(0)
                    return Constraint.Skip
                elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                    return (
                        b.cold_ch.mass_transfer_term[t, x, p, j]
                        == -b.hot_ch.mass_transfer_term[t, x, "Liq", j]
                    )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.hot_ch.property_package.phase_list,
            self.config.hot_ch.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, x, p, j):
            if j == "H2O":
                return (
                    b.hot_ch.mass_transfer_term[t, x, p, j]
                    == -b.width * b.flux_mass[t, x]
                )
            else:
                b.hot_ch.mass_transfer_term[t, x, p, j].fix(0)
                return Constraint.Skip

    def _add_heat_transfer(self):
        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Conductive heat transfer to cold channel",
        )
        def eq_conductive_heat_transfer_hot(b, t, x):
            return b.hot_ch.heat[t, x] == -b.width * b.flux_conduction_heat[t, x]

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Conductive heat transfer to cold channel",
        )
        def eq_conductive_heat_transfer_cold(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.DCMD:
                return b.cold_ch.heat[t, x] == -b.hot_ch.heat[t, x]
            elif self.config.MD_configuration_Type == MDconfigurationType.PGMD_CGMD:
                return (
                    b.cold_ch.heat[t, x]
                    == -b.hot_ch.heat[t, x]
                    - b.hot_ch.enthalpy_transfer[t, x]
                    - b.gap_ch.enthalpy_transfer[t, x]
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.VMD:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Connecting the cold channel conductive heat transfer to conductive heat across the gap",
        )
        def eq_conductive_heat_transfer_gap(b, t, x):
            if self.config.MD_configuration_Type == MDconfigurationType.PGMD_CGMD:
                return (
                    b.cold_ch.heat[t, x] == b.flux_conduction_heat_gap[t, x] * b.width
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Enthalpy heat transfer from the hot channel",
        )
        def eq_enthalpy_transfer_hot(b, t, x):
            return b.hot_ch.enthalpy_transfer[t, x] == -b.width * b.flux_enth_hot[t, x]

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Enthalpy heat transfer to the cold channel",
        )
        def eq_enthalpy_transfer_cold(b, t, x):
            if self.config.MD_configuration_Type in [
                MDconfigurationType.DCMD,
                MDconfigurationType.VMD,
            ]:
                return (
                    b.cold_ch.enthalpy_transfer[t, x]
                    == b.width * b.flux_enth_cold[t, x]
                )
            elif self.config.MD_configuration_Type == MDconfigurationType.PGMD_CGMD:
                return (
                    b.gap_ch.enthalpy_transfer[t, x] == b.width * b.flux_enth_cold[t, x]
                )
            else:
                return Constraint.Skip
