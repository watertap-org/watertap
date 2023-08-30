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
)

from idaes.core import declare_process_block_class
from idaes.core.util import scaling as iscale
from .MD_channel_base import (
    TemperaturePolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from .MD_channel_0D import MDChannel0DBlock
from .membrane_distillation_base import (
    MembraneDistillationBaseData,
)
from idaes.core.util.config import is_physical_parameter_block, DefaultBool

# __author__ = "Elmira Shamlou"


@declare_process_block_class("MembraneDistillation0D")
class MembraneDistillationData(MembraneDistillationBaseData):
    """
    Standard DCMD Unit Model Class:
    - zero dimensional model
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

    def _make_MD_channel_control_volume(self, name_ch, config):

        if not isinstance(name_ch, str):
            raise TypeError(
                f"{name_ch} is not a string. Please provide a string for the name_ch argument."
            )

        # Build membrane channel control volume
        self.add_component(
            name_ch,
            MDChannel0DBlock(
                dynamic=False,
                has_holdup=False,
                property_package=config.property_package,
                property_package_args=config.property_package_args,
            ),
        )

        channel = getattr(self, name_ch)

        if (
            (config.pressure_change_type != PressureChangeType.fixed_per_stage)
            or (config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
            or (
                config.temperature_polarization_type
                == TemperaturePolarizationType.calculated
            )
        ):

            if not hasattr(self, "length") and not hasattr(self, "width"):
                self._add_length_and_width()
            channel.add_geometry(
                length_var=self.length,
                width_var=self.width,
            )
            if not hasattr(self, "eq_area"):
                add_eq_area = True
            else:
                add_eq_area = False
            self._add_area(include_constraint=add_eq_area)
        else:
            channel.add_geometry(length_var=None, width_var=None)
            self._add_area(include_constraint=False)

    def _add_mass_transfer(self):

        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )
        """
        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(
                self.hot_ch.properties_in[t].get_material_flow_terms("Liq", "H2O")
                * self.recovery_mass_phase_comp[t]
            )

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.hot_ch.property_package.phase_list,
            self.config.hot_ch.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="Mass transfer to cold channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.hot_ch.property_package.phase_list,
            self.config.hot_ch.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(self, t, p, j):
            return (
                self.mass_transfer_phase_comp[t, p, j]
                == -self.hot_ch.mass_transfer_term[t, p, j]
            )
            """

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.cold_ch.property_package.phase_list,
            self.config.cold_ch.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, p, j):
            return (
                b.cold_ch.mass_transfer_term[t, p, j]
                == -b.hot_ch.mass_transfer_term[t, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.cold_ch.property_package.phase_list,
            self.config.cold_ch.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, p, j):
            if j == "H2O":
                return (
                    b.cold_ch.mass_transfer_term[t, p, j] == b.area * b.flux_mass_avg[t]
                )
            else:
                b.cold_ch.mass_transfer_term[t, p, j].fix(0)
                return Constraint.Skip

    def _add_heat_transfer(self):
        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )

        # Enthalpy Transfer Initialization
        # todo: change the following to general format to different types of vapor
        """
                def enthalpy_transfer_initialize(b, t):
                    return value(
                        self.hot_ch.properties_in[t].get_material_flow_terms("Liq", "H2O")
                        * self.recovery_mass_phase_comp[t, "Liq", "H2O"]
                        * self.hot_ch.properties_vapor[t, 0].enth_mass_phase["Vap"]
                    )

                self.enthalpy_transfer_var = Var(
                    self.flowsheet().config.time,
                    initialize=enthalpy_transfer_initialize,
                    domain=NonNegativeReals,
                    units=pyunits.J * pyunits.s**-1,
                    doc="Intermediate variable for enthalpy transfer",
                )

                @self.Constraint(self.flowsheet().config.time)
                def eq_enthalpy_transfer_hot_initialize(b, t):
                    return b.enthalpy_transfer_var[t] == b.hot_ch.enthalpy_transfer[t]

                @self.Constraint(self.flowsheet().config.time)
                def eq_enthalpy_transfer_cold_initialize(b, t):
                    return b.enthalpy_transfer_var[t] == b.cold_ch.enthalpy_transfer[t]

                # Heat Transfer Initialization
            """

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Conductive heat transfer to cold channel",
        )
        def eq_conductive_heat_transfer_hot(b, t):
            return b.hot_ch.heat[t] == -b.area * b.flux_conduction_heat_avg[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Conductive heat transfer to cold channel",
        )
        def eq_conductive_heat_transfer_cold(b, t):
            return b.cold_ch.heat[t] == -b.hot_ch.heat[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Enthalpy heat transfer from the hot channel",
        )
        def eq_enthalpy_transfer_hot(b, t):
            return b.hot_ch.enthalpy_transfer[t] == -b.area * b.flux_enth_hot_avg[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Enthalpy heat transfer to the cold channel",
        )
        def eq_enthalpy_transfer_cold(b, t):
            return b.cold_ch.enthalpy_transfer[t] == b.area * b.flux_enth_cold_avg[t]

    def calculate_scaling_factors(self):
        # First, we'll call the `calculate_scaling_factors` method of the super class
        super().calculate_scaling_factors()

        # Next, we'll define the scaling factors for the variables and constraints
        # that we've defined in this class.
        """
        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(
                self.hot_ch.properties_in[t].get_material_flow_terms(p, j)
            )
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.hot_ch.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.cold_ch.mass_transfer_term[t, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
        """
        # Get scaling factor from vapor enthalpy flow

        # Scaling factors for heat transfer
        for t in self.flowsheet().config.time:
            sf_vap = 1e4

            # If scaling factor has not been set for heat transfer, set it to sf_vap
            if iscale.get_scaling_factor(self.hot_ch.heat[t]) is None:
                iscale.set_scaling_factor(self.hot_ch.heat[t], sf_vap)
            if iscale.get_scaling_factor(self.cold_ch.heat[t]) is None:
                iscale.set_scaling_factor(self.cold_ch.heat[t], sf_vap)
            if iscale.get_scaling_factor(self.hot_ch.enthalpy_transfer[t]) is None:
                iscale.set_scaling_factor(self.hot_ch.enthalpy_transfer[t], sf_vap)
            if iscale.get_scaling_factor(self.cold_ch.enthalpy_transfer[t]) is None:
                iscale.set_scaling_factor(self.cold_ch.enthalpy_transfer[t], sf_vap)
            if iscale.get_scaling_factor(self.cold_ch.enthalpy_transfer[t]) is None:
                iscale.set_scaling_factor(self.cold_ch.enthalpy_transfer_var[t], sf_vap)
            # if iscale.get_scaling_factor(self.enthalpy_transfer_var[t]) is None:
            # iscale.set_scaling_factor(self.enthalpy_transfer_var[t], sf_vap)
