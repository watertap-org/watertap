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

from copy import deepcopy
from enum import Enum, auto

from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import (
    Constraint,
    Expression,
    Param,
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


class ConcentrationPolarizationType(Enum):
    """
    none: no concentration polarization
    fixed: concentration polarization modulus is a user specified value
    calculated: calculate concentration polarization (concentration at membrane interface)
    """

    none = auto()
    fixed = auto()
    calculated = auto()


class MassTransferCoefficient(Enum):
    """
    none: mass transfer coefficient not utilized for concentration polarization effect
    fixed: mass transfer coefficient is a user specified value
    calculated: mass transfer coefficient is calculated
    """

    none = auto()
    fixed = auto()
    calculated = auto()
    # TODO: add option for users to define their own relationship?


class TransportModel(Enum):
    """
    SD: Solvent and solute mass flux is calculated using the Solution-Diffusion model
    SKK: Solvent and solute mass flux is calculated using the Spiegler-Kedem-Katchalsky model
    """

    SD = auto()
    SKK = auto()


class PressureChangeType(Enum):
    """
    fixed_per_stage: pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length: pressure drop per unit length across membrane channel is a user-specified value
    calculated: pressure drop across membrane channel is calculated
    """

    fixed_per_stage = auto()
    fixed_per_unit_length = auto()
    calculated = auto()


class FrictionFactor(Enum):
    """
    flat_sheet: Darcy's friction factor correlation by Guillen & Hoek
    spiral_wound: Darcy's friction factor correlation by Schock & Miquel, 1987
    """

    flat_sheet = auto()
    spiral_wound = auto()


CONFIG_Template = ConfigDict()

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
    "transport_model",
    ConfigValue(
        default=TransportModel.SD,
        domain=In(TransportModel),
        description="Mass transfer model in RO feed channel",
        doc="""
        Options to account for mass transfer model.

        **default** - ``TransportModel.SD``

        "``TransportModel.SD``", "Solution-diffusion model for describing water and salt transport for most membrane types"
        "``TransportModel.SKK``", "Speigler-Kedem-Katchalsky model for describing water and salt transport"

        """,
    ),
)

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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

CONFIG_Template.declare(
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


class MembraneChannelMixin:
    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        raise NotImplementedError()

    def _skip_element(self, x):
        raise NotImplementedError()

    def _add_var_reference(self, pyomo_var, reference_name, param_name):
        if pyomo_var is not None:
            # Validate pyomo_var and add a reference
            if not isinstance(pyomo_var, (Var, Param, Expression)):
                raise ConfigurationError(
                    f"{self.name} {param_name} must be a Pyomo Var, Param or "
                    "Expression."
                )
            elif pyomo_var.is_indexed():
                raise ConfigurationError(
                    f"{self.name} {param_name} must be a scalar (unindexed) "
                    "component."
                )
            add_object_reference(self, reference_name, pyomo_var)

    def _set_nfe(self):
        self.first_element = self.length_domain.first()
        self.last_element = self.length_domain.last()

        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

    def add_total_pressure_balances(
        self,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        custom_term=None,
        friction_factor=FrictionFactor.flat_sheet,
    ):
        super().add_total_pressure_balances(
            has_pressure_change=has_pressure_change, custom_term=custom_term
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Pressure at interface",
        )
        def eq_equal_pressure_interface(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return b.properties_interface[t, x].pressure == b.properties[t, x].pressure

        if has_pressure_change:
            self._add_pressure_change(pressure_change_type=pressure_change_type)

        if pressure_change_type == PressureChangeType.calculated:
            self._add_calculated_pressure_change(friction_factor=friction_factor)

    def add_interface_isothermal_conditions(self):

        # ==========================================================================
        # Bulk and interface connections on the feed-side
        # TEMPERATURE
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Temperature at interface",
        )
        def eq_equal_temp_interface(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.properties[t, x].temperature
                == b.properties_interface[t, x].temperature
            )

    def add_control_volume_isothermal_conditions(self):

        ## ==========================================================================
        # Feed-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isothermal assumption for feed channel",
        )
        def eq_feed_isothermal(b, t, x):
            if x == b.length_domain.first():
                return Constraint.Skip
            return (
                b.properties[t, b.length_domain.first()].temperature
                == b.properties[t, x].temperature
            )

    def add_extensive_flow_to_interface(self):
        # VOLUMETRIC FLOWRATE
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Volumetric flow at interface of inlet",
        )
        def eq_equal_flow_vol_interface(b, t, x):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.properties_interface[t, x].flow_vol_phase["Liq"]
                == b.properties[t, x].flow_vol_phase["Liq"]
            )

    def add_concentration_polarization(
        self,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    ):

        solute_set = self.config.property_package.solute_set
        units_meta = self.config.property_package.get_metadata().get_derived_units

        if concentration_polarization_type == ConcentrationPolarizationType.none:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                doc="Unit concentration polarization modulus",
            )
            def eq_cp_modulus(b, t, x, j):
                if b._skip_element(x):
                    return Constraint.Skip
                return (
                    b.properties_interface[t, x].conc_mass_phase_comp["Liq", j]
                    == b.properties[t, x].conc_mass_phase_comp["Liq", j]
                )

            return self.eq_cp_modulus

        if concentration_polarization_type not in (
            ConcentrationPolarizationType.fixed,
            ConcentrationPolarizationType.calculated,
        ):
            raise ConfigurationError(
                f"Unrecognized concentration_polarization_type {concentration_polarization_type}"
            )

        self.cp_modulus = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=1.1,
            bounds=(0.1, 3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Concentration polarization modulus",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            doc="Concentration polarization modulus",
        )
        def eq_cp_modulus(b, t, x, j):  # pylint: disable=function-redefined
            if b._skip_element(x):
                return Constraint.Skip
            return (
                self.properties_interface[t, x].conc_mass_phase_comp["Liq", j]
                == self.properties[t, x].conc_mass_phase_comp["Liq", j]
                * self.cp_modulus[t, x, j]
            )

        if concentration_polarization_type == ConcentrationPolarizationType.calculated:
            if mass_transfer_coefficient == MassTransferCoefficient.none:
                raise ConfigurationError()

            # mass_transfer_coefficient is either calculated or fixed
            self.K = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta("length") * units_meta("time") ** -1,
                doc="Membrane channel mass transfer coefficient",
            )

            if mass_transfer_coefficient == MassTransferCoefficient.calculated:
                self._add_calculated_mass_transfer_coefficient()

        return self.eq_cp_modulus

    def add_expressions(self):
        """
        Generate expressions for additional results desired for full report
        """

        if hasattr(self, "N_Re"):

            @self.Expression(
                self.flowsheet().config.time, doc="Average Reynolds Number expression"
            )
            def N_Re_avg(b, t):
                return sum(b.N_Re[t, x] for x in self.length_domain) / self.nfe

        if hasattr(self, "K"):

            @self.Expression(
                self.flowsheet().config.time,
                self.config.property_package.solute_set,
                doc="Average mass transfer coefficient expression",
            )
            def K_avg(b, t, j):
                return sum(b.K[t, x, j] for x in self.difference_elements) / self.nfe

    ## should be called by add concentration polarization
    def _add_calculated_mass_transfer_coefficient(self):
        self._add_calculated_pressure_change_mass_transfer_components()

        solute_set = self.config.property_package.solute_set

        self.N_Sc_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=5e2,
            bounds=(1e2, 2e3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Schmidt number in membrane channel",
        )
        self.N_Sh_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            solute_set,
            initialize=1e2,
            bounds=(1, 3e2),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Sherwood number in membrane channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Mass transfer coefficient in membrane channel",
        )
        def eq_K(b, t, x, j):
            if b._skip_element(x):
                return Constraint.Skip
            return (
                b.K[t, x, j] * b.dh
                # TODO: add diff coefficient to SW prop and consider multi-components
                == b.properties[t, x].diffus_phase_comp["Liq", j] * b.N_Sh_comp[t, x, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Sherwood number",
        )
        def eq_N_Sh_comp(b, t, x, j):
            return (
                b.N_Sh_comp[t, x, j]
                == 0.46 * (b.N_Re[t, x] * b.N_Sc_comp[t, x, j]) ** 0.36
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.solute_set,
            doc="Schmidt number",
        )
        def eq_N_Sc_comp(b, t, x, j):
            return (
                b.N_Sc_comp[t, x, j]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.properties[t, x].diffus_phase_comp["Liq", j]
                == b.properties[t, x].visc_d_phase["Liq"]
            )

        return self.eq_K

    def _add_calculated_pressure_change_mass_transfer_components(self):
        # NOTE: This function could be called by either
        # `_add_calculated_pressure_change` *and/or*
        # `_add_calculated_mass_transfer_coefficient`.
        # Therefore, we add this simple guard against it being called twice.
        if hasattr(self, "channel_height"):
            return

        if not hasattr(self, "width"):
            raise ConfigurationError(
                f"Due to either a calculated mass transfer coefficient or a calculated pressure change, a ``width`` variable needs to be supplied to `add_geometry` for this MembraneChannel"
            )

        units_meta = self.config.property_package.get_metadata().get_derived_units

        if not hasattr(self, "area"):
            # comes from ControlVolume1D for 1DMC
            self.area = Var(
                initialize=1e-3 * 1 * 0.95,
                bounds=(0, 1e3),
                domain=NonNegativeReals,
                units=units_meta("length") ** 2,
                doc="Cross sectional area",
            )

        self.channel_height = Var(
            initialize=1e-3,
            bounds=(1e-4, 5e-3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="membrane-channel height",
        )

        self.dh = Var(
            initialize=1e-3,
            bounds=(1e-4, 5e-3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Hydraulic diameter of membrane channel",
        )

        self.spacer_porosity = Var(
            initialize=0.95,
            bounds=(0.1, 0.99),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="membrane-channel spacer porosity",
        )

        self.N_Re = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=5e2,
            bounds=(10, 5e3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Reynolds number in membrane channel",
        )

        @self.Constraint(doc="Hydraulic diameter")  # eqn. 17 in Schock & Miquel, 1987
        def eq_dh(b):
            return b.dh == 4 * b.spacer_porosity / (
                2 / b.channel_height + (1 - b.spacer_porosity) * 8 / b.channel_height
            )

        @self.Constraint(doc="Cross-sectional area")
        def eq_area(b):
            return b.area == b.channel_height * b.width * b.spacer_porosity

        @self.Constraint(
            self.flowsheet().config.time, self.length_domain, doc="Reynolds number"
        )
        def eq_N_Re(b, t, x):
            return (
                b.N_Re[t, x] * b.area * b.properties[t, x].visc_d_phase["Liq"]
                == sum(
                    b.properties[t, x].flow_mass_phase_comp["Liq", j]
                    for j in b.config.property_package.component_list
                )
                * b.dh
            )

    def _add_interface_stateblock(self, has_phase_equilibrium=None):
        """
        This method constructs the interface state blocks for the
        control volume.

        Args:
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
        Returns:
            None
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets

        self.properties_interface = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of feed-side membrane interface",
            **tmp_dict,
        )

    def _add_calculated_pressure_change(
        self, friction_factor=FrictionFactor.flat_sheet
    ):
        self._add_calculated_pressure_change_mass_transfer_components()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.velocity = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=0.5,
            bounds=(1e-2, 5),
            domain=NonNegativeReals,
            units=units_meta("length") / units_meta("time"),
            doc="Crossflow velocity in feed channel",
        )
        self.friction_factor_darcy = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=0.5,
            bounds=(1e-2, 5),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Darcy friction factor in feed channel",
        )

        # Constraints active when MassTransferCoefficient.calculated
        # Mass transfer coefficient calculation

        ## ==========================================================================
        # Crossflow velocity
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Crossflow velocity constraint",
        )
        def eq_velocity(b, t, x):
            return b.velocity[t, x] * b.area == b.properties[t, x].flow_vol_phase["Liq"]

        ## ==========================================================================
        # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
        if friction_factor == FrictionFactor.flat_sheet:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Darcy friction factor constraint for flat sheet membranes",
            )
            def eq_friction_factor(b, t, x):
                return (b.friction_factor_darcy[t, x] - 0.42) * b.N_Re[t, x] == 189.3

        # Darcy friction factor based on eq. 24 in Mass transfer and pressure loss in spiral wound modules (Schock & Miquel, 1987)
        elif friction_factor == FrictionFactor.spiral_wound:

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Darcy friction factor constraint for spiral-wound membranes",
            )
            def eq_friction_factor(b, t, x):
                return b.friction_factor_darcy[t, x] == 6.23 * b.N_Re[t, x] ** -0.3

        else:
            raise ConfigurationError(
                f"Unrecognized friction_factor type {friction_factor}"
            )

        ## ==========================================================================
        # Pressure change per unit length due to friction
        # -1/2*f/dh*density*velocity^2
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="pressure change per unit length due to friction",
        )
        def eq_dP_dx(b, t, x):
            return (
                b.dP_dx[t, x] * b.dh
                == -0.5
                * b.friction_factor_darcy[t, x]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.velocity[t, x] ** 2
            )

    def _get_state_args(self, initialize_guess, state_args):
        """
        Arguments:
            initialize_guess : a dict of guesses for solvent_recovery, solute_recovery,
                               and cp_modulus. These guesses offset the initial values
                               for the retentate, permeate, and membrane interface
                               state blocks from the inlet feed
                               (default =
                               {'deltaP': -1e4,
                               'solvent_recovery': 0.5,
                               'solute_recovery': 0.01,
                               'cp_modulus': 1.1})
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for the inlet
                         feed side state block (see documentation of the specific
                         property package).
        """

        # assumptions
        if initialize_guess is None:
            initialize_guess = {}
        if "deltaP" not in initialize_guess:
            initialize_guess["deltaP"] = -1e4
        if "solvent_recovery" not in initialize_guess:
            initialize_guess["solvent_recovery"] = 0.5
        if "solute_recovery" not in initialize_guess:
            initialize_guess["solute_recovery"] = 0.01

        # Get source block
        # TODO: need to re-visit for counterflow
        if self._flow_direction == FlowDirection.forward:
            source_idx = self.length_domain.first()
        else:
            source_idx = self.length_domain.last()
        source = self.properties[self.flowsheet().config.time.first(), source_idx]

        if state_args is None:
            state_args = {}
            state_dict = source.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        if "flow_mass_phase_comp" not in state_args.keys():
            raise ConfigurationError(
                f"{self.__class__.__name__} initialization routine expects "
                "flow_mass_phase_comp as a state variable. Check "
                "that the property package supports this state "
                "variable or that the state_args provided to the "
                "initialize call includes this state variable"
            )

        # slightly modify initial values for other state blocks
        state_args_retentate = deepcopy(state_args)

        state_args_retentate["pressure"] += initialize_guess["deltaP"]
        for j in self.config.property_package.solvent_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solvent_recovery"]
            )
        for j in self.config.property_package.solute_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solute_recovery"]
            )

        # slightly modify initial values for other state blocks
        state_args_permeate = deepcopy(state_args)

        state_args_permeate["pressure"] = 101325  # 1 bar
        for j in self.config.property_package.solvent_set:
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solvent_recovery"
            ]
        for j in self.config.property_package.solute_set:
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solute_recovery"
            ]

        return {
            "feed_side": state_args,
            "retentate": state_args_retentate,
            "permeate": state_args_permeate,
        }

    def _get_state_args_interface(self, initialize_guess, prop_in, prop_out):
        if initialize_guess is None:
            initialize_guess = {}
        if "cp_modulus" not in initialize_guess:
            if hasattr(self, "cp_modulus"):
                if self._flow_direction == FlowDirection.forward:
                    initialize_guess["cp_modulus"] = 1.1
                else:
                    initialize_guess["cp_modulus"] = 0.9
            else:
                initialize_guess["cp_modulus"] = 1

        state_args_interface_in = deepcopy(prop_in)
        state_args_interface_out = deepcopy(prop_out)

        for j in self.config.property_package.solute_set:
            state_args_interface_in["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]
            state_args_interface_out["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]

        x = 0.5
        state_args_tx = {}
        for k in state_args_interface_in:
            if isinstance(state_args_interface_in[k], dict):
                if k not in state_args_tx:
                    state_args_tx[k] = {}
                for index in state_args_interface_in[k]:
                    state_args_tx[k][index] = (1.0 - x) * state_args_interface_in[k][
                        index
                    ] + x * state_args_interface_out[k][index]
            else:
                state_args_tx[k] = (1.0 - x) * state_args_interface_in[
                    k
                ] + x * state_args_interface_out[k]
        state_args_interface = state_args_tx

        return state_args_interface

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "channel_height"):
            if iscale.get_scaling_factor(self.channel_height) is None:
                iscale.set_scaling_factor(self.channel_height, 1e3)

        if hasattr(self, "spacer_porosity"):
            if iscale.get_scaling_factor(self.spacer_porosity) is None:
                iscale.set_scaling_factor(self.spacer_porosity, 1)

        if hasattr(self, "dh"):
            if iscale.get_scaling_factor(self.dh) is None:
                iscale.set_scaling_factor(self.dh, 1e3)

        if hasattr(self, "K"):
            for v in self.K.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e4)

        if hasattr(self, "N_Re"):
            for t, x in self.N_Re.keys():
                if iscale.get_scaling_factor(self.N_Re[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Re[t, x], 1e-2)

        if hasattr(self, "N_Sc_comp"):
            for v in self.N_Sc_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, "N_Sh_comp"):
            for v in self.N_Sh_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, "velocity"):
            for v in self.velocity.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "friction_factor_darcy"):
            for v in self.friction_factor_darcy.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "cp_modulus"):
            for v in self.cp_modulus.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)


# helper for validating configuration arguments for this CV
def validate_membrane_config_args(unit):

    if (
        unit.config.pressure_change_type is not PressureChangeType.fixed_per_stage
        and unit.config.has_pressure_change is False
    ):
        raise ConfigurationError(
            "\nConflict between configuration options:\n"
            "'has_pressure_change' cannot be False "
            "while 'pressure_change_type' is set to {}.\n\n"
            "'pressure_change_type' must be set to PressureChangeType.fixed_per_stage\nor "
            "'has_pressure_change' must be set to True".format(
                unit.config.pressure_change_type
            )
        )

    if (
        unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
        and unit.config.mass_transfer_coefficient == MassTransferCoefficient.none
    ):
        raise ConfigurationError(
            "\n'mass_transfer_coefficient' and 'concentration_polarization_type' options configured incorrectly:\n"
            "'mass_transfer_coefficient' cannot be set to MassTransferCoefficient.none "
            "while 'concentration_polarization_type' is set to ConcentrationPolarizationType.calculated.\n "
            "\n\nSet 'mass_transfer_coefficient' to MassTransferCoefficient.fixed or "
            "MassTransferCoefficient.calculated "
            "\nor set 'concentration_polarization_type' to ConcentrationPolarizationType.fixed or "
            "ConcentrationPolarizationType.none"
        )

    if (
        unit.config.concentration_polarization_type
        != ConcentrationPolarizationType.calculated
        and unit.config.mass_transfer_coefficient != MassTransferCoefficient.none
    ):
        raise ConfigurationError(
            "\nConflict between configuration options:\n"
            "'mass_transfer_coefficient' cannot be set to {} "
            "while 'concentration_polarization_type' is set to {}.\n\n"
            "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
            "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated".format(
                unit.config.mass_transfer_coefficient,
                unit.config.concentration_polarization_type,
            )
        )

    if (
        unit.config.pressure_change_type != PressureChangeType.calculated
        and unit.config.friction_factor != unit.config.get("friction_factor")._default
    ):
        raise ConfigurationError(
            "\nChanging the 'friction_factor' will have no effect if the 'pressure_change_type' is not `PressureChangeType.calculated`"
        )
