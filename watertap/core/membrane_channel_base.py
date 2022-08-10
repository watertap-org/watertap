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

from pyomo.common.config import Bool, ConfigBlock, In

from idaes.core import (
    FlowDirection,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.core.util.exceptions import ConfigurationError


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


class PressureChangeType(Enum):
    """
    fixed_per_stage: pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length: pressure drop per unit length across membrane channel is a user-specified value
    calculated: pressure drop across membrane channel is calculated
    """
    fixed_per_stage = auto()
    fixed_per_unit_length = auto()
    calculated = auto()


CONFIG = ConfigBlock()

CONFIG.declare(
    "material_balance_type",
    ConfigValue(
        default=useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
    ),
)

CONFIG.declare(
    "energy_balance_type",
    ConfigValue(
        default=useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed.
**default** - useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
    ),
)

CONFIG.declare(
    "momentum_balance_type",
    ConfigValue(
        default=useDefault,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - useDefault.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
    ),
)

CONFIG.declare(
    "concentration_polarization_type",
    ConfigValue(
        default=useDefault,
        domain=In(ConcentrationPolarizationType),
        description="External concentration polarization effect in RO",
        doc="""
        Options to account for concentration polarization.

        **default** - useDefault

    .. csv-table::
        :header: "Configuration Options", "Description"

        "``ConcentrationPolarizationType.none``", "Simplifying assumption to ignore concentration polarization"
        "``ConcentrationPolarizationType.fixed``", "Specify an estimated value for the concentration polarization modulus"
        "``ConcentrationPolarizationType.calculated``", "Allow model to perform calculation of membrane-interface concentration"
    """,
    ),
)

CONFIG.declare(
    "mass_transfer_coefficient",
    ConfigValue(
        default=useDefault,
        domain=In(MassTransferCoefficient),
        description="Mass transfer coefficient in RO feed channel",
        doc="""
        Options to account for mass transfer coefficient.

        **default** - useDefault

    .. csv-table::
        :header: "Configuration Options", "Description"

        "``MassTransferCoefficient.none``", "Mass transfer coefficient not used in calculations"
        "``MassTransferCoefficient.fixed``", "Specify an estimated value for the mass transfer coefficient in the feed channel"
        "``MassTransferCoefficient.calculated``", "Allow model to perform calculation of mass transfer coefficient"
    """,
    ),
)

CONFIG.declare(
    "has_pressure_change",
    ConfigValue(
        default=useDefault,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - useDefault.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
    ),
)

   CONFIG.declare(
       "pressure_change_type",
       ConfigValue(
           default=useDefault,
           domain=In(PressureChangeType),
           description="Pressure change term construction flag",
           doc="""
       Indicates what type of pressure change calculation will be made. To use any of the
       ``pressure_change_type`` options to account for pressure drop, the configuration keyword
       ``has_pressure_change`` must also be set to ``True``. Also, if a value is specified for pressure
       change, it should be negative to represent pressure drop.

       **default** - useDefault


   .. csv-table::
       :header: "Configuration Options", "Description"

       "``PressureChangeType.fixed_per_stage``", "Specify an estimated value for pressure drop across the membrane feed channel"
       "``PressureChangeType.fixed_per_unit_length``", "Specify an estimated value for pressure drop per unit length across the membrane feed channel"
       "``PressureChangeType.calculated``", "Allow model to perform calculation of pressure drop across the membrane feed channel"

""",
    ),
)

class MembraneChannelMixin:

    def add_total_pressure_balances(
        self,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        custom_term=None
    ):
        super().add_total_pressure_balances(
            has_pressure_change=has_pressure_change,
            custom_term=custom_term)

        self._add_membrane_pressure_balances()

        if pressure_change_type == PressureChangeType.calculated:
            self._add_calculated_pressure_change()

    def add_flux_balance(self):

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        self.A_comp = Var(
            self.flowsheet().config.time,
            solvent_set,
            initialize=1e-12,
            bounds=(1e-18, 1e-6),
            domain=NonNegativeReals,
            units=units_meta("length")
            * units_meta("pressure") ** -1
            * units_meta("time") ** -1,
            doc="Solvent permeability coeff.",
        )

        self.B_comp = Var(
            self.flowsheet().config.time,
            solute_set,
            initialize=1e-8,
            bounds=(1e-11, 1e-5),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Solute permeability coeff.",
        )

        # TODO: add water density to NaCl prop model and remove here (or use IDAES version)
        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        self.flux_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, x, p, j: 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b, t, x, p, j: (1e-4, 3e-2)
            if j in solvent_set
            else (1e-8, 1e-3),
            units=units_meta("mass")
            * units_meta("length") ** -2
            * units_meta("time") ** -1,
            doc="Mass flux across membrane at inlet and outlet",
        )


        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Solvent and solute mass flux",
        )
        def eq_flux_mass(b, t, x, p, j):
            prop_feed = b.properties[t, x]
            prop_perm = b.permeate_side[t, x]
            interface = b.properties_interface[t, x]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return b.flux_mass_phase_comp[t, x, p, j] == b.A_comp[
                    t, j
                ] * b.dens_solvent * (
                    (prop_feed.pressure - prop_perm.pressure)
                    - (
                        interface.pressure_osm_phase[p]
                        - prop_perm.pressure_osm_phase[p]
                    )
                )
            elif comp.is_solute():
                return b.flux_mass_phase_comp[t, x, p, j] == b.B_comp[t, j] * (
                    interface.conc_mass_phase_comp[p, j]
                    - prop_perm.conc_mass_phase_comp[p, j]
                )

        return self.eq_flux_mass

    def add_concentration_polarization(
        self,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    ):
        if concentration_polarization_type == ConcentrationPolarizationType.none:
            @self.Constraint(
               self.flowsheet().config.time,
               self.difference_elements,
               solute_set,
               doc="Concentration polarization",
            )   
            def eq_concentration_polarization(b, t, x, j):
                return (
                    b.properties_interface[t,x].conc_mass_phase_comp["Liq", j]
                    == b.properties[t,x].conc_mass_phase_comp["Liq", j]
                )

        elif concentration_polarization_type == ConcentrationPolarizationType.fixed:

            self.cp_modulus = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=1.1,
                bounds=(0.9, 3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Concentration polarization modulus",
            )

            @self.Constraint(
               self.flowsheet().config.time,
               self.difference_elements,
               solute_set,
               doc="Concentration polarization",
            )   
            def eq_concentration_polarization(b, t, x, j):
                return (
                    self.properties_interface[t,x].conc_mass_phase_comp["Liq", j]
                    == self.properties[t,x].conc_mass_phase_comp["Liq", j] * self.cp_modulus[t, x, j]
                )

        elif concentration_polarization_type == ConcentrationPolarizationType.calculated:
            if mass_transfer_coefficient == MassTransferCoefficient.none:
                raise ConfigurationError()

            self.K = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta("length") * units_meta("time") ** -1,
                doc="Mass transfer coefficient in membrane channel",
            )

            @self.Constraint(
               self.flowsheet().config.time,
               self.difference_elements,
               solute_set,
               doc="Concentration polarization",
            )   
            def eq_concentration_polarization(b, t, x, j):
                jw = self.flux_mass_phase_comp[t, x, "Liq", "H2O"] / self.dens_solvent
                js = self.flux_mass_phase_comp[t, x, "Liq", j]
                return self.properties_interface[t,x].conc_mass_phase_comp[
                    "Liq", j
                ] == self.properties[t,x].conc_mass_phase_comp["Liq", j] * exp(
                    jw / self.K[t, x, j]
                ) - js / jw * (
                    exp(jw / self.K[t, x, j]) - 1
                )

            if mass_transfer_coefficient == MassTransferCoefficient.calculated:
                self._add_calculated_mass_transfer_coefficient()

        else:
            raise ConfigurationError()
            
        return self.eq_concentration_polarization


    ## should be called by add concentration polarization
    def _add_calculated_mass_transfer_coefficient(self):
        self._add_calculated_pressure_change_mass_transfer_components()

        self.N_Sc = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=5e2,
            bounds=(1e2, 2e3),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Schmidt number in membrane channel",
        )
        self.N_Sh = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=1e2,
            bounds=(1, 3e2),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Sherwood number in membrane channel",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            solute_set,
            doc="Mass transfer coefficient in membrane channel",
        )
        def eq_K(b, t, x, j):
            return (
                b.K[t, x, j] * b.dh
                # TODO: add diff coefficient to SW prop and consider multi-components
                == b.properties[t, x].diffus_phase_comp["Liq", j] * b.N_Sh[t, x]
            )

        @self.Constraint(
            self.flowsheet().config.time, self.length_domain, doc="Sherwood number"
        )
        def eq_N_Sh(b, t, x):
            return b.N_Sh[t, x] == 0.46 * (b.N_Re[t, x] * b.N_Sc[t, x]) ** 0.36

        @self.Constraint(
            self.flowsheet().config.time, self.length_domain, doc="Schmidt number"
        )
        def eq_N_Sc(b, t, x):
            bulk = b.feed_side.properties[t, x]
            # # TODO: This needs to be revisted. Diffusion is now by component, but
            #   not H2O and this var should also be by component, but the implementation
            #   is not immediately clear.
            return (
                b.N_Sc[t, x]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.properties[t, x].diffus_phase_comp["Liq", b.properties[t, x].params.component_list.last()]
                == b.properties[t, x].visc_d_phase["Liq"]
            )

        return self.eq_K

    def _add_calculated_pressure_change_mass_transfer_components(self):
        if hasattr(self, "channel_height"):
            return

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

        @self.Constraint(
            doc="Hydraulic diameter"
        )  # eqn. 17 in Schock & Miquel, 1987
        def eq_dh(b):
            return b.dh == 4 * b.spacer_porosity / (
                2 / b.channel_height
                + (1 - b.spacer_porosity) * 8 / b.channel_height
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

    def _add_interface_blocks(
        self, information_flow=FlowDirection.forward, has_phase_equilibrium=None
    ):
        """
        This method constructs the interface state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments: dict-like object of arguments to be passed to
                                state blocks as construction arguments
        Returns:
            None
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets

        self.permeate_side = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict,
        )
        self.mixed_permeate = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict,
        )

        self.properties_interface = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of feed-side membrane interface",
            default=tmp_dict,
        )

    def _add_membrane_pressure_balances(self):

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isobaric assumption for permeate out",
        )
        def eq_permeate_outlet_isobaric(b, t, x):
            return b.permeate_side[t, x].pressure == b.mixed_permeate[t].pressure

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Pressure at interface",
        )
        def eq_equal_pressure_interface(b, t, x):
            return b.properties_interface[t, x].pressure == b.properties[t, x].pressure

    def _add_calculated_pressure_change(self):
        self._add_calculated_pressure_change_mass_transfer_components()

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
            return (
                b.velocity[t, x] * b.area
                == b.properties[t, x].flow_vol_phase["Liq"]
            )

        ## ==========================================================================
        # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
        # TODO: this relationship for friction factor is specific to a particular spacer geometry. Add alternatives.
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Darcy friction factor constraint",
        )
        def eq_friction_factor_darcy(b, t, x):
            return (b.friction_factor_darcy[t, x] - 0.42) * b.N_Re[t, x] == 189.3

        ## ==========================================================================
        # Pressure change per unit length due to friction
        # -1/2*f/dh*density*velocity^2
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="pressure change per unit length due to friction",
        )
        def eq_deltaP(b, t, x):
            return (
                b.deltaP[t, x] * b.dh
                == -0.5
                * b.friction_factor_darcy[t, x]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.velocity[t, x] ** 2
            )


    def _validate_membrane_config_args(self):

        if (
            self.config.pressure_change_type is not PressureChangeType.fixed_per_stage
            and self.config.has_pressure_change is False
        ):
            raise ConfigurationError(
                "\nConflict between configuration options:\n"
                "'has_pressure_change' cannot be False "
                "while 'pressure_change_type' is set to {}.\n\n"
                "'pressure_change_type' must be set to PressureChangeType.fixed_per_stage\nor "
                "'has_pressure_change' must be set to True".format(
                    self.config.pressure_change_type
                )
            )

        if (
            self.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
            and self.config.mass_transfer_coefficient == MassTransferCoefficient.none
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
            self.config.concentration_polarization_type
            != ConcentrationPolarizationType.calculated
            and self.config.mass_transfer_coefficient != MassTransferCoefficient.none
        ):
            raise ConfigurationError(
                "\nConflict between configuration options:\n"
                "'mass_transfer_coefficient' cannot be set to {} "
                "while 'concentration_polarization_type' is set to {}.\n\n"
                "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
                "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated".format(
                    self.config.mass_transfer_coefficient,
                    self.config.concentration_polarization_type,
                )
            )
