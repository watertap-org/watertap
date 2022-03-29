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

from copy import deepcopy
from enum import Enum, auto
from pyomo.environ import Block, exp, NonNegativeReals, Param, Suffix, Var, units as pyunits
from pyomo.common.collections import ComponentSet
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import UnitModelBlockData, useDefault, MaterialBalanceType,\
        EnergyBalanceType, MomentumBalanceType
from idaes.core.util import scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class ConcentrationPolarizationType(Enum):
    none = auto()                    # simplified assumption: no concentration polarization
    fixed = auto()                   # simplified assumption: concentration polarization modulus is a user specified value
    calculated = auto()              # calculate concentration polarization (concentration at membrane interface)


class MassTransferCoefficient(Enum):
    none = auto()                    # mass transfer coefficient not utilized for concentration polarization effect
    fixed = auto()                   # mass transfer coefficient is a user specified value
    calculated = auto()              # mass transfer coefficient is calculated
    # TODO: add option for users to define their own relationship?


class PressureChangeType(Enum):
    fixed_per_stage = auto()         # pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length = auto()   # pressure drop per unit length across membrane channel is a user-specified value
    calculated = auto()              # pressure drop across membrane channel is calculated


class _ReverseOsmosisBaseData(UnitModelBlockData):
    """
    Reverse Osmosis base class
    """

    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not.
    **default** = False. Membrane units do not yet support dynamic
    behavior."""))

    CONFIG.declare("has_holdup", ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Membrane units do not have defined volume, thus
    this must be False."""))

    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))

    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these.
    **default** - None.
    **Valid values:** {
    see property package for documentation.}"""))

    CONFIG.declare("material_balance_type", ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}"""))

    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed.
    **default** - EnergyBalanceType.useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))

    CONFIG.declare("momentum_balance_type", ConfigValue(
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
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))

    CONFIG.declare("concentration_polarization_type", ConfigValue(
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
        """))

    CONFIG.declare("mass_transfer_coefficient", ConfigValue(
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
        """))

    CONFIG.declare("has_pressure_change", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}"""))

    CONFIG.declare("pressure_change_type", ConfigValue(
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

"""))

    CONFIG.declare("has_full_reporting", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Level of reporting results",
            doc="""Level of reporting results.
            **default** - False.
            **Valid values:** {
            **False** - include minimal reporting of results,
            **True** - report additional properties of interest that aren't constructed by
            the unit model by default. Also, report averaged expression values"""))


    def build(self):
        """
        Common variables and constraints for an RO unit model
        """
        super().build()

        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError("Membrane models only support one solvent component,"
                                     "the provided property package has specified {} solvent components"
                                     .format(len(self.config.property_package.solvent_set)))

        if (len(self.config.property_package.phase_list) > 1
                or 'Liq' not in [p for p in self.config.property_package.phase_list]):
            raise ConfigurationError(
                "Membrane models only support one liquid phase ['Liq'],"
                "the property package has specified the following phases {}"
                .format([p for p in self.config.property_package.phase_list]))

        if (self.config.pressure_change_type is not PressureChangeType.fixed_per_stage
                and self.config.has_pressure_change is False):
            raise ConfigurationError(
                "\nConflict between configuration options:\n"
                "'has_pressure_change' cannot be False "
                "while 'pressure_change_type' is set to {}.\n\n"
                "'pressure_change_type' must be set to PressureChangeType.fixed_per_stage\nor "
                "'has_pressure_change' must be set to True"
                .format(self.config.pressure_change_type))

        if (self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated
                and self.config.mass_transfer_coefficient == MassTransferCoefficient.none):
            raise ConfigurationError(
                "\n'mass_transfer_coefficient' and 'concentration_polarization_type' options configured incorrectly:\n"
                "'mass_transfer_coefficient' cannot be set to MassTransferCoefficient.none "
                "while 'concentration_polarization_type' is set to ConcentrationPolarizationType.calculated.\n "
                "\n\nSet 'mass_transfer_coefficient' to MassTransferCoefficient.fixed or "
                "MassTransferCoefficient.calculated "
                "\nor set 'concentration_polarization_type' to ConcentrationPolarizationType.fixed or "
                "ConcentrationPolarizationType.none")

        if (self.config.concentration_polarization_type != ConcentrationPolarizationType.calculated
                and self.config.mass_transfer_coefficient != MassTransferCoefficient.none):
            raise ConfigurationError(
                "\nConflict between configuration options:\n"
                "'mass_transfer_coefficient' cannot be set to {} "
                "while 'concentration_polarization_type' is set to {}.\n\n"
                "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
                "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated"
                    .format(self.config.mass_transfer_coefficient, self.config.concentration_polarization_type))

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

    def _make_performance(self):

        # For permeate-specific scaling in calculate_scaling_factors
        self._permeate_scaled_properties = ComponentSet()

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        units_meta = \
            self.config.property_package.get_metadata().get_derived_units

        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements")

        """ Unit model variables"""
        self.A_comp = Var(
            self.flowsheet().config.time,
            solvent_set,
            initialize=1e-12,
            bounds=(1e-18, 1e-6),
            domain=NonNegativeReals,
            units=units_meta('length') * units_meta('pressure') ** -1 * units_meta('time') ** -1,
            doc='Solvent permeability coeff.')

        self.B_comp = Var(
            self.flowsheet().config.time,
            solute_set,
            initialize=1e-8,
            bounds=(1e-11, 1e-5),
            domain=NonNegativeReals,
            units=units_meta('length')*units_meta('time')**-1,
            doc='Solute permeability coeff.')

        # TODO: add water density to NaCl prop model and remove here (or use IDAES version)
        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta('mass')*units_meta('length')**-3,
            doc='Pure water density')

        self.area = Var(
            initialize=10,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta('length')**2,
            doc='Membrane area')

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b,t,p,j : 0.4037 if j in solvent_set else 0.0033,
            bounds=lambda b,t,p,j : (1e-2, 1 - 1e-6) if j in solvent_set else (1e-5, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Mass-based component recovery')

        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.4,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Volumetric recovery rate')

        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Observed solute rejection')

        self.flux_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b,t,x,p,j : 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b,t,x,p,j : (1e-4, 3e-2) if j in solvent_set else (1e-8, 1e-3),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Mass flux across membrane at inlet and outlet')

        if ((self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
                or self.config.pressure_change_type == PressureChangeType.calculated):

            self.channel_height = Var(
                initialize=1e-3,
                bounds=(1e-4, 5e-3),
                domain=NonNegativeReals,
                units=units_meta('length'),
                doc='Feed-channel height')

            self.dh = Var(
                initialize=1e-3,
                bounds=(1e-4, 5e-3),
                domain=NonNegativeReals,
                units=units_meta('length'),
                doc='Hydraulic diameter of feed channel')

            self.spacer_porosity = Var(
                initialize=0.95,
                bounds=(0.1, 0.99),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc='Feed-channel spacer porosity')

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated:
            self.Kf = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta('length') * units_meta('time')**-1,
                doc='Mass transfer coefficient in feed channel')
        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
                or self.config.pressure_change_type == PressureChangeType.calculated):
            self.N_Re = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=5e2,
                bounds=(10, 5e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Reynolds number in feed channel")
        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            self.N_Sc = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=5e2,
                bounds=(1e2, 2e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Schmidt number in feed channel")
            self.N_Sh = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=1e2,
                bounds=(1, 3e2),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Sherwood number in feed channel")

        if self.config.pressure_change_type == PressureChangeType.calculated:
            self.velocity = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=units_meta('length')/units_meta('time'),
                doc="Crossflow velocity in feed channel")
            self.friction_factor_darcy = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Darcy friction factor in feed channel")

        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            or self.config.pressure_change_type == PressureChangeType.calculated):

            @self.Constraint(doc="Hydraulic diameter")  # eqn. 17 in Schock & Miquel, 1987
            def eq_dh(b):
                return (b.dh ==
                        4 * b.spacer_porosity
                        / (2 / b.channel_height
                           + (1 - b.spacer_porosity) * 8 / b.channel_height))

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.fixed:
            self.cp_modulus = Var(
                self.flowsheet().config.time,
                self.length_domain,
                solute_set,
                initialize=1.1,
                bounds=(0.9, 3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc='Concentration polarization modulus')

        # ==========================================================================
        # Volumetric Recovery rate
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (b.recovery_vol_phase[t, 'Liq'] ==
                    b.mixed_permeate[t].flow_vol_phase['Liq'] /
                    b.feed_side.properties[t, self.first_element].flow_vol_phase['Liq'])

        # ==========================================================================
        # Mass-based Component Recovery rate
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list)
        def eq_recovery_mass_phase_comp(b, t, j):
            return (b.recovery_mass_phase_comp[t, 'Liq', j]
                    * b.feed_side.properties[t, b.first_element].flow_mass_phase_comp['Liq', j] ==
                    b.mixed_permeate[t].flow_mass_phase_comp['Liq', j])

        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Solvent and solute mass flux")
        def eq_flux_mass(b, t, x, p, j):
            prop_feed = b.feed_side.properties[t, x]
            prop_perm = b.permeate_side[t, x]
            interface = b.feed_side.properties_interface[t, x]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return (b.flux_mass_phase_comp[t, x, p, j] == b.A_comp[t, j] * b.dens_solvent
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - (interface.pressure_osm - prop_perm.pressure_osm)))
            elif comp.is_solute():
                return (b.flux_mass_phase_comp[t, x, p, j] == b.B_comp[t, j]
                        * (interface.conc_mass_phase_comp[p, j] - prop_perm.conc_mass_phase_comp[p, j]))

        # # ==========================================================================
        # Concentration polarization
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.difference_elements,
                                   solute_set,
                                   doc="Concentration polarization")
        def eq_concentration_polarization(b, t, x, j):
            bulk = b.properties[t, x]
            interface = b.properties_interface[t, x]
            if self.config.concentration_polarization_type == ConcentrationPolarizationType.none:
                return interface.conc_mass_phase_comp['Liq', j] == \
                       bulk.conc_mass_phase_comp['Liq', j]
            elif self.config.concentration_polarization_type == ConcentrationPolarizationType.fixed:
                return (interface.conc_mass_phase_comp['Liq', j] ==
                        bulk.conc_mass_phase_comp['Liq', j]
                        * self.cp_modulus[t, x, j])
            elif self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated:
                jw = self.flux_mass_phase_comp[t, x, 'Liq', 'H2O'] / self.dens_solvent
                js = self.flux_mass_phase_comp[t, x, 'Liq', j]
                return (interface.conc_mass_phase_comp['Liq', j] ==
                        bulk.conc_mass_phase_comp['Liq', j] * exp(jw / self.Kf[t, x, j])
                        - js / jw * (exp(jw / self.Kf[t, x, j]) - 1))

        # # ==========================================================================
        # Feed and permeate-side isothermal conditions
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t, x):
            return b.feed_side.properties[t, x].temperature == \
                   b.permeate_side[t, x].temperature

        # ==========================================================================
        # isothermal conditions at permeate outlet
        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate out")
        def eq_permeate_outlet_isothermal(b, t):
            return b.feed_side.properties[t, b.length_domain.last()].temperature == \
                   b.mixed_permeate[t].temperature

        # ==========================================================================
        # isobaric conditions across permeate channel and at permeate outlet
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Isobaric assumption for permeate out")
        def eq_permeate_outlet_isobaric(b, t, x):
            return b.permeate_side[t, x].pressure == \
                   b.mixed_permeate[t].pressure

        # rejection
        @self.Constraint(self.flowsheet().config.time,
                         solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return (b.rejection_phase_comp[t, 'Liq', j] ==
                    1 - (b.mixed_permeate[t].conc_mass_phase_comp['Liq', j] /
                         b.feed_side.properties[t, self.first_element].conc_mass_phase_comp['Liq', j]))

        # ==========================================================================
        # Membrane area equation
        if hasattr(self, 'length') or hasattr(self, 'width'):
            @self.Constraint(doc="Membrane area")
            def eq_area(b):
                return b.area == b.length * b.width

        # # ==========================================================================
        # Constraints active when MassTransferCoefficient.calculated
        # Mass transfer coefficient calculation
        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            @self.Constraint(self.flowsheet().config.time,
                             self.difference_elements,
                             solute_set,
                             doc="Mass transfer coefficient in feed channel")
            def eq_Kf(b, t, x, j):
                bulk = b.feed_side.properties[t, x]
                return (b.Kf[t, x, j] * b.dh ==
                        bulk.diffus_phase['Liq']  # TODO: add diff coefficient to SW prop and consider multi-components
                        * b.N_Sh[t, x])

            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Sherwood number")
            def eq_N_Sh(b, t, x):
                return (b.N_Sh[t, x] ==
                        0.46 * (b.N_Re[t, x] * b.N_Sc[t, x])**0.36)

            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Schmidt number")
            def eq_N_Sc(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.N_Sc[t, x] * bulk.dens_mass_phase['Liq'] * bulk.diffus_phase['Liq'] ==
                        bulk.visc_d_phase['Liq'])

        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            or self.config.pressure_change_type == PressureChangeType.calculated):
            @self.Constraint(doc="Cross-sectional area")
            def eq_area_cross(b):
                return b.area_cross == b.channel_height * b.width * b.spacer_porosity

            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Reynolds number")
            def eq_N_Re(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.N_Re[t, x] * b.area_cross * bulk.visc_d_phase['Liq'] ==
                        sum(bulk.flow_mass_phase_comp['Liq', j] for j in b.config.property_package.component_list)
                        * b.dh)

        if self.config.pressure_change_type == PressureChangeType.calculated:
            ## ==========================================================================
            # Crossflow velocity
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Crossflow velocity constraint")
            def eq_velocity(b, t, x):
                return b.velocity[t, x] * b.area_cross == b.feed_side.properties[t, x].flow_vol_phase['Liq']

            ## ==========================================================================
            # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
            # TODO: this relationship for friction factor is specific to a particular spacer geometry. Add alternatives.
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Darcy friction factor constraint")
            def eq_friction_factor_darcy(b, t, x):
                return (b.friction_factor_darcy[t, x] - 0.42) * b.N_Re[t, x] == 189.3

            ## ==========================================================================
            # Pressure change per unit length due to friction
            # -1/2*f/dh*density*velocity^2
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="pressure change per unit length due to friction")
            def eq_dP_dx(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.dP_dx[t, x] * b.dh ==
                        -0.5 * b.friction_factor_darcy[t, x]
                        * b.feed_side.properties[t, x].dens_mass_phase['Liq'] * b.velocity[t, x]**2)

        # ==========================================================================
        # Bulk and interface connections on the feed-side
        # TEMPERATURE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.difference_elements,
                                   doc="Temperature at interface")
        def eq_equal_temp_interface(b, t, x):
            return b.properties_interface[t, x].temperature == \
                   b.properties[t, x].temperature
        # PRESSURE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.difference_elements,
                                   doc="Pressure at interface")
        def eq_equal_pressure_interface(b, t, x):
            return b.properties_interface[t, x].pressure == \
                   b.properties[t, x].pressure
        # VOLUMETRIC FLOWRATE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.difference_elements,
                                   doc="Volumetric flow at interface of inlet")
        def eq_equal_flow_vol_interface(b, t, x):
            return b.properties_interface[t, x].flow_vol_phase['Liq'] ==\
                   b.properties[t, x].flow_vol_phase['Liq']

        @self.Expression(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Average flux expression")
        def flux_mass_phase_comp_avg(b, t, p, j):
            return sum(b.flux_mass_phase_comp[t, x, p, j]
                       for x in self.difference_elements) / self.nfe

    def _add_expressions(self):
        """
        Generate expressions for additional results desired for full report
        """

        solute_set = self.config.property_package.solute_set
        if hasattr(self, 'N_Re'):
            @self.Expression(self.flowsheet().config.time,
                             doc="Average Reynolds Number expression")
            def N_Re_avg(b, t):
                return sum(b.N_Re[t, x]
                           for x in self.length_domain) / self.nfe
        if hasattr(self, 'Kf'):
            @self.Expression(self.flowsheet().config.time,
                             solute_set,
                             doc="Average mass transfer coefficient expression")
            def Kf_avg(b, t, j):
                return sum(b.Kf[t, x, j]
                           for x in self.difference_elements) / self.nfe

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )

    def get_costing(self, module=None, **kwargs):
        self.costing = Block()
        module.ReverseOsmosis_costing(self.costing, **kwargs)

    def _get_state_args(self, source, mixed_permeate_properties, initialize_guess, state_args): 
        '''
        Arguments:
            source : property model containing inlet feed
            mixed_permeate_properties : mixed permeate property block
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
        '''

        # assumptions
        if initialize_guess is None:
            initialize_guess = {}
        if 'deltaP' not in initialize_guess:
            initialize_guess['deltaP'] = -1e4
        if 'solvent_recovery' not in initialize_guess:
            initialize_guess['solvent_recovery'] = 0.5
        if 'solute_recovery' not in initialize_guess:
            initialize_guess['solute_recovery'] = 0.01
        if 'cp_modulus' not in initialize_guess:
            initialize_guess['cp_modulus'] = 1.1

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

        if 'flow_mass_phase_comp' not in state_args.keys():
            raise ConfigurationError(f'{self.__class__.__name__} initialization routine expects '
                                     'flow_mass_phase_comp as a state variable. Check '
                                     'that the property package supports this state '
                                     'variable or that the state_args provided to the '
                                     'initialize call includes this state variable')

        # slightly modify initial values for other state blocks
        state_args_retentate = deepcopy(state_args)
        state_args_permeate = deepcopy(state_args)

        state_args_retentate['pressure'] += initialize_guess['deltaP']
        state_args_permeate['pressure'] = mixed_permeate_properties.pressure.value
        for j in self.config.property_package.solvent_set:
            state_args_retentate['flow_mass_phase_comp'][('Liq', j)] *= (1 - initialize_guess['solvent_recovery'])
            state_args_permeate['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['solvent_recovery']
        for j in self.config.property_package.solute_set:
            state_args_retentate['flow_mass_phase_comp'][('Liq', j)] *= (1 - initialize_guess['solute_recovery'])
            state_args_permeate['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['solute_recovery']

        state_args_interface_in = deepcopy(state_args)
        state_args_interface_out = deepcopy(state_args_retentate)

        for j in self.config.property_package.solute_set:
            state_args_interface_in['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['cp_modulus']
            state_args_interface_out['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['cp_modulus']

        # TODO: I think this is what we'd like to do, but IDAES needs to be changed
        state_args_interface = {}
        for t in self.flowsheet().config.time:
            for x in self.length_domain:
                assert 0. <= x <= 1.
                state_args_tx = {}
                for k in state_args_interface_in:
                    if isinstance(state_args_interface_in[k], dict):
                        if k not in state_args_tx:
                            state_args_tx[k] = {}
                        for index in state_args_interface_in[k]:
                            state_args_tx[k][index] = (1.-x)*state_args_interface_in[k][index] + x*state_args_interface_out[k][index]
                    else:
                        state_args_tx[k] = (1.-x)*state_args_interface_in[k] + x*state_args_interface_out[k]
                state_args_interface[t,x] = state_args_tx

        x = 0.5
        state_args_tx = {}
        for k in state_args_interface_in:
            if isinstance(state_args_interface_in[k], dict):
                if k not in state_args_tx:
                    state_args_tx[k] = {}
                for index in state_args_interface_in[k]:
                    state_args_tx[k][index] = (1.-x)*state_args_interface_in[k][index] + x*state_args_interface_out[k][index]
            else:
                state_args_tx[k] = (1.-x)*state_args_interface_in[k] + x*state_args_interface_out[k]
        state_args_interface = state_args_tx

        return {'feed_side' : state_args,
                'retentate' : state_args_retentate,
                'permeate' : state_args_permeate,
                'interface' : state_args_interface,
               }

    def _get_performance_contents(self, time_point=0):
        x_in = self.first_element
        x_interface_in = self.difference_elements.first()
        x_out = self.length_domain.last()
        feed_inlet = self.feed_side.properties[time_point, x_in]
        feed_outlet = self.feed_side.properties[time_point, x_out]
        interface_inlet = self.feed_side.properties_interface[time_point, x_interface_in]
        interface_outlet = self.feed_side.properties_interface[time_point, x_out]
        permeate = self.mixed_permeate[time_point]
        var_dict = {}
        expr_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[time_point, 'Liq']
        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[time_point, 'Liq', 'H2O']
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "length") and self.config.has_full_reporting:
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width") and self.config.has_full_reporting:
            var_dict["Membrane Width"] = self.width
        if hasattr(self, "deltaP") and self.config.has_full_reporting:
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if hasattr(self, "N_Re") and self.config.has_full_reporting:
            var_dict["Reynolds Number @Inlet"] = self.N_Re[time_point, x_in]
            var_dict["Reynolds Number @Outlet"] = self.N_Re[time_point, x_out]
        if hasattr(self, "velocity") and self.config.has_full_reporting:
            var_dict["Velocity @Inlet"] = self.velocity[time_point, x_in]
            var_dict["Velocity @Outlet"] = self.velocity[time_point, x_out]
        for j in self.config.property_package.solute_set:
            if interface_inlet.is_property_constructed('conc_mass_phase_comp') and self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Inlet,Membrane-Interface '] = (
                    interface_inlet.conc_mass_phase_comp['Liq', j])
            if interface_outlet.is_property_constructed('conc_mass_phase_comp') and self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Outlet,Membrane-Interface '] = (
                    interface_outlet.conc_mass_phase_comp['Liq', j])
            if feed_inlet.is_property_constructed('conc_mass_phase_comp') and self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Inlet,Bulk'] = (
                    feed_inlet.conc_mass_phase_comp['Liq', j])
            if feed_outlet.is_property_constructed('conc_mass_phase_comp') and self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Outlet,Bulk'] = (
                    feed_outlet.conc_mass_phase_comp['Liq', j])
            if permeate.is_property_constructed('conc_mass_phase_comp') and self.config.has_full_reporting:
                var_dict[f'{j} Permeate Concentration'] = (
                    permeate.conc_mass_phase_comp['Liq', j])
        if interface_outlet.is_property_constructed('pressure_osm') and self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Outlet,Membrane-Interface '] = (
                interface_outlet.pressure_osm)
        if feed_outlet.is_property_constructed('pressure_osm') and self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Outlet,Bulk'] = (
                feed_outlet.pressure_osm)
        if interface_inlet.is_property_constructed('pressure_osm') and self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Inlet,Membrane-Interface'] = (
                interface_inlet.pressure_osm)
        if feed_inlet.is_property_constructed('pressure_osm') and self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Inlet,Bulk'] = (
                feed_inlet.pressure_osm)
        if feed_inlet.is_property_constructed('flow_vol_phase') and self.config.has_full_reporting:
            var_dict['Volumetric Flowrate @Inlet'] = (
                feed_inlet.flow_vol_phase['Liq'])
        if feed_outlet.is_property_constructed('flow_vol_phase') and self.config.has_full_reporting:
            var_dict['Volumetric Flowrate @Outlet'] = (
                feed_outlet.flow_vol_phase['Liq'])
        if hasattr(self, 'dh') and self.config.has_full_reporting:
            var_dict["Hydraulic Diameter"] = self.dh

        if self.config.has_full_reporting:
            expr_dict['Average Solvent Flux (LMH)'] = self.flux_mass_phase_comp_avg[time_point, 'Liq', 'H2O'] * 3.6e3
            expr_dict['Average Reynolds Number'] = self.N_Re_avg[time_point]
            for j in self.config.property_package.solute_set:
                expr_dict[f'{j} Average Solute Flux (GMH)'] = self.flux_mass_phase_comp_avg[time_point, 'Liq', j] * 3.6e6
                expr_dict[f'{j} Average Mass Transfer Coefficient (mm/h)'] = self.Kf_avg[time_point, j] * 3.6e6


        # TODO: add more vars
        return {"vars": var_dict, "exprs": expr_dict}

    # permeate properties need to rescale solute values by 100
    def _rescale_permeate_variable(self, var, factor=100):
        if var not in self._permeate_scaled_properties:
            sf = iscale.get_scaling_factor(var)
            iscale.set_scaling_factor(var, sf * factor)
            self._permeate_scaled_properties.add(var)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=10, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        if iscale.get_scaling_factor(self.A_comp) is None:
            iscale.set_scaling_factor(self.A_comp, 1e12)

        if iscale.get_scaling_factor(self.B_comp) is None:
            iscale.set_scaling_factor(self.B_comp, 1e8)

        if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
            iscale.set_scaling_factor(self.recovery_vol_phase, 1)

        for (t, p, j), v in self.recovery_mass_phase_comp.items():
            if j in self.config.property_package.solvent_set:
                sf = 1
            elif j in self.config.property_package.solute_set:
                sf = 100
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        for v in self.rejection_phase_comp.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'channel_height'):
            if iscale.get_scaling_factor(self.channel_height) is None:
                iscale.set_scaling_factor(self.channel_height, 1e3)

        if hasattr(self, 'spacer_porosity'):
            if iscale.get_scaling_factor(self.spacer_porosity) is None:
                iscale.set_scaling_factor(self.spacer_porosity, 1)

        if hasattr(self, 'dh'):
            if iscale.get_scaling_factor(self.dh) is None:
                iscale.set_scaling_factor(self.dh, 1e3)

        if hasattr(self, 'Kf'):
            for v in self.Kf.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e4)

        if hasattr(self, 'N_Re'):
            for t, x in self.N_Re.keys():
                if iscale.get_scaling_factor(self.N_Re[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Re[t, x], 1e-2)

        if hasattr(self, 'N_Sc'):
            for t, x in self.N_Sc.keys():
                if iscale.get_scaling_factor(self.N_Sc[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sc[t, x], 1e-2)

        if hasattr(self, 'N_Sh'):
            for t, x in self.N_Sh.keys():
                if iscale.get_scaling_factor(self.N_Sh[t, x]) is None:
                     iscale.set_scaling_factor(self.N_Sh[t, x], 1e-2)

        if hasattr(self, 'velocity'):
            for v in self.velocity.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'friction_factor_darcy'):
            for v in self.friction_factor_darcy.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'cp_modulus'):
            for v in self.cp_modulus.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        for sb in (self.permeate_side, self.mixed_permeate):
            for blk in sb.values():
                for j in self.config.property_package.solute_set:
                    self._rescale_permeate_variable(blk.flow_mass_phase_comp['Liq', j])
                    if blk.is_property_constructed('mass_frac_phase_comp'):
                        self._rescale_permeate_variable(blk.mass_frac_phase_comp['Liq', j])
                    if blk.is_property_constructed('conc_mass_phase_comp'):
                        self._rescale_permeate_variable(blk.conc_mass_phase_comp['Liq', j])
                    if blk.is_property_constructed('mole_frac_phase_comp'):
                        self._rescale_permeate_variable(blk.mole_frac_phase_comp[j])
                    if blk.is_property_constructed('molality_comp'):
                        self._rescale_permeate_variable(blk.molality_comp[j])
                if blk.is_property_constructed('pressure_osm'):
                    self._rescale_permeate_variable(blk.pressure_osm)

        for (t, x, p, j), v in self.flux_mass_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():  # scaling based on solvent flux equation
                    sf = (iscale.get_scaling_factor(self.A_comp[t, j])
                          * iscale.get_scaling_factor(self.dens_solvent)
                          * iscale.get_scaling_factor(self.feed_side.properties[t, x].pressure))
                    iscale.set_scaling_factor(v, sf)
                elif comp.is_solute():  # scaling based on solute flux equation
                    sf = (iscale.get_scaling_factor(self.B_comp[t, j])
                          * iscale.get_scaling_factor(self.feed_side.properties[t, x].conc_mass_phase_comp[p, j]))
                    iscale.set_scaling_factor(v, sf)
