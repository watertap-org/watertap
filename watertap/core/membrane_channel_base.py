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

from pyomo.common.collections import ComponentSet
from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import (
    Param,
    NegativeReals,
    NonNegativeReals,
    Var,
    check_optimal_termination,
    exp,
    units as pyunits,
)

from idaes.core import (
    FlowDirection,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError

import idaes.logger as idaeslog


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


class PressureChangeType(Enum):
    """
    fixed_per_stage: pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length: pressure drop per unit length across membrane channel is a user-specified value
    calculated: pressure drop across membrane channel is calculated
    """
    fixed_per_stage = auto()
    fixed_per_unit_length = auto()
    calculated = auto()


CONFIG_Template = ConfigDict()

CONFIG_Template.declare(
    "dynamic",
    ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not.
**default** = False. Membrane units do not yet support dynamic
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
**MaterialBalanceType.useDefault - refer to property package for default
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
**EnergyBalanceType.useDefault - refer to property package for default
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

class MembraneChannelMixin:

    def add_mass_transfer(self):
        raise NotImplementedError()

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        raise NotImplementedError()

    def _add_pressure_change_equation(self):
        @self.Constraint(
            self.flowsheet().config.time, doc="Total Pressure drop across channel"
        )
        def eq_pressure_change(b, t):
            return b.pressure_change_total[t] == sum(
                b.dP_dx[t, x] * b.length / b.nfe for x in b.difference_elements
            )

    def _add_recovery_rejection(self,**kwrags):

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, p, j: 0.4037 if j in solvent_set else 0.0033,
            bounds=lambda b, t, p, j: (1e-2, 1 - 1e-6)
            if j in solvent_set
            else (1e-5, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Mass-based component recovery",
        )

        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Observed solute rejection",
        )

        # ==========================================================================
        # Mass-based Component Recovery rate
        @self.Constraint(
            self.flowsheet().config.time, self.config.property_package.component_list
        )
        def eq_recovery_mass_phase_comp(b, t, j):
            return (
                b.recovery_mass_phase_comp[t, "Liq", j]
                * b.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", j
                ]
                == b.mixed_permeate[t].flow_mass_phase_comp["Liq", j]
            )

        # rejection
        @self.Constraint(self.flowsheet().config.time, solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return b.rejection_phase_comp[t, "Liq", j] == 1 - (
                b.mixed_permeate[t].conc_mass_phase_comp["Liq", j]
                / b.properties[t, self.first_element].conc_mass_phase_comp[
                    "Liq", j
                ]
            )

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

        if has_pressure_change:
            self._add_pressure_change(pressure_change_type=pressure_change_type)

        if pressure_change_type == PressureChangeType.calculated:
            self._add_calculated_pressure_change()

    def add_isothermal_conditions(self):

        # # ==========================================================================
        # Feed and permeate-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isothermal assumption for permeate",
        )
        def eq_permeate_isothermal(b, t, x):
            return (
                b.properties[t, x].temperature
                == b.permeate_side[t, x].temperature
            )

        # ==========================================================================
        # Bulk and interface connections on the feed-side
        # TEMPERATURE
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Temperature at interface",
        )
        def eq_equal_temp_interface(b, t, x):
            return (
                b.properties[t, x].temperature
                == b.properties_interface[t, x].temperature
            )

        # ==========================================================================
        # isothermal conditions at permeate outlet
        @self.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for permeate out"
        )
        def eq_permeate_outlet_isothermal(b, t):
            return (
                b.properties[t, b.length_domain.last()].temperature
                == b.mixed_permeate[t].temperature
            )


    def add_volumetric_flowrate_balance(self):
        # VOLUMETRIC FLOWRATE
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Volumetric flow at interface of inlet",
        )
        def eq_equal_flow_vol_interface(b, t, x):
            return (
                b.properties_interface[t, x].flow_vol_phase["Liq"]
                == b.properties[t, x].flow_vol_phase["Liq"]
            )


    def add_flux_balance(self):

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        units_meta = self.config.property_package.get_metadata().get_derived_units

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

        @self.Expression(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Average flux expression",
        )
        def flux_mass_phase_comp_avg(b, t, p, j):
            return (
                sum(
                    b.flux_mass_phase_comp[t, x, p, j] for x in self.difference_elements
                )
                / self.nfe
            )

        return self.eq_flux_mass

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

    def add_recovery_vol_phase(self):

        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.4,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Volumetric recovery rate",
        )

        # ==========================================================================
        # Volumetric Recovery rate
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (
                b.recovery_vol_phase[t, "Liq"]
                == b.mixed_permeate[t].flow_vol_phase["Liq"]
                / b.properties[t, self.first_element].flow_vol_phase["Liq"]
            )

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

    def _add_area_total(self):
        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.area_total = Var(
            initialize=10,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Total Membrane Channel area",
        )

    def _add_area_total_equation(self):
        # ==========================================================================
        # Membrane area equation
        @self.Constraint(doc="Total Membrane area")
        def eq_area_total(b):
            return b.area_total == b.length * b.width

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
            self.config.property_package.solute_set,
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
            bulk = b.properties[t, x]
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
        # NOTE: This function could be called by either
        # `_add_calculated_pressure_change` *and/or* 
        # `_add_calculated_mass_transfer_coefficient`.
        # Therefore, we add this simple gaurd against it being called twice.
        if hasattr(self, "channel_height"):
            return

        units_meta = self.config.property_package.get_metadata().get_derived_units

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
        self, has_phase_equilibrium=None
    ):
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
        def eq_dP_dx(b, t, x):
            return (
                b.dP_dx[t, x] * b.dh
                == -0.5
                * b.friction_factor_darcy[t, x]
                * b.properties[t, x].dens_mass_phase["Liq"]
                * b.velocity[t, x] ** 2
            )


    def _get_state_args(
        self, initialize_guess, state_args
    ):
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

        source = self.properties[
            self.flowsheet().config.time.first(), self.first_element
        ]
        mixed_permeate_properties = self.mixed_permeate[self.flowsheet().config.time.first()]

        # assumptions
        if initialize_guess is None:
            initialize_guess = {}
        if "deltaP" not in initialize_guess:
            initialize_guess["deltaP"] = -1e4
        if "solvent_recovery" not in initialize_guess:
            initialize_guess["solvent_recovery"] = 0.5
        if "solute_recovery" not in initialize_guess:
            initialize_guess["solute_recovery"] = 0.01
        if "cp_modulus" not in initialize_guess:
            if (
                self.config.concentration_polarization_type
                != ConcentrationPolarizationType.none
            ):
                initialize_guess["cp_modulus"] = 1.1
            else:
                initialize_guess["cp_modulus"] = 1

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
        state_args_permeate = deepcopy(state_args)

        state_args_retentate["pressure"] += initialize_guess["deltaP"]
        state_args_permeate["pressure"] = mixed_permeate_properties.pressure.value
        for j in self.config.property_package.solvent_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solvent_recovery"]
            )
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solvent_recovery"
            ]
        for j in self.config.property_package.solute_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solute_recovery"]
            )
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solute_recovery"
            ]

        state_args_interface_in = deepcopy(state_args)
        state_args_interface_out = deepcopy(state_args_retentate)

        for j in self.config.property_package.solute_set:
            state_args_interface_in["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]
            state_args_interface_out["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]

        # TODO: I think this is what we'd like to do, but IDAES needs to be changed
        state_args_interface = {}
        for t in self.flowsheet().config.time:
            for x in self.length_domain:
                assert 0.0 <= x <= 1.0
                state_args_tx = {}
                for k in state_args_interface_in:
                    if isinstance(state_args_interface_in[k], dict):
                        if k not in state_args_tx:
                            state_args_tx[k] = {}
                        for index in state_args_interface_in[k]:
                            state_args_tx[k][index] = (
                                1.0 - x
                            ) * state_args_interface_in[k][
                                index
                            ] + x * state_args_interface_out[
                                k
                            ][
                                index
                            ]
                    else:
                        state_args_tx[k] = (1.0 - x) * state_args_interface_in[
                            k
                        ] + x * state_args_interface_out[k]
                state_args_interface[t, x] = state_args_tx

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

        return {
            "feed_side": state_args,
            "retentate": state_args_retentate,
            "permeate": state_args_permeate,
            "interface": state_args_interface,
        }

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
        if optarg is None:
            optarg = {}
        # Create solver
        opt = get_solver(solver, optarg)

        # Get inlet state if not provided
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="membrane_channel")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="membrane_channel")

        state_args = self._get_state_args(initialize_guess, state_args)

        # intialize self.properties
        source_flags = super().initialize(state_args=state_args["feed_side"], outlvl=outlvl, optarg=optarg, solver=solver, hold_state=True)

        self.properties_interface.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["interface"],
        )

        self.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["permeate"],
        )

        self.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["permeate"],
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            # occasionally it might be worth retrying a solve
            if not check_optimal_termination(res):
                init_log.warning(
                    "Trouble solving ReverseOsmosis1D unit model, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)
        if not check_optimal_termination(res):
            raise InitializationError("Membrane channel {self.name} failed to initailize")

        init_log.info("Initialization Complete")

        self.release_state(source_flags, outlvl)
    


    def _get_performance_contents(self, time_point=0):
        x_in = self.first_element
        x_interface_in = self.difference_elements.first()
        x_out = self.length_domain.last()
        feed_inlet = self.properties[time_point, x_in]
        feed_outlet = self.properties[time_point, x_out]
        interface_inlet = self.properties_interface[
            time_point, x_interface_in
        ]
        interface_outlet = self.properties_interface[time_point, x_out]
        permeate = self.mixed_permeate[time_point]
        var_dict = {}
        expr_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[
            time_point, "Liq"
        ]
        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[
            time_point, "Liq", "H2O"
        ]
        var_dict["Membrane Area"] = self.area_total
        if hasattr(self, "length"):
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width"):
            var_dict["Membrane Width"] = self.width
        if hasattr(self, "pressure_change_total"):
            var_dict["Pressure Change"] = self.pressure_change_total[time_point]
        if hasattr(self, "N_Re"):
            var_dict["Reynolds Number @Inlet"] = self.N_Re[time_point, x_in]
            var_dict["Reynolds Number @Outlet"] = self.N_Re[time_point, x_out]
        if hasattr(self, "velocity"):
            var_dict["Velocity @Inlet"] = self.velocity[time_point, x_in]
            var_dict["Velocity @Outlet"] = self.velocity[time_point, x_out]
        for j in self.config.property_package.solute_set:
            if interface_inlet.is_property_constructed("conc_mass_phase_comp"):
                var_dict[
                    f"{j} Concentration @Inlet,Membrane-Interface "
                ] = interface_inlet.conc_mass_phase_comp["Liq", j]
            if interface_outlet.is_property_constructed("conc_mass_phase_comp"):
                var_dict[
                    f"{j} Concentration @Outlet,Membrane-Interface "
                ] = interface_outlet.conc_mass_phase_comp["Liq", j]
            if feed_inlet.is_property_constructed("conc_mass_phase_comp"):
                var_dict[
                    f"{j} Concentration @Inlet,Bulk"
                ] = feed_inlet.conc_mass_phase_comp["Liq", j]
            if feed_outlet.is_property_constructed("conc_mass_phase_comp"):
                var_dict[
                    f"{j} Concentration @Outlet,Bulk"
                ] = feed_outlet.conc_mass_phase_comp["Liq", j]
            if permeate.is_property_constructed("conc_mass_phase_comp"):
                var_dict[f"{j} Permeate Concentration"] = permeate.conc_mass_phase_comp[
                    "Liq", j
                ]
        if interface_outlet.is_property_constructed("pressure_osm_phase"):
            var_dict[
                "Osmotic Pressure @Outlet,Membrane-Interface "
            ] = interface_outlet.pressure_osm_phase["Liq"]
        if feed_outlet.is_property_constructed("pressure_osm_phase"):
            var_dict["Osmotic Pressure @Outlet,Bulk"] = feed_outlet.pressure_osm_phase[
                "Liq"
            ]
        if interface_inlet.is_property_constructed("pressure_osm_phase"):
            var_dict[
                "Osmotic Pressure @Inlet,Membrane-Interface"
            ] = interface_inlet.pressure_osm_phase["Liq"]
        if feed_inlet.is_property_constructed("pressure_osm_phase"):
            var_dict["Osmotic Pressure @Inlet,Bulk"] = feed_inlet.pressure_osm_phase[
                "Liq"
            ]
        if feed_inlet.is_property_constructed("flow_vol_phase"):
            var_dict["Volumetric Flowrate @Inlet"] = feed_inlet.flow_vol_phase["Liq"]
        if feed_outlet.is_property_constructed("flow_vol_phase"):
            var_dict["Volumetric Flowrate @Outlet"] = feed_outlet.flow_vol_phase["Liq"]
        if hasattr(self, "dh"):
            var_dict["Hydraulic Diameter"] = self.dh

        expr_dict["Average Solvent Flux (LMH)"] = (
            self.flux_mass_phase_comp_avg[time_point, "Liq", "H2O"] * 3.6e3
        )
        expr_dict["Average Reynolds Number"] = self.N_Re_avg[time_point]
        for j in self.config.property_package.solute_set:
            expr_dict[f"{j} Average Solute Flux (GMH)"] = (
                self.flux_mass_phase_comp_avg[time_point, "Liq", j] * 3.6e6
            )
            expr_dict[f"{j} Average Mass Transfer Coefficient (mm/h)"] = (
                self.K_avg[time_point, j] * 3.6e6
            )

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
        if iscale.get_scaling_factor(self.area_total) is None:
            sf = iscale.get_scaling_factor(self.area_total, default=10, warning=True)
            iscale.set_scaling_factor(self.area_total, sf)

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

        if hasattr(self, "N_Sc"):
            for t, x in self.N_Sc.keys():
                if iscale.get_scaling_factor(self.N_Sc[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sc[t, x], 1e-2)

        if hasattr(self, "N_Sh"):
            for t, x in self.N_Sh.keys():
                if iscale.get_scaling_factor(self.N_Sh[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sh[t, x], 1e-2)

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

        # for self._rescale_permeate_variable
        self._permeate_scaled_properties = ComponentSet()

        for sb in (self.permeate_side, self.mixed_permeate):
            for blk in sb.values():
                for j in self.config.property_package.solute_set:
                    self._rescale_permeate_variable(blk.flow_mass_phase_comp["Liq", j])
                    if blk.is_property_constructed("mass_frac_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.mass_frac_phase_comp["Liq", j]
                        )
                    if blk.is_property_constructed("conc_mass_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.conc_mass_phase_comp["Liq", j]
                        )
                    if blk.is_property_constructed("mole_frac_phase_comp"):
                        self._rescale_permeate_variable(blk.mole_frac_phase_comp[j])
                    if blk.is_property_constructed("molality_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.molality_phase_comp["Liq", j]
                        )
                if blk.is_property_constructed("pressure_osm_phase"):
                    self._rescale_permeate_variable(blk.pressure_osm_phase["Liq"])

        for (t, x, p, j), v in self.flux_mass_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():  # scaling based on solvent flux equation
                    sf = (
                        iscale.get_scaling_factor(self.A_comp[t, j])
                        * iscale.get_scaling_factor(self.dens_solvent)
                        * iscale.get_scaling_factor(
                            self.properties[t, x].pressure
                        )
                    )
                    iscale.set_scaling_factor(v, sf)
                elif comp.is_solute():  # scaling based on solute flux equation
                    sf = iscale.get_scaling_factor(
                        self.B_comp[t, j]
                    ) * iscale.get_scaling_factor(
                        self.properties[t, x].conc_mass_phase_comp[p, j]
                    )
                    iscale.set_scaling_factor(v, sf)


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
