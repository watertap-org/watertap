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


# Import Pyomo libraries
from pyomo.environ import (Var,
                           Param,
                           NonNegativeReals,
                           NegativeReals,
                           units as pyunits,
                           exp,
                           value,
                           Constraint,
                          )
from pyomo.common.config import ConfigValue, In
# Import IDAES cores
from idaes.core import (ControlVolume1DBlock,
                        declare_process_block_class,
                        MomentumBalanceType,
                        useDefault)
from idaes.core.control_volume1d import DistributedVars
from idaes.core.util.misc import add_object_reference
from idaes.core.util import get_solver, scaling as iscale
from idaes.core.util.initialization import solve_indexed_blocks
from watertap.core.util.initialization import check_solve, check_dof
from watertap.unit_models._reverse_osmosis_base import (ConcentrationPolarizationType,
        MassTransferCoefficient,
        PressureChangeType,
        _ReverseOsmosisBaseData)
import idaes.logger as idaeslog


__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis1D")
class ReverseOsmosis1DData(_ReverseOsmosisBaseData):
    """Standard 1D Reverse Osmosis Unit Model Class."""

    CONFIG = _ReverseOsmosisBaseData.CONFIG()

    CONFIG.declare("area_definition", ConfigValue(
            default=DistributedVars.uniform,
            domain=In(DistributedVars),
            description="Argument for defining form of area variable",
            doc="""Argument defining whether area variable should be spatially
    variant or not. **default** - DistributedVars.uniform.
    **Valid values:** {
    DistributedVars.uniform - area does not vary across spatial domain,
    DistributedVars.variant - area can vary over the domain and is indexed
    by time and space.}"""))

    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
    documentation for supported transformations."""))

    CONFIG.declare("transformation_scheme", ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See
    Pyomo documentation for supported schemes."""))

    CONFIG.declare("finite_elements", ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements in length domain",
            doc="""Number of finite elements to use when discretizing length 
            domain (default=20)"""))

    CONFIG.declare("collocation_points", ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=5)"""))

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


    def _process_config(self):
        #TODO: add config errors here:
        if self.config.transformation_method is useDefault:
            _log.warning(
                "Discretization method was "
                "not specified for the "
                "reverse osmosis module. "
                "Defaulting to finite "
                "difference method."
            )
            self.config.transformation_method = "dae.finite_difference"

        if self.config.transformation_scheme is useDefault:
            _log.warning(
                "Discretization scheme was "
                "not specified for the "
                "reverse osmosis module."
                "Defaulting to backward finite "
                "difference."
            )
            self.config.transformation_scheme = "BACKWARD"

    def build(self):
        """
        Build 1D RO model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # Check configuration errors
        self._process_config()

        # Build 1D Control volume for feed side
        self.feed_side = feed_side = ControlVolume1DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": self.config.area_definition,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args,
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points
        })

        # Add geometry to feed side
        feed_side.add_geometry()
        # Add state blocks to feed side
        feed_side.add_state_blocks(has_phase_equilibrium=False)
        # Populate feed side
        feed_side.add_material_balances(balance_type=self.config.material_balance_type,
                                        has_mass_transfer=True)
        feed_side.add_momentum_balances(balance_type=self.config.momentum_balance_type,
                                        has_pressure_change=self.config.has_pressure_change)
        # Apply transformation to feed side
        feed_side.apply_transformation()
        # Add inlet/outlet ports for feed side
        self.add_inlet_port(name="inlet", block=feed_side)
        self.add_outlet_port(name="retentate", block=feed_side)
        # Make indexed stateblock and separate stateblock for permeate-side and permeate outlet, respectively.
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets
        self.permeate_side = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.feed_side.length_domain,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict)
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict)

        # Membrane interface: indexed state block
        self.feed_side.properties_interface = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.feed_side.length_domain,
            doc="Material properties of feed-side membrane interface",
            default=tmp_dict)

        # Add port to mixed_permeate
        self.add_port(name="permeate", block=self.mixed_permeate)

        # ==========================================================================
        """ Add references to control volume geometry."""
        add_object_reference(self, 'length', feed_side.length)
        add_object_reference(self, 'area_cross', feed_side.area)


        # Add reference to pressure drop for feed side only
        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != MomentumBalanceType.none):
            add_object_reference(self, 'deltaP', feed_side.deltaP)
            self.deltaP.setub(0)

        self._make_performance()

        self._add_expressions()

        self._get_performance_contents()

    def _add_expressions(self):
        """
        Generate expressions for additional results desired for full report
        """
        solute_set = self.config.property_package.solute_set

        if self.config.has_full_reporting is False:
            pass
        else:
            @self.Expression(self.flowsheet().config.time,
                             self.config.property_package.phase_list,
                             self.config.property_package.component_list,
                             doc="Average flux expression")
            def flux_mass_phase_comp_avg(b, t, p, j):
                return sum(b.flux_mass_phase_comp[t, x, p, j]
                           for x in self.feed_side.length_domain
                           if x > 0) / self.nfe
            if hasattr(self, 'N_Re'):
                @self.Expression(self.flowsheet().config.time,
                                 doc="Average Reynolds Number expression")
                def N_Re_avg(b, t):
                    return sum(b.N_Re[t, x]
                               for x in self.feed_side.length_domain) / self.nfe
            if hasattr(self, 'Kf'):
                @self.Expression(self.flowsheet().config.time,
                                 solute_set,
                                 doc="Average mass transfer coefficient expression")
                def Kf_avg(b, t, j):
                    return sum(b.Kf[t, x, j]
                               for x in self.feed_side.length_domain
                               if x > 0) / self.nfe

    def _make_performance(self):
        """
        Variables and constraints for unit model.

        Args:
            None

        Returns:
            None
        """

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        # Units
        units_meta = \
            self.config.property_package.get_metadata().get_derived_units

        self.nfe = Param(
            initialize=(len(self.feed_side.length_domain)-1),
            units=pyunits.dimensionless,
            doc="Number of finite elements")

        # ==========================================================================
        self.flux_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.feed_side.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b,t,x,p,j : 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b,t,x,p,j : (1e-4, 3e-2) if j in solvent_set else (1e-8, 1e-3),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Mass flux across membrane')

        self.width = Var(
            initialize=1,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta('length'),
            doc='Membrane width')

        # mass transfer
        # TODO: replace self.recovery_vol_phase[t, 'Liq'] w/self.recovery_mass_phase_comp[t, 'Liq', j])
        def mass_transfer_phase_comp_initialize(b, t, x, p, j):
            return value(self.feed_side.properties[t, x].get_material_flow_terms('Liq', j)
                         * self.recovery_mass_phase_comp[t, 'Liq', j])

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.feed_side.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass') * units_meta('time')**-1 * units_meta('length')**-1,
            doc='Mass transfer to permeate')

        if self.config.has_pressure_change:
            self.deltaP_stage = Var(
                self.flowsheet().config.time,
                initialize=-1e5,
                bounds=(-1e6, 0),
                domain=NegativeReals,
                units=units_meta('pressure'),
                doc='Pressure drop across unit')

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.fixed:
            self.cp_modulus = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                solute_set,
                initialize=1.1,
                bounds=(0.9, 3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc='Concentration polarization modulus')

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated:
            self.Kf = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta('length') * units_meta('time')**-1,
                doc='Mass transfer coefficient in feed channel')

        if ((self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
                or (self.config.pressure_change_type == PressureChangeType.calculated
                    and self.config.has_pressure_change)):
            self.N_Re = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                initialize=5e2,
                bounds=(10, 5e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Reynolds number in feed channel")
        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            self.N_Sc = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                initialize=5e2,
                bounds=(1e2, 2e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Schmidt number in feed channel")
            self.N_Sh = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                initialize=1e2,
                bounds=(1, 3e2),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Sherwood number in feed channel")
        if (self.config.pressure_change_type == PressureChangeType.calculated
                and self.config.has_pressure_change):
            self.velocity = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=units_meta('length')/units_meta('time'),
                doc="Crossflow velocity in feed channel")
            self.friction_factor_darcy = Var(
                self.flowsheet().config.time,
                self.feed_side.length_domain,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Darcy friction factor in feed channel")

        # ==========================================================================
        # Volumetric Recovery rate

        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (b.recovery_vol_phase[t, 'Liq'] ==
                    b.mixed_permeate[t].flow_vol_phase['Liq'] /
                    b.feed_side.properties[t, self.feed_side.length_domain.first()].flow_vol_phase['Liq'])

        # ==========================================================================
        # Mass-based Component Recovery rate

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list)
        def eq_recovery_mass_phase_comp(b, t, j):
            return (b.recovery_mass_phase_comp[t, 'Liq', j]
                    * b.feed_side.properties[t, b.feed_side.length_domain.first()].flow_mass_phase_comp['Liq', j] ==
                    b.mixed_permeate[t].flow_mass_phase_comp['Liq', j])

        # ==========================================================================
        # Mass transfer term equation

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(b, t, x, p, j):
            if x == b.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                return b.mass_transfer_phase_comp[t, x, p, j] == -b.feed_side.mass_transfer_term[t, x, p, j]
        # ==========================================================================
        # Membrane area equation

        @self.Constraint(doc="Membrane area")
        def eq_area(b):
            return b.area == b.length * b.width
        # ==========================================================================
        # Mass flux = feed mass transfer equation

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_flux_equal_mass_transfer(b, t, x, p, j):
            if x == b.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                return b.flux_mass_phase_comp[t, x, p, j] * b.width == -b.feed_side.mass_transfer_term[t, x, p, j]
        # ==========================================================================
        # Mass flux equations (Jw and Js)

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Solvent and solute mass flux")
        def eq_flux_mass(b, t, x, p, j):
            if x == b.feed_side.length_domain.first():
                return Constraint.Skip
            else:
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

        # ==========================================================================
        # Final permeate mass flow rate (of solvent and solute) --> Mp,j, final = sum(Mp,j)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Permeate mass flow rates exiting unit")
        def eq_permeate_production(b, t, p, j):
            return (b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == sum(b.permeate_side[t, x].get_material_flow_terms(p, j)
                           for x in b.feed_side.length_domain if x != 0))
        # ==========================================================================
        # Feed and permeate-side mass transfer connection --> Mp,j = Mf,transfer = Jj * W * L/n

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, x, p, j):
            if x == b.feed_side.length_domain.first():
                return b.permeate_side[t, x].get_material_flow_terms(p, j) == 0
            else:
                return (b.permeate_side[t, x].get_material_flow_terms(p, j)
                        == -b.feed_side.mass_transfer_term[t, x, p, j] * b.length / b.nfe)

        # # ==========================================================================
        # Concentration polarization

        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.feed_side.length_domain,
                                   solute_set,
                                   doc="Concentration polarization")
        def eq_concentration_polarization(b, t, x, j):
            if x == self.feed_side.length_domain.first():
                return Constraint.Skip
            else:
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
        # Constraints active when MassTransferCoefficient.calculated
        # Mass transfer coefficient calculation
        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             solute_set,
                             doc="Mass transfer coefficient in feed channel")
            def eq_Kf(b, t, x, j):
                if x == self.feed_side.length_domain.first():
                    return Constraint.Skip
                else:
                    bulk = b.feed_side.properties[t, x]
                    return (b.Kf[t, x, j] * b.dh ==
                            bulk.diffus_phase['Liq']  # TODO: add diff coefficient to SW prop and consider multi-components
                            * b.N_Sh[t, x])

            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="Sherwood number")
            def eq_N_Sh(b, t, x):
                return (b.N_Sh[t, x] ==
                        0.46 * (b.N_Re[t, x] * b.N_Sc[t, x])**0.36)

            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="Schmidt number")
            def eq_N_Sc(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.N_Sc[t, x] * bulk.dens_mass_phase['Liq'] * bulk.diffus_phase['Liq'] ==
                        bulk.visc_d_phase['Liq'])

        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            or (self.config.pressure_change_type == PressureChangeType.calculated
                and self.config.has_pressure_change)):
            @self.Constraint(doc="Cross-sectional area")
            def eq_area_cross(b):
                return b.area_cross == b.channel_height * b.width * b.spacer_porosity

            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="Reynolds number")
            def eq_N_Re(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.N_Re[t, x] * b.area_cross * bulk.visc_d_phase['Liq'] ==
                        sum(bulk.flow_mass_phase_comp['Liq', j] for j in b.config.property_package.component_list)
                        * b.dh)

        ## ==========================================================================
        # Pressure drop
        if ((self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length
             or self.config.pressure_change_type == PressureChangeType.calculated)
                and self.config.has_pressure_change):
            @self.Constraint(self.flowsheet().config.time,
                             doc='Pressure drop across unit')
            def eq_pressure_drop(b, t):
                return (b.deltaP_stage[t] ==
                        sum(b.deltaP[t, x] * b.length / b.nfe
                            for x in b.feed_side.length_domain if x != 0))

        if (self.config.pressure_change_type == PressureChangeType.fixed_per_stage
                and self.config.has_pressure_change):
            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc='Fixed pressure drop across unit')
            def eq_pressure_drop(b, t, x):
                return b.deltaP_stage[t] == b.length * b.deltaP[t, x]

        if (self.config.pressure_change_type == PressureChangeType.calculated
                and self.config.has_pressure_change):
            ## ==========================================================================
            # Crossflow velocity
            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="Crossflow velocity constraint")
            def eq_velocity(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return b.velocity[t, x] * b.area_cross == bulk.flow_vol_phase['Liq']
            ## ==========================================================================
            # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
            # TODO: this relationship for friction factor is specific to a particular spacer geometry. Add alternatives.

            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="Darcy friction factor constraint")
            def eq_friction_factor_darcy(b, t, x):
                return (b.friction_factor_darcy[t, x] - 0.42) * b.N_Re[t, x] == 189.3
            ## ==========================================================================
            # Pressure change per unit length due to friction
            # -1/2*f/dh*density*velocity^2

            @self.Constraint(self.flowsheet().config.time,
                             self.feed_side.length_domain,
                             doc="pressure change per unit length due to friction")
            def eq_dP_dx(b, t, x):
                bulk = b.feed_side.properties[t, x]
                return (b.deltaP[t, x] * b.dh ==
                        -0.5 * b.friction_factor_darcy[t, x]
                        * bulk.dens_mass_phase['Liq'] * b.velocity[t, x]**2)
        ## ==========================================================================
        # Feed-side isothermal conditions

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         doc="Isothermal assumption for permeate")
        def eq_feed_isothermal(b, t, x):
            if x == b.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                return b.feed_side.properties[t, b.feed_side.length_domain.first()].temperature == \
                       b.feed_side.properties[t, x].temperature
        # # ==========================================================================
        # Feed and permeate-side isothermal conditions

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t, x):
            return b.feed_side.properties[t, x].temperature == \
                   b.permeate_side[t, x].temperature
        # ==========================================================================
        # isothermal conditions at permeate outlet

        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate out")
        def eq_permeate_outlet_isothermal(b, t):
            return b.feed_side.properties[t, b.feed_side.length_domain.first()].temperature == \
                   b.mixed_permeate[t].temperature
        # ==========================================================================
        # isobaric conditions across permeate channel and at permeate outlet

        @self.Constraint(self.flowsheet().config.time,
                         self.feed_side.length_domain,
                         doc="Isobaric assumption for permeate out")
        def eq_permeate_outlet_isobaric(b, t, x):
            return b.permeate_side[t, x].pressure == \
                   b.mixed_permeate[t].pressure
        # ==========================================================================
        # Bulk and interface connections on the feed-side
        # TEMPERATURE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.feed_side.length_domain,
                                   doc="Temperature at interface")
        def eq_equal_temp_interface(b, t, x):
            if x == self.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                bulk = b.properties[t, x]
                interface = b.properties_interface[t, x]
            return interface.temperature == \
                   bulk.temperature
        # PRESSURE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.feed_side.length_domain,
                                   doc="Pressure at interface")
        def eq_equal_pressure_interface(b, t, x):
            if x == self.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                bulk = b.properties[t, x]
                interface = b.properties_interface[t, x]
            return interface.pressure == \
                   bulk.pressure
        # VOLUMETRIC FLOWRATE
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.feed_side.length_domain,
                                   doc="Volumetric flow at interface of inlet")
        def eq_equal_flow_vol_interface(b, t, x):
            if x == self.feed_side.length_domain.first():
                return Constraint.Skip
            else:
                bulk = b.properties[t, x]
                interface = b.properties_interface[t, x]
            return interface.flow_vol_phase['Liq'] ==\
                   bulk.flow_vol_phase['Liq']

    def initialize(blk,
                   feed_side_args=None,
                   permeate_side_args=None,
                   permeate_block_args=None,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None,
                   fail_on_warning=False,
                   ignore_dof=False):
        """
        Initialization routine for 1D-RO unit.

        Keyword Arguments:
            feed_side_args : a dict of arguments to be passed to the property
             package(s) of the feed_side to provide an initial state for
             initialization (see documentation of the specific
             property package)
            permeate_side_args : a dict of arguments to be passed to the property
             package(s) of the permeate_side to provide an initial state for
             initialization (see documentation of the specific
             property package)
            permeate_block_args : a dict of arguments to be passed to the property
             package(s) of the final permeate StateBlock to provide an initial state for
             initialization (see documentation of the specific
             property package)
            outlvl : sets output level of initialization routine
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            optarg : solver options dictionary object (default=None, use default solver options)
            fail_on_warning : boolean argument to fail or only produce  warning upon unsuccessful solve (default=False)
            ignore_dof : boolean argument to ignore when DOF != 0 (default=False)
        Returns:
            None

        """

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        if optarg is None:
            optarg = {'bound_push': 1e-8}

        opt = get_solver(solver, optarg)

        init_log.info('Starting Initialization Step 1: initialize blocks.')
        # ---------------------------------------------------------------------
        # Step 1: Initialize feed_side, permeate_side, and mixed_permeate blocks
        flags_feed_side = blk.feed_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=feed_side_args)
        init_log.info('Feed-side initialization complete. Initialize permeate-side. ')
        flags_permeate_side = blk.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=permeate_side_args)
        init_log.info('Permeate-side initialization complete. Initialize permeate outlet. ')
        flags_mixed_permeate = blk.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=permeate_block_args)
        init_log.info('Permeate outlet initialization complete. Initialize permeate outlet. ')

        if not ignore_dof:
           check_dof(blk, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Step 2: Solve unit
        init_log.info('Initialization Step 1 complete: all state blocks initialized.'
                      'Starting Initialization Step 2: solve indexed blocks.')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        # only fail on the final solve
        check_solve(results, logger=init_log, fail_flag=False, checkpoint='Initialization Step 2: solve indexed blocks')
        init_log.info('Starting Initialization Step 3: perform final solve.')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        check_solve(res, logger=init_log, fail_flag=fail_on_warning, checkpoint='Initialization Step 3: final solve')
        # Release Inlet state
        blk.feed_side.release_state(flags_feed_side, outlvl)

    def _get_performance_contents(self, time_point=0):
        x_in = self.feed_side.length_domain.first()
        x_interface_in = self.feed_side.length_domain.at(2)
        x_out = self.feed_side.length_domain.last()
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
        if hasattr(self, "length") or self.config.has_full_reporting:
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width") or self.config.has_full_reporting:
            var_dict["Membrane Width"] = self.width
        if hasattr(self, "deltaP_stage") or self.config.has_full_reporting:
            var_dict["Pressure Change"] = self.deltaP_stage[time_point]
        if hasattr(self, "N_Re") or self.config.has_full_reporting:
            var_dict["Reynolds Number @Inlet"] = self.N_Re[time_point, x_in]
            var_dict["Reynolds Number @Outlet"] = self.N_Re[time_point, x_out]
        if hasattr(self, "velocity") or self.config.has_full_reporting:
            var_dict["Velocity @Inlet"] = self.velocity[time_point, x_in]
            var_dict["Velocity @Outlet"] = self.velocity[time_point, x_out]
        for j in self.config.property_package.solute_set:
            if interface_inlet.is_property_constructed('conc_mass_phase_comp') or self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Inlet,Membrane-Interface '] = (
                    interface_inlet.conc_mass_phase_comp['Liq', j])
            if interface_outlet.is_property_constructed('conc_mass_phase_comp') or self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Outlet,Membrane-Interface '] = (
                    interface_outlet.conc_mass_phase_comp['Liq', j])
            if feed_inlet.is_property_constructed('conc_mass_phase_comp') or self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Inlet,Bulk'] = (
                    feed_inlet.conc_mass_phase_comp['Liq', j])
            if feed_outlet.is_property_constructed('conc_mass_phase_comp') or self.config.has_full_reporting:
                var_dict[f'{j} Concentration @Outlet,Bulk'] = (
                    feed_outlet.conc_mass_phase_comp['Liq', j])
            if permeate.is_property_constructed('conc_mass_phase_comp') or self.config.has_full_reporting:
                var_dict[f'{j} Permeate Concentration'] = (
                    permeate.conc_mass_phase_comp['Liq', j])
        if interface_outlet.is_property_constructed('pressure_osm') or self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Outlet,Membrane-Interface '] = (
                interface_outlet.pressure_osm)
        if feed_outlet.is_property_constructed('pressure_osm') or self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Outlet,Bulk'] = (
                feed_outlet.pressure_osm)
        if interface_inlet.is_property_constructed('pressure_osm') or self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Inlet,Membrane-Interface'] = (
                interface_inlet.pressure_osm)
        if feed_inlet.is_property_constructed('pressure_osm') or self.config.has_full_reporting:
            var_dict['Osmotic Pressure @Inlet,Bulk'] = (
                feed_inlet.pressure_osm)
        if feed_inlet.is_property_constructed('flow_vol_phase') or self.config.has_full_reporting:
            var_dict['Volumetric Flowrate @Inlet'] = (
                feed_inlet.flow_vol_phase['Liq'])
        if feed_outlet.is_property_constructed('flow_vol_phase') or self.config.has_full_reporting:
            var_dict['Volumetric Flowrate @Outlet'] = (
                feed_outlet.flow_vol_phase['Liq'])
        if hasattr(self, 'dh') or self.config.has_full_reporting:
            var_dict["Hydraulic Diameter"] = self.dh

        if self.config.has_full_reporting:
            expr_dict['Average Solvent Flux (LMH)'] = self.flux_mass_phase_comp_avg[time_point, 'Liq', 'H2O'] * 3.6e3
            expr_dict['Average Reynolds Number'] = self.N_Re_avg[time_point]
            for j in self.config.property_package.solute_set:
                expr_dict[f'{j} Average Solute Flux (GMH)'] = self.flux_mass_phase_comp_avg[time_point, 'Liq', j] * 3.6e6
                expr_dict[f'{j} Average Mass Transfer Coefficient (mm/h)'] = self.Kf_avg[time_point, j] * 3.6e6


        # TODO: add more vars
        return {"vars": var_dict, "exprs": expr_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # setting scaling factors for variables
        for j in self.config.property_package.component_list:
            iscale.set_scaling_factor(self.permeate_side[0, 0].flow_mass_phase_comp['Liq', j], 1e+5)

        # these variables should have user input, if not there will be a warning

        if iscale.get_scaling_factor(self.width) is None:
            sf = iscale.get_scaling_factor(self.width, default=1, warning=True)
            iscale.set_scaling_factor(self.width, sf)

        if iscale.get_scaling_factor(self.length) is None:
            sf = iscale.get_scaling_factor(self.length, default=10, warning=True)
            iscale.set_scaling_factor(self.length, sf)

        # will not override if the user provides the scaling factor
        ## default of 1 set by ControlVolume1D
        if iscale.get_scaling_factor(self.area_cross) == 1:
            iscale.set_scaling_factor(self.area_cross, 100)

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties[0, 0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        for (t, x, p, j), v in self.flux_mass_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if x == self.feed_side.length_domain.first():
                    if comp.is_solvent():
                        iscale.set_scaling_factor(v, 5e4)  # inverse of initial value from flux_mass_phase_comp_initialize
                    elif comp.is_solute():
                        iscale.set_scaling_factor(v, 1e6)  # inverse of initial value from flux_mass_phase_comp_initialize
                else:
                    if comp.is_solvent():  # scaling based on solvent flux equation
                        sf = (iscale.get_scaling_factor(self.A_comp[t, j])
                              * iscale.get_scaling_factor(self.dens_solvent)
                              * iscale.get_scaling_factor(self.feed_side.properties[t, x].pressure))
                        iscale.set_scaling_factor(v, sf)
                    elif comp.is_solute():  # scaling based on solute flux equation
                        sf = (iscale.get_scaling_factor(self.B_comp[t, j])
                              * iscale.get_scaling_factor(self.feed_side.properties[t, x].conc_mass_phase_comp[p, j]))
                        iscale.set_scaling_factor(v, sf)

        if hasattr(self, 'cp_modulus'):
            for v in self.cp_modulus.values():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(v)
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, 'Kf'):
            for v in self.Kf.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e5)

        if hasattr(self, 'N_Re'):
            for t, x in self.N_Re.keys():
                if iscale.get_scaling_factor(self.N_Re[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Re[t, x], 1e-1)

        if hasattr(self, 'N_Sc'):
            for t, x in self.N_Sc.keys():
                if iscale.get_scaling_factor(self.N_Sc[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sc[t, x], 1e-1)

        if hasattr(self, 'N_Sh'):
            for t, x in self.N_Sh.keys():
                if iscale.get_scaling_factor(self.N_Sh[t, x]) is None:
                     iscale.set_scaling_factor(self.N_Sh[t, x], 1e-2)

        if hasattr(self, 'deltaP_stage'):
            for v in self.deltaP_stage.values():
                if iscale.get_scaling_factor(v) is None:
                     iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, 'velocity'):
            for v in self.velocity.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'friction_factor_darcy'):
            for v in self.friction_factor_darcy.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        for (t, x, p, j), v in self.eq_mass_flux_equal_mass_transfer.items():
            if iscale.get_scaling_factor(v) is None:
                if x == self.feed_side.length_domain.first():
                    pass
                else:
                    sf = iscale.get_scaling_factor(self.flux_mass_phase_comp[t, x, p, j])\
                         * iscale.get_scaling_factor(self.width)
                    iscale.set_scaling_factor(v, sf/100.)

        for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties[t, x].get_material_flow_terms(p, j)) \
                     / iscale.get_scaling_factor(self.feed_side.length)
                if x == 0:
                    sf *= 10.
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, 'deltaP'):
            for v in self.feed_side.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e-5)
        else:
            for v in self.feed_side.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e5)

        # Scale constraints
        for ind, c in self.eq_mass_transfer_term.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_connect_mass_transfer.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf*10.)

        sf = iscale.get_scaling_factor(self.area)
        iscale.constraint_scaling_transform(self.eq_area, sf/10.)

        for ind, c in self.eq_permeate_production.items():
            # TODO: revise this scaling factor; setting to 100 for now
            iscale.constraint_scaling_transform(c, 100.)

        for ind, c in self.eq_flux_mass.items():
            sf = iscale.get_scaling_factor(self.flux_mass_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.eq_feed_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties[t, x].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.eq_permeate_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties[t, x].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_permeate_outlet_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties[t, 0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.eq_permeate_outlet_isobaric.items():
            sf = iscale.get_scaling_factor(self.permeate_side[t, x].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_recovery_vol_phase.items():
            iscale.constraint_scaling_transform(self.eq_recovery_vol_phase[t], 1)

        for (t, j), c in self.eq_recovery_mass_phase_comp.items():
            sf = (iscale.get_scaling_factor(self.recovery_mass_phase_comp[t, 'Liq', j])
                  * iscale.get_scaling_factor(self.inlet.flow_mass_phase_comp[0, 'Liq', j]))
            iscale.constraint_scaling_transform(c, sf)

        for (t, x, j), c in self.feed_side.eq_concentration_polarization.items():
            prop_interface = self.feed_side.properties_interface[t, x]
            sf = iscale.get_scaling_factor(prop_interface.conc_mass_phase_comp['Liq', j])
            iscale.constraint_scaling_transform(c, sf*10.)

        for (t, x), c in self.feed_side.eq_equal_temp_interface.items():
            prop_interface = self.feed_side.properties_interface[t, x]
            sf = iscale.get_scaling_factor(prop_interface.temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.feed_side.eq_equal_pressure_interface.items():
            prop_interface = self.feed_side.properties_interface[t, x]
            sf = iscale.get_scaling_factor(prop_interface.pressure)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.feed_side.eq_equal_flow_vol_interface.items():
            prop_interface = self.feed_side.properties_interface[t, x]
            sf = iscale.get_scaling_factor(prop_interface.flow_vol_phase['Liq'])
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_Kf'):
            for ind, c in self.eq_Kf.items():
                sf = iscale.get_scaling_factor(self.Kf[ind])
                iscale.constraint_scaling_transform(c, sf*10.)

        if hasattr(self, 'eq_N_Re'):
            for ind, c in self.eq_N_Re.items():
                sf = iscale.get_scaling_factor(self.N_Re[ind])
                iscale.constraint_scaling_transform(c, sf*1e4)

        if hasattr(self, 'eq_N_Sc'):
            for ind, c in self.eq_N_Sc.items():
                sf = iscale.get_scaling_factor(self.N_Sc[ind])
                iscale.constraint_scaling_transform(c, sf*1e4)

        if hasattr(self, 'eq_N_Sh'):
            for ind, c in self.eq_N_Sh.items():
                sf = iscale.get_scaling_factor(self.N_Sh[ind])
                iscale.constraint_scaling_transform(c, sf*1e1)

        if hasattr(self, 'eq_area_cross'):
            sf = iscale.get_scaling_factor(self.area_cross)
            iscale.constraint_scaling_transform(self.eq_area_cross, sf*10.)

        if hasattr(self, 'eq_pressure_drop'):
            if (self.config.pressure_change_type == PressureChangeType.calculated
                or self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length):
                for t, c in self.eq_pressure_drop.items():
                    sf = iscale.get_scaling_factor(self.deltaP_stage[t])
                    iscale.constraint_scaling_transform(c, sf)
            elif self.config.pressure_change_type == PressureChangeType.fixed_per_stage:
                for (t, x), c in self.eq_pressure_drop.items():
                    sf = iscale.get_scaling_factor(self.deltaP_stage[t])
                    iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_velocity'):
            for ind, c in self.eq_velocity.items():
                sf = iscale.get_scaling_factor(self.velocity[ind])
                iscale.constraint_scaling_transform(c, sf*1e4)

        if hasattr(self, 'eq_friction_factor_darcy'):
            for ind, c in self.eq_friction_factor_darcy.items():
                sf = iscale.get_scaling_factor(self.friction_factor_darcy[ind])
                iscale.constraint_scaling_transform(c, sf/10.)

        if hasattr(self, 'eq_dP_dx'):
            for ind, c in self.eq_dP_dx.items():
                sf = (iscale.get_scaling_factor(self.deltaP[ind])
                      * iscale.get_scaling_factor(self.dh))
                iscale.constraint_scaling_transform(c, sf)
