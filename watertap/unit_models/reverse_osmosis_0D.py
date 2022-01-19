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
                           Set,
                           NonNegativeReals,
                           NegativeReals,
                           Reference,
                           units as pyunits,
                           exp,
                           value,
                           check_optimal_termination)
# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import get_solver
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
from watertap.core.util.initialization import check_solve, check_dof
from watertap.unit_models._reverse_osmosis_base import (ConcentrationPolarizationType,
        MassTransferCoefficient,
        PressureChangeType,
        _ReverseOsmosisBaseData)
import idaes.logger as idaeslog


__author__ = "Tim Bartholomew, Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis0D")
class ReverseOsmosisData(_ReverseOsmosisBaseData):
    """
    Standard RO Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """
    CONFIG = _ReverseOsmosisBaseData.CONFIG()

    def build(self):
        """
        Build the RO model.
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        self.io_list = Set(ordered=True, initialize=('in', 'out'))  # inlet/outlet set
        add_object_reference(self, 'length_domain', self.io_list)

        # Build control volume for feed side
        self.feed_side = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.feed_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True)

        self.feed_side.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=True)

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Add additional state blocks
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets
        # Build permeate side
        self.permeate_side = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.io_list,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict)
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict)
        # Interface properties
        self.feed_side.properties_interface_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed-side interface at inlet",
            default=tmp_dict)
        self.feed_side.properties_interface_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed-side interface at outlet",
            default=tmp_dict)

        # Add Ports
        self.add_inlet_port(name='inlet', block=self.feed_side)
        self.add_outlet_port(name='retentate', block=self.feed_side)
        self.add_port(name='permeate', block=self.mixed_permeate)

        # References for control volume
        # pressure change
        if (self.config.has_pressure_change and
                self.config.momentum_balance_type != 'none'):
            self.deltaP = Reference(self.feed_side.deltaP)

        self._make_performance()

        self._add_expressions()


    def _make_performance(self):
        super()._make_performance()

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.flux_mass_io_phase_comp = Var(
            self.flowsheet().config.time,
            self.io_list,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b,t,x,p,j : 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b,t,x,p,j : (1e-4, 3e-2) if j in solvent_set else (1e-8, 1e-3),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Mass flux across membrane at inlet and outlet')

        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Observed solute rejection')

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(self.feed_side.properties_in[t].get_material_flow_terms('Liq', j)
                         * self.recovery_mass_phase_comp[t, 'Liq', j])
        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass') * units_meta('time')**-1,
            doc='Mass transfer to permeate')

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.fixed:
            self.cp_modulus = Var(
                self.flowsheet().config.time,
                solute_set,
                initialize=1.1,
                bounds=(0.9, 3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc='Concentration polarization modulus')

        if self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated:
            self.Kf_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                solute_set,
                initialize=5e-5,
                bounds=(1e-6, 1e-3),
                domain=NonNegativeReals,
                units=units_meta('length') * units_meta('time')**-1,
                doc='Mass transfer coefficient in feed channel at inlet and outlet')

        if ((self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
                or self.config.pressure_change_type == PressureChangeType.calculated):
            self.N_Re_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=5e2,
                bounds=(10, 5e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Reynolds number at inlet and outlet")

        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            self.N_Sc_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=5e2,
                bounds=(1e2, 2e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Schmidt number at inlet and outlet")
            self.N_Sh_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=1e2,
                bounds=(1, 3e2),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Sherwood number at inlet and outlet")

        if self.config.pressure_change_type == PressureChangeType.calculated:
            self.velocity_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=units_meta('length')/units_meta('time'),
                doc="Crossflow velocity in feed channel at inlet and outlet")
            self.friction_factor_darcy_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=0.5,
                bounds=(1e-2, 5),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Darcy friction factor in feed channel at inlet and outlet")

            # different representation in 1DRO
            self.dP_dx_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta('pressure')*units_meta('length')**-1,
                doc="Pressure drop per unit length of feed channel at inlet and outlet")

        if ((self.config.pressure_change_type != PressureChangeType.fixed_per_stage)
                or (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)):
            # comes from ControlVolume1D in 1DRO
            self.length = Var(
                initialize=10,
                bounds=(0.1, 5e2),
                domain=NonNegativeReals,
                units=units_meta('length'),
                doc='Effective membrane length')
            # not optional in 1DRO
            self.width = Var(
                initialize=1,
                bounds=(0.1, 5e2),
                domain=NonNegativeReals,
                units=units_meta('length'),
                doc='Effective feed-channel width')

        if self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length:
            self.dP_dx = Var(
                self.flowsheet().config.time,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta('pressure')*units_meta('length')**-1,
                doc="pressure drop per unit length across feed channel")

        # constraints for additional variables (i.e. variables not used in other constraints)
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (b.recovery_vol_phase[t, 'Liq'] ==
                    b.mixed_permeate[t].flow_vol_phase['Liq'] /
                    b.feed_side.properties_in[t].flow_vol_phase['Liq'])

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list)
        def eq_recovery_mass_phase_comp(b, t, j):
            return (b.recovery_mass_phase_comp[t, 'Liq', j] ==
                    b.mixed_permeate[t].flow_mass_phase_comp['Liq', j] /
                    b.feed_side.properties_in[t].flow_mass_phase_comp['Liq', j])

        @self.Constraint(self.flowsheet().config.time,
                         solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return (b.rejection_phase_comp[t, 'Liq', j] ==
                    1 - (b.mixed_permeate[t].conc_mass_phase_comp['Liq', j] /
                         b.feed_side.properties_in[t].conc_mass_phase_comp['Liq', j]))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(self, t, p, j):
            return self.mass_transfer_phase_comp[t, p, j] == -self.feed_side.mass_transfer_term[t, p, j]

        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Water and salt flux")
        def eq_flux_io(b, t, x, p, j):
            prop_perm = b.permeate_side[t, x]
            if x == 'in':
                prop_feed = b.feed_side.properties_in[t]
                prop_feed_inter = b.feed_side.properties_interface_in[t]
            elif x == 'out':
                prop_feed = b.feed_side.properties_out[t]
                prop_feed_inter = b.feed_side.properties_interface_out[t]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return (b.flux_mass_io_phase_comp[t, x, p, j] == b.A_comp[t, j] * b.dens_solvent
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - (prop_feed_inter.pressure_osm - prop_perm.pressure_osm)))
            elif comp.is_solute():
                return (b.flux_mass_io_phase_comp[t, x, p, j] == b.B_comp[t, j]
                        * (prop_feed_inter.conc_mass_phase_comp[p, j] - prop_perm.conc_mass_phase_comp[p, j]))

        # RO performance equations (not in 1DRO)
        @self.Expression(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Average flux expression")
        def flux_mass_phase_comp_avg(b, t, p, j):
            return 0.5 * sum(b.flux_mass_io_phase_comp[t, x, p, j] for x in self.io_list)

        # Difference expression in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Permeate production")
        def eq_permeate_production(b, t, p, j):
            return (b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == b.area * b.flux_mass_phase_comp_avg[t, p, j])

        # Feed and permeate-side connection
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, p, j):
            return (b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == -b.feed_side.mass_transfer_term[t, p, j])

        # Non-existent in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         doc="Enthalpy transfer from feed to permeate")
        def eq_connect_enthalpy_transfer(b, t):
            return (b.mixed_permeate[t].get_enthalpy_flow_terms('Liq')
                    == -b.feed_side.enthalpy_transfer[t])

        # # Permeate-side stateblocks
        # Nor in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         solute_set,
                         doc="Permeate mass fraction")
        def eq_mass_frac_permeate_io(b, t, x, j):
            prop_io = b.permeate_side[t,x]
            return (prop_io.mass_frac_phase_comp['Liq', j]
                    * sum(self.flux_mass_io_phase_comp[t, x, 'Liq', jj]
                          for jj in self.config.property_package.component_list)
                    == self.flux_mass_io_phase_comp[t, x, 'Liq', j])
        # # ==========================================================================
        # Feed and permeate-side isothermal conditions

        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t, x):
            prop_io = b.feed_side.properties_in[t] if x == 'in' else b.feed_side.properties_out[t]
            return prop_io.temperature == \
                   b.permeate_side[t, x].temperature
        # ==========================================================================
        # isothermal conditions at permeate outlet (uses properties[t,0] in 1DRO)

        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate out")
        def eq_permeate_outlet_isothermal(b, t):
            return b.feed_side.properties_out[t].temperature == \
                   b.mixed_permeate[t].temperature
        # ==========================================================================
        # isobaric conditions across permeate channel and at permeate outlet

        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         doc="Isobaric assumption for permeate out")
        def eq_permeate_outlet_isobaric(b, t, x):
            return b.permeate_side[t, x].pressure == \
                   b.mixed_permeate[t].pressure

        # not in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         doc="Permeate flowrate")
        def eq_flow_vol_permeate_io(b, t, x):
            prop_io = b.permeate_side[t,x]
            return prop_io.flow_vol_phase['Liq'] == b.mixed_permeate[t].flow_vol_phase['Liq']

        # Concentration polarization (slightly different form in 1DRO)
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   solute_set,
                                   doc="Concentration polarization")
        def eq_concentration_polarization_io(b, t, x, j):
            if x == 'in':
                prop_io = b.properties_in[t]
                prop_interface_io = b.properties_interface_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
                prop_interface_io = b.properties_interface_out[t]
            if self.config.concentration_polarization_type == ConcentrationPolarizationType.none:
                return prop_interface_io.conc_mass_phase_comp['Liq', j] == \
                       prop_io.conc_mass_phase_comp['Liq', j]
            elif self.config.concentration_polarization_type == ConcentrationPolarizationType.fixed:
                return (prop_interface_io.conc_mass_phase_comp['Liq', j] ==
                        prop_io.conc_mass_phase_comp['Liq', j]
                        * self.cp_modulus[t, j])
            elif self.config.concentration_polarization_type == ConcentrationPolarizationType.calculated:
                jw = self.flux_mass_io_phase_comp[t, x, 'Liq', 'H2O'] / self.dens_solvent
                js = self.flux_mass_io_phase_comp[t, x, 'Liq', j]
                return (prop_interface_io.conc_mass_phase_comp['Liq', j] ==
                        prop_io.conc_mass_phase_comp['Liq', j] * exp(jw / self.Kf_io[t, x, j])
                        - js / jw * (exp(jw / self.Kf_io[t, x, j]) - 1))

        # Mass transfer coefficient calculation
        if self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated:
            @self.Constraint(self.flowsheet().config.time,
                                       self.io_list,
                                       solute_set,
                                       doc="Mass transfer coefficient in feed channel")
            # no "in" in 1DRO
            def eq_Kf_io(b, t, x, j):
                if x == 'in':
                    prop_io = b.feed_side.properties_in[t]
                elif x == 'out':
                    prop_io = b.feed_side.properties_out[t]
                return (b.Kf_io[t, x, j] * b.dh ==
                        prop_io.diffus_phase['Liq']
                        * b.N_Sh_io[t, x])

            @self.Constraint(self.flowsheet().config.time,
                                       self.io_list,
                                       doc="Sherwood number")
            def eq_N_Sh_io(b, t, x):
                return (b.N_Sh_io[t, x] ==
                        0.46 * (b.N_Re_io[t, x] * b.N_Sc_io[t, x])**0.36)

            @self.Constraint(self.flowsheet().config.time,
                                       self.io_list,
                                       doc="Schmidt number")
            def eq_N_Sc_io(b, t, x):
                if x == 'in':
                    prop_io = b.feed_side.properties_in[t]
                elif x == 'out':
                    prop_io = b.feed_side.properties_out[t]
                return (b.N_Sc_io[t, x] * prop_io.dens_mass_phase['Liq'] * prop_io.diffus_phase['Liq'] ==
                        prop_io.visc_d_phase['Liq'])

        if hasattr(self, 'length') or hasattr(self, 'width'):
            @self.Constraint(doc="Membrane area")
            def eq_area(b):
                return b.area == b.length * b.width

        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            or self.config.pressure_change_type == PressureChangeType.calculated):

            # var in 1DRO
            @self.Expression(doc="Cross-sectional area")
            def area_cross(b):
                return b.channel_height * b.width * b.spacer_porosity

            @self.Constraint(self.flowsheet().config.time,
                                       self.io_list,
                                       doc="Reynolds number")
            def eq_N_Re_io(b, t, x):
                if x == 'in':
                    prop_io = b.feed_side.properties_in[t]
                elif x == 'out':
                    prop_io = b.feed_side.properties_out[t]
                return (b.N_Re_io[t, x] * b.area_cross * prop_io.visc_d_phase['Liq'] ==
                        sum(prop_io.flow_mass_phase_comp['Liq', j] for j in b.config.property_package.component_list)
                        * b.dh)

        if self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length:
            # Pressure change equation when dP/dx = user-specified constant,
            @self.Constraint(self.flowsheet().config.time,
                             doc="pressure change due to friction")
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx[t] * b.length

        elif self.config.pressure_change_type == PressureChangeType.calculated:
            # Crossflow velocity at inlet and outlet
            @self.Constraint(self.flowsheet().config.time,
                             self.io_list,
                             doc="Crossflow velocity constraint")
            def eq_velocity_io(b, t, x):
                if x == 'in':
                    prop_io = b.feed_side.properties_in[t]
                elif x == 'out':
                    prop_io = b.feed_side.properties_out[t]
                return b.velocity_io[t, x] * b.area_cross == prop_io.flow_vol_phase['Liq']

            # Darcy friction factor based on eq. S27 in SI for Cost Optimization of Osmotically Assisted Reverse Osmosis
            # TODO: this relationship for friction factor is specific to a particular spacer geometry. Add alternatives.
            @self.Constraint(self.flowsheet().config.time,
                             self.io_list,
                             doc="Darcy friction factor constraint")
            def eq_friction_factor_darcy_io(b, t, x):
                return (b.friction_factor_darcy_io[t, x] - 0.42) * b.N_Re_io[t, x] == 189.3

            # Pressure change per unit length due to friction,
            # -1/2*f/dh*density*velocity^2
            @self.Constraint(self.flowsheet().config.time,
                             self.io_list,
                             doc="pressure change per unit length due to friction")
            def eq_dP_dx_io(b, t, x):
                if x == 'in':
                    prop_io = b.feed_side.properties_in[t]
                elif x == 'out':
                    prop_io = b.feed_side.properties_out[t]
                return (b.dP_dx_io[t, x] * b.dh ==
                        -0.5 * b.friction_factor_darcy_io[t, x]
                        * prop_io.dens_mass_phase['Liq'] * b.velocity_io[t, x]**2)

            # Average pressure change per unit length due to friction
            @self.Expression(self.flowsheet().config.time,
                             doc="expression for average pressure change per unit length due to friction")
            def dP_dx_avg(b, t):
                return 0.5 * sum(b.dP_dx_io[t, x] for x in b.io_list)

            # Pressure change equation
            @self.Constraint(self.flowsheet().config.time,
                             doc="pressure change due to friction")
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx_avg[t] * b.length

        # Bulk and interface connection on the feed-side
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Temperature at interface")
        def eq_equal_temp_interface_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
                prop_interface_io = b.properties_interface_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
                prop_interface_io = b.properties_interface_out[t]
            return prop_interface_io.temperature == \
                   prop_io.temperature

        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Pressure at interface")
        def eq_equal_pressure_interface_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
                prop_interface_io = b.properties_interface_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
                prop_interface_io = b.properties_interface_out[t]
            return prop_interface_io.pressure == \
                   prop_io.pressure

        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Volumetric flow at interface of inlet")
        def eq_equal_flow_vol_interface_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
                prop_interface_io = b.properties_interface_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
                prop_interface_io = b.properties_interface_out[t]
            return prop_interface_io.flow_vol_phase['Liq'] ==\
                   prop_io.flow_vol_phase['Liq']

    def _add_expressions(self):

        @self.Expression(self.flowsheet().config.time,
                         doc='Over pressure ratio')
        def over_pressure_ratio(b, t):
            return (b.feed_side.properties_out[t].pressure_osm
                    - b.permeate_side[t,'out'].pressure_osm) / \
                    b.feed_side.properties_out[t].pressure


    def initialize(blk,
                   initialize_guess=None,
                   state_args=None,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None,
                   fail_on_warning=False,
                   ignore_dof=False):
        """
        General wrapper for RO initialization routines

        Keyword Arguments:

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
                         property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : solver object or string indicating which solver to use during
                     initialization, if None provided the default solver will be used
                     (default = None)
            fail_on_warning : boolean argument to fail or only produce  warning upon unsuccessful solve (default=False)
            ignore_dof : boolean argument to ignore when DOF != 0 (default=False)
        Returns:
            None
        """

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver and options
        opt = get_solver(solver, optarg)


        # ---------------------------------------------------------------------
        # Extract initial state of inlet feed
        source = blk.feed_side.properties_in[blk.flowsheet().config.time.first()]
        state_args = blk._get_state_args(source, blk.mixed_permeate[0], initialize_guess, state_args)

        # Initialize feed inlet state block
        flags_feed_side = blk.feed_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['feed_side'],
            hold_state=True)

        init_log.info("Initialization Step 1 Complete.")
        if not ignore_dof:
            check_dof(blk, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # base properties on inlet state block

        blk.feed_side.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['retentate'],)
        blk.feed_side.properties_interface_in.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args['interface_in'],)
        blk.feed_side.properties_interface_out.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args['interface_out'],)
        blk.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['permeate'],)
        blk.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['permeate'],)
        init_log.info("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
            # occasionally it might be worth retrying a solve
            if not check_optimal_termination(res):
                init_log.warn("Trouble solving ReverseOsmosis0D unit model, trying one more time")
                res = opt.solve(blk, tee=slc.tee)
        check_solve(res, checkpoint='Initialization Step 3', logger=init_log, fail_flag=fail_on_warning)
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_side.release_state(flags_feed_side, outlvl)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[time_point, 'Liq']
        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[time_point, 'Liq', 'H2O']
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "length"):
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width"):
            var_dict["Membrane Width"] = self.width
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if hasattr(self, "N_Re_io"):
            var_dict["Reynolds Number @Inlet"] = self.N_Re_io[time_point, 'in']
            var_dict["Reynolds Number @Outlet"] = self.N_Re_io[time_point, 'out']
        if hasattr(self, "velocity_io"):
            var_dict["Velocity @Inlet"] = self.velocity_io[time_point, 'in']
            var_dict["Velocity @Outlet"] = self.velocity_io[time_point, 'out']
        for j in self.config.property_package.solute_set:
            cin_mem_name=f'{j} Concentration @Inlet,Membrane-Interface '
            var_dict[cin_mem_name] = (
                        self.feed_side.properties_interface_in[time_point].conc_mass_phase_comp['Liq', j])
            cout_mem_name=f'{j} Concentration @Outlet,Membrane-Interface '
            var_dict[cout_mem_name] = (
                self.feed_side.properties_interface_out[time_point].conc_mass_phase_comp['Liq', j])
            cin_bulk_name=f'{j} Concentration @Inlet,Bulk '
            var_dict[cin_bulk_name] = (
                        self.feed_side.properties_in[time_point].conc_mass_phase_comp['Liq', j])
            cout_bulk_name=f'{j} Concentration @Outlet,Bulk '
            var_dict[cout_bulk_name] = (
                self.feed_side.properties_out[time_point].conc_mass_phase_comp['Liq', j])
            cp_name=f'{j} Permeate Concentration '
            var_dict[cp_name] = (
                self.mixed_permeate[time_point].conc_mass_phase_comp['Liq', j])
        if self.feed_side.properties_interface_out[time_point].is_property_constructed('pressure_osm'):
            var_dict['Osmotic Pressure @Outlet,Membrane-Interface '] = (
                self.feed_side.properties_interface_out[time_point].pressure_osm)
        if self.feed_side.properties_out[time_point].is_property_constructed('pressure_osm'):
            var_dict['Osmotic Pressure @Outlet,Bulk'] = (
                self.feed_side.properties_out[time_point].pressure_osm)
        if self.feed_side.properties_interface_in[time_point].is_property_constructed('pressure_osm'):
            var_dict['Osmotic Pressure @Inlet,Membrane-Interface'] = (
                self.feed_side.properties_interface_in[time_point].pressure_osm)
        if self.feed_side.properties_in[time_point].is_property_constructed('pressure_osm'):
            var_dict['Osmotic Pressure @Inlet,Bulk'] = (
                self.feed_side.properties_in[time_point].pressure_osm)
        if self.feed_side.properties_in[time_point].is_property_constructed('flow_vol_phase'):
            var_dict['Volumetric Flowrate @Inlet'] = (
                self.feed_side.properties_in[time_point].flow_vol_phase['Liq'])
        if self.feed_side.properties_out[time_point].is_property_constructed('flow_vol_phase'):
            var_dict['Volumetric Flowrate @Outlet'] = (
                self.feed_side.properties_out[time_point].flow_vol_phase['Liq'])

        # TODO: (1) add more vars, (2) would be nice to add units to output, and (3) should be able to report output of
        #  "NaN" or "Not Reported", mainly for properties that exist but are not necessarily constructed within model
        #  constraints. One example in this case is the osmotic pressure of the bulk feed, which would certainly be of
        #  interest to users with a background in desalination (it is currently not reported because it is not directly
        #  used in any model constraints). Furthermore, the report() method seems to be limited to Pyomo Var objects for
        #  which the pyomo value() method is applied to. That is, a Pyomo Var object must be used; e.g., providing a
        #  list as output would yield an error.

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t,p,j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        for (t, x, p, j), v in self.flux_mass_io_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():  # scaling based on solvent flux equation
                    if x == 'in':
                        prop_io = self.feed_side.properties_in[t]
                    elif x == 'out':
                        prop_io = self.feed_side.properties_out[t]
                        prop_interface_io = self.feed_side.properties_interface_out[t]
                    sf = (iscale.get_scaling_factor(self.A_comp[t, j])
                          * iscale.get_scaling_factor(self.dens_solvent)
                          * iscale.get_scaling_factor(prop_io.pressure))
                    iscale.set_scaling_factor(v, sf)
                elif comp.is_solute():  # scaling based on solute flux equation
                    sf = (iscale.get_scaling_factor(self.B_comp[t, j])
                          * iscale.get_scaling_factor(self.feed_side.properties_in[t].conc_mass_phase_comp[p, j]))
                    iscale.set_scaling_factor(v, sf)

        for v in self.rejection_phase_comp.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'cp_modulus'):
            if iscale.get_scaling_factor(self.cp_modulus) is None:
                iscale.set_scaling_factor(self.cp_modulus, 1)

        if hasattr(self, 'Kf_io'):
            for v in self.Kf_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e4)

        if hasattr(self, 'N_Re_io'):
            for v in self.N_Re_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, 'N_Sc_io'):
            for v in self.N_Sc_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, 'N_Sh_io'):
            for v in self.N_Sh_io.values():
                if iscale.get_scaling_factor(v) is None:
                     iscale.set_scaling_factor(v, 1e-2)

        if hasattr(self, 'length'):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, 'width'):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)

        if hasattr(self, 'dP_dx'):
            for v in self.dP_dx.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, 'velocity_io'):
            for v in self.velocity_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'friction_factor_darcy_io'):
            for v in self.friction_factor_darcy_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'dP_dx_io'):
            for v in self.dP_dx_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
