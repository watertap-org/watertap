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
# Import Pyomo libraries
from pyomo.environ import (Var,
                           Set,
                           NonNegativeReals,
                           NegativeReals,
                           Reference,
                           units as pyunits,
                           exp,
                           value)
from pyomo.common.collections import ComponentSet
# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import get_solver
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

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.io_list = Set(initialize=['in', 'out'])  # inlet/outlet set

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        # For permeate-specific scaling in calculate_scaling_factors
        self._permeate_scaled_properties = ComponentSet()

        # Add unit variables
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

        if ((self.config.pressure_change_type != PressureChangeType.fixed_per_stage)
                or (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)):
            self.length = Var(
                initialize=10,
                bounds=(0.1, 5e2),
                domain=NonNegativeReals,
                units=units_meta('length'),
                doc='Effective membrane length')
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
            self.dP_dx_io = Var(
                self.flowsheet().config.time,
                self.io_list,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta('pressure')*units_meta('length')**-1,
                doc="Pressure drop per unit length of feed channel at inlet and outlet")

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

        # Build control volume for permeate side
        self.permeate_side = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.permeate_side.add_state_blocks(
            has_phase_equilibrium=False)

        # Add additional state blocks
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets
        # Interface properties
        self.feed_side.properties_interface_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed-side interface at inlet",
            default=tmp_dict)
        self.feed_side.properties_interface_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed-side interface at outlet",
            default=tmp_dict)
        self.permeate_side.properties_mixed = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate",
            default=tmp_dict)

        # Add Ports
        self.add_inlet_port(name='inlet', block=self.feed_side)
        self.add_outlet_port(name='retentate', block=self.feed_side)
        self.add_port(name='permeate', block=self.permeate_side.properties_mixed)

        # References for control volume
        # pressure change
        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != 'none'):
            self.deltaP = Reference(self.feed_side.deltaP)

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

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(self, t, p, j):
            return self.mass_transfer_phase_comp[t, p, j] == -self.feed_side.mass_transfer_term[t, p, j]

        # RO performance equations
        @self.Expression(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Average flux expression")
        def flux_mass_phase_comp_avg(b, t, p, j):
            return 0.5 * sum(b.flux_mass_io_phase_comp[t, x, p, j] for x in self.io_list)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Permeate production")
        def eq_permeate_production(b, t, p, j):
            return (b.permeate_side.properties_mixed[t].get_material_flow_terms(p, j)
                    == b.area * b.flux_mass_phase_comp_avg[t, p, j])

        @self.Constraint(self.flowsheet().config.time,
                         self.io_list,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Water and salt flux")
        def eq_flux_io(b, t, x, p, j):
            if x == 'in':
                prop_feed = b.feed_side.properties_in[t]
                prop_feed_inter = b.feed_side.properties_interface_in[t]
                prop_perm = b.permeate_side.properties_in[t]
            elif x == 'out':
                prop_feed = b.feed_side.properties_out[t]
                prop_feed_inter = b.feed_side.properties_interface_out[t]
                prop_perm = b.permeate_side.properties_out[t]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return (b.flux_mass_io_phase_comp[t, x, p, j] == b.A_comp[t, j] * b.dens_solvent
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - (prop_feed_inter.pressure_osm - prop_perm.pressure_osm)))
            elif comp.is_solute():
                return (b.flux_mass_io_phase_comp[t, x, p, j] == b.B_comp[t, j]
                        * (prop_feed_inter.conc_mass_phase_comp[p, j] - prop_perm.conc_mass_phase_comp[p, j]))

        # Feed and permeate-side connection
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, p, j):
            return (b.permeate_side.properties_mixed[t].get_material_flow_terms(p, j)
                    == -b.feed_side.mass_transfer_term[t, p, j])

        @self.Constraint(self.flowsheet().config.time,
                         doc="Enthalpy transfer from feed to permeate")
        def eq_connect_enthalpy_transfer(b, t):
            return (b.permeate_side.properties_mixed[t].get_enthalpy_flow_terms('Liq')
                    == -b.feed_side.enthalpy_transfer[t])

        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t):
            return b.feed_side.properties_out[t].temperature == \
                   b.permeate_side.properties_mixed[t].temperature

        # # Permeate-side stateblocks
        @self.permeate_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   solute_set,
                                   doc="Permeate mass fraction")
        def eq_mass_frac_permeate_io(b, t, x, j):
            if x == 'in':
                prop_io = b.properties_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
            return (prop_io.mass_frac_phase_comp['Liq', j]
                    * sum(self.flux_mass_io_phase_comp[t, x, 'Liq', jj]
                          for jj in self.config.property_package.component_list)
                    == self.flux_mass_io_phase_comp[t, x, 'Liq', j])
        @self.permeate_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Permeate temperature")
        def eq_temperature_permeate_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
            return prop_io.temperature == b.properties_mixed[t].temperature

        @self.permeate_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Permeate pressure")
        def eq_pressure_permeate_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
            return prop_io.pressure == b.properties_mixed[t].pressure

        @self.permeate_side.Constraint(self.flowsheet().config.time,
                                   self.io_list,
                                   doc="Permeate flowrate")
        def eq_flow_vol_permeate_io(b, t, x):
            if x == 'in':
                prop_io = b.properties_in[t]
            elif x == 'out':
                prop_io = b.properties_out[t]
            return prop_io.flow_vol_phase['Liq'] == b.properties_mixed[t].flow_vol_phase['Liq']

        # Concentration polarization
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

        # constraints for additional variables (i.e. variables not used in other constraints)
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (b.recovery_vol_phase[t, 'Liq'] ==
                    b.permeate_side.properties_mixed[t].flow_vol_phase['Liq'] /
                    b.feed_side.properties_in[t].flow_vol_phase['Liq'])

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list)
        def eq_recovery_mass_phase_comp(b, t, j):
            return (b.recovery_mass_phase_comp[t, 'Liq', j] ==
                    b.permeate_side.properties_mixed[t].flow_mass_phase_comp['Liq', j] /
                    b.feed_side.properties_in[t].flow_mass_phase_comp['Liq', j])

        @self.Constraint(self.flowsheet().config.time,
                         solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return (b.rejection_phase_comp[t, 'Liq', j] ==
                    1 - (b.permeate_side.properties_mixed[t].conc_mass_phase_comp['Liq', j] /
                         b.feed_side.properties_in[t].conc_mass_phase_comp['Liq', j]))

        @self.Expression(self.flowsheet().config.time,
                         doc='Over pressure ratio')
        def over_pressure_ratio(b, t):
            return (b.feed_side.properties_out[t].pressure_osm
                    - b.permeate_side.properties_out[t].pressure_osm) / \
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
        if optarg is None:
            optarg = {'bound_push': 1e-8}
        opt = get_solver(solver, optarg)

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

        # ---------------------------------------------------------------------
        # Extract initial state of inlet feed
        if state_args is None:
            state_args = {}
            state_dict = blk.feed_side.properties_in[
                blk.flowsheet().config.time.first()].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value


        # Initialize feed inlet state block
        flags = blk.feed_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True)

        init_log.info_high("Initialization Step 1 Complete.")
        if not ignore_dof:
            check_dof(blk, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # base properties on inlet state block

        if 'flow_mass_phase_comp' not in state_args.keys():
            raise ConfigurationError('ReverseOsmosis0D initialization routine expects '
                                     'flow_mass_phase_comp as a state variable. Check '
                                     'that the property package supports this state '
                                     'variable or that the state_args provided to the '
                                     'initialize call includes this state variable')

        # slightly modify initial values for other state blocks
        state_args_retentate = deepcopy(state_args)
        state_args_permeate = deepcopy(state_args)

        state_args_retentate['pressure'] += initialize_guess['deltaP']
        state_args_permeate['pressure'] = blk.permeate_side.properties_mixed[0].pressure.value
        for j in blk.config.property_package.solvent_set:
            state_args_retentate['flow_mass_phase_comp'][('Liq', j)] *= (1 - initialize_guess['solvent_recovery'])
            state_args_permeate['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['solvent_recovery']
        for j in blk.config.property_package.solute_set:
            state_args_retentate['flow_mass_phase_comp'][('Liq', j)] *= (1 - initialize_guess['solute_recovery'])
            state_args_permeate['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['solute_recovery']

        state_args_interface_in = deepcopy(state_args)
        state_args_interface_out = deepcopy(state_args_retentate)

        for j in blk.config.property_package.solute_set:
            state_args_interface_in['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['cp_modulus']
            state_args_interface_out['flow_mass_phase_comp'][('Liq', j)] *= initialize_guess['cp_modulus']

        blk.feed_side.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_retentate,)
        blk.feed_side.properties_interface_in.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_interface_in,)
        blk.feed_side.properties_interface_out.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_interface_out,)
        blk.permeate_side.properties_mixed.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,)
        blk.permeate_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,)
        blk.permeate_side.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,)
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        check_solve(res, checkpoint='Initialization Step 3', logger=init_log, fail_flag=fail_on_warning)
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_side.release_state(flags, outlvl)
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
                self.permeate_side.properties_mixed[time_point].conc_mass_phase_comp['Liq', j])
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

    # permeate properties need to rescale solute values by 100
    def _rescale_permeate_variable(self, var, factor=100):
        if var not in self._permeate_scaled_properties:
            sf = iscale.get_scaling_factor(var)
            iscale.set_scaling_factor(var, sf * factor)
            self._permeate_scaled_properties.add(var)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for sb in (self.permeate_side.properties_in, self.permeate_side.properties_out,
                self.permeate_side.properties_mixed):
            for t in self.flowsheet().config.time:
                for j in self.config.property_package.solute_set:
                    self._rescale_permeate_variable(sb[t].flow_mass_phase_comp['Liq', j])
                    if sb[t].is_property_constructed('mass_frac_phase_comp'):
                        self._rescale_permeate_variable(sb[t].mass_frac_phase_comp['Liq', j])
                    if sb[t].is_property_constructed('conc_mass_phase_comp'):
                        self._rescale_permeate_variable(sb[t].conc_mass_phase_comp['Liq', j])
                    if sb[t].is_property_constructed('mole_frac_phase_comp'):
                        self._rescale_permeate_variable(sb[t].mole_frac_phase_comp[j])
                    if sb[t].is_property_constructed('molality_comp'):
                        self._rescale_permeate_variable(sb[t].molality_comp[j])
                if sb[t].is_property_constructed('pressure_osm'):
                    self._rescale_permeate_variable(sb[t].pressure_osm)

        # TODO: require users to set scaling factor for area or calculate it based on mass transfer and flux
        iscale.set_scaling_factor(self.area, 1e-1)

        # setting scaling factors for variables

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        for v in self.rejection_phase_comp.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, 'cp_modulus'):
            if iscale.get_scaling_factor(self.cp_modulus) is None:
                sf = iscale.get_scaling_factor(self.cp_modulus)
                iscale.set_scaling_factor(self.cp_modulus, sf)

        if hasattr(self, 'Kf_io'):
            for v in self.Kf_io.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e5)

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

        for (t, p, j), v in self.feed_side.mass_transfer_term.items():
            # already scaled by control volume with the default based on properties_in flow
            # solute typically has mass transfer 2 orders magnitude less than flow
            if j in self.config.property_package.solute_set:
                self._rescale_permeate_variable(v)

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
                comp = self.config.property_package.get_component(j)
                if comp.is_solute:
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        # transforming constraints
        for ind, c in self.eq_mass_transfer_term.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf/10.)

        for ind, c in self.eq_permeate_production.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf/10.)

        for ind, c in self.eq_flux_io.items():
            sf = iscale.get_scaling_factor(self.flux_mass_io_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_connect_mass_transfer.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf/10.)

        for ind, c in self.eq_connect_enthalpy_transfer.items():
            sf = iscale.get_scaling_factor(self.feed_side.enthalpy_transfer[ind])
            iscale.constraint_scaling_transform(c, sf*10.)

        for t, c in self.eq_permeate_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x, j), c in self.permeate_side.eq_mass_frac_permeate_io.items():
            if x == 'in':
                prop_io = self.permeate_side.properties_in[t]
            elif x == 'out':
                prop_io = self.permeate_side.properties_out[t]
            sf = iscale.get_scaling_factor(prop_io.mass_frac_phase_comp['Liq', j])
            iscale.constraint_scaling_transform(c, sf*100.)

        for (t, x), c in self.permeate_side.eq_temperature_permeate_io.items():
            if x == 'in':
                prop_io = self.permeate_side.properties_in[t]
            elif x == 'out':
                prop_io = self.permeate_side.properties_out[t]
            sf = iscale.get_scaling_factor(prop_io.temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.permeate_side.eq_pressure_permeate_io.items():
            if x == 'in':
                prop_io = self.permeate_side.properties_in[t]
            elif x == 'out':
                prop_io = self.permeate_side.properties_out[t]
            sf = iscale.get_scaling_factor(prop_io.pressure)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.permeate_side.eq_flow_vol_permeate_io.items():
            if x == 'in':
                prop_io = self.permeate_side.properties_in[t]
            elif x == 'out':
                prop_io = self.permeate_side.properties_out[t]
            sf = iscale.get_scaling_factor(prop_io.flow_vol_phase['Liq'])
            iscale.constraint_scaling_transform(c, sf)

        for (t, x, j), c in self.feed_side.eq_concentration_polarization_io.items():
            if x == 'in':
                prop_interface_io = self.feed_side.properties_interface_in[t]
            elif x == 'out':
                prop_interface_io = self.feed_side.properties_interface_out[t]
            sf = iscale.get_scaling_factor(prop_interface_io.conc_mass_phase_comp['Liq', j])
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_Kf_io'):
            for ind, c in self.eq_Kf_io.items():
                sf = iscale.get_scaling_factor(self.Kf_io[ind])
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_N_Re_io'):
            for ind, c in self.eq_N_Re_io.items():
                sf = iscale.get_scaling_factor(self.N_Re_io[ind])
                iscale.constraint_scaling_transform(c, sf*1e4)

        if hasattr(self, 'eq_N_Sc_io'):
            for ind, c in self.eq_N_Sc_io.items():
                sf = iscale.get_scaling_factor(self.N_Sc_io[ind])
                iscale.constraint_scaling_transform(c, sf*1e3)

        if hasattr(self, 'eq_N_Sh_io'):
            for ind, c in self.eq_N_Sh_io.items():
                sf = iscale.get_scaling_factor(self.N_Sh_io[ind])
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_area'):
            sf = iscale.get_scaling_factor(self.area)
            iscale.constraint_scaling_transform(self.eq_area, sf)

        if hasattr(self, 'eq_pressure_change'):
            for ind, c in self.eq_pressure_change.items():
                sf = iscale.get_scaling_factor(self.deltaP[ind])
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, 'eq_velocity_io'):
            for ind, c in self.eq_velocity_io.items():
                sf = iscale.get_scaling_factor(self.velocity_io[ind])
                iscale.constraint_scaling_transform(c, sf*100.)

        if hasattr(self, 'eq_friction_factor_darcy_io'):
            for ind, c in self.eq_friction_factor_darcy_io.items():
                sf = iscale.get_scaling_factor(self.friction_factor_darcy_io[ind])
                iscale.constraint_scaling_transform(c, sf/10.)

        if hasattr(self, 'eq_dP_dx_io'):
            for ind, c in self.eq_dP_dx_io.items():
                sf = iscale.get_scaling_factor(self.dP_dx_io[ind])
                iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.feed_side.eq_equal_temp_interface_io.items():
            if x == 'in':
                prop_interface_io = self.feed_side.properties_interface_in[t]
            elif x == 'out':
                prop_interface_io = self.feed_side.properties_interface_out[t]
            sf = iscale.get_scaling_factor(prop_interface_io.temperature)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.feed_side.eq_equal_pressure_interface_io.items():
            if x == 'in':
                prop_interface_io = self.feed_side.properties_interface_in[t]
            elif x == 'out':
                prop_interface_io = self.feed_side.properties_interface_out[t]
            sf = iscale.get_scaling_factor(prop_interface_io.pressure)
            iscale.constraint_scaling_transform(c, sf)

        for (t, x), c in self.feed_side.eq_equal_flow_vol_interface_io.items():
            if x == 'in':
                prop_interface_io = self.feed_side.properties_interface_in[t]
            elif x == 'out':
                prop_interface_io = self.feed_side.properties_interface_out[t]
            sf = iscale.get_scaling_factor(prop_interface_io.flow_vol_phase['Liq'])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_recovery_vol_phase.items():
            sf = iscale.get_scaling_factor(self.recovery_vol_phase[t, 'Liq'])
            iscale.constraint_scaling_transform(c, sf)

        for (t, j), c in self.eq_recovery_mass_phase_comp.items():
            sf = iscale.get_scaling_factor(self.recovery_mass_phase_comp[t, 'Liq', j])
            iscale.constraint_scaling_transform(c, sf)

        for (t, j), c in self.eq_rejection_phase_comp.items():
            sf = iscale.get_scaling_factor(self.rejection_phase_comp[t, 'Liq', j])
            iscale.constraint_scaling_transform(c, sf)
