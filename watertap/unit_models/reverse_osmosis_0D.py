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

        # for quacking like 1D model -> 0. is "in", 1. is "out"
        self.length_domain = Set(ordered=True, initialize=(0., 1.))  # inlet/outlet set
        add_object_reference(self, 'difference_elements', self.length_domain)
        self.first_element = self.length_domain.first()

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

        # for quacking like 1D model
        add_object_reference(self.feed_side,'properties',
                {**{(t,0.) : self.feed_side.properties_in[t] for t in self.flowsheet().config.time},
                 **{(t,1.) : self.feed_side.properties_out[t] for t in self.flowsheet().config.time}}
                )

        # Add additional state blocks
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets
        # Build permeate side
        self.permeate_side = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict)
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict)
        # Interface properties
        self.feed_side.properties_interface = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of feed-side membrane interface",
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

        units_meta = self.config.property_package.get_metadata().get_derived_units

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        if self.config.pressure_change_type == PressureChangeType.calculated:
            self.dP_dx = Var(
                self.flowsheet().config.time,
                self.length_domain,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta('pressure')*units_meta('length')**-1,
                doc="Pressure drop per unit length of feed channel at inlet and outlet")
        elif self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length:
            self.dP_dx = Var(
                self.flowsheet().config.time,
                initialize=-5e4,
                bounds=(-2e5, -1e3),
                domain=NegativeReals,
                units=units_meta('pressure')*units_meta('length')**-1,
                doc="pressure drop per unit length across feed channel")

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

        if (self.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            or self.config.pressure_change_type == PressureChangeType.calculated):
            self.area_cross = Var(
                    initialize=1e-3*1*0.95,
                    bounds=(0, 1e3),
                    domain=NonNegativeReals,
                    units=units_meta('length')**2,
                    doc='Cross sectional area')

        super()._make_performance()

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

        # constraints for additional variables (i.e. variables not used in other constraints)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(self, t, p, j):
            return self.mass_transfer_phase_comp[t, p, j] == -self.feed_side.mass_transfer_term[t, p, j]

        # Different expression in 1DRO
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
        # Not in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         solute_set,
                         doc="Permeate mass fraction")
        def eq_mass_frac_permeate(b, t, x, j):
            return (b.permeate_side[t, x].mass_frac_phase_comp['Liq', j]
                    * sum(self.flux_mass_phase_comp[t, x, 'Liq', jj]
                          for jj in self.config.property_package.component_list)
                    == self.flux_mass_phase_comp[t, x, 'Liq', j])

        # not in 1DRO
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Permeate flowrate")
        def eq_flow_vol_permeate(b, t, x):
            return b.permeate_side[t, x].flow_vol_phase['Liq'] == b.mixed_permeate[t].flow_vol_phase['Liq']

        if self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length:
            # Pressure change equation when dP/dx = user-specified constant,
            @self.Constraint(self.flowsheet().config.time,
                             doc="pressure change due to friction")
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx[t] * b.length

        elif self.config.pressure_change_type == PressureChangeType.calculated:

            # Average pressure change per unit length due to friction
            @self.Expression(self.flowsheet().config.time,
                             doc="expression for average pressure change per unit length due to friction")
            def dP_dx_avg(b, t):
                return 0.5 * sum(b.dP_dx[t, x] for x in b.length_domain)

            # Pressure change equation
            @self.Constraint(self.flowsheet().config.time,
                             doc="pressure change due to friction")
            def eq_pressure_change(b, t):
                return b.deltaP[t] == b.dP_dx_avg[t] * b.length


    def _add_expressions(self):
        super()._add_expressions()

        @self.Expression(self.flowsheet().config.time,
                         doc='Over pressure ratio')
        def over_pressure_ratio(b, t):
            return (b.feed_side.properties_out[t].pressure_osm
                    - b.permeate_side[t,1.].pressure_osm) / \
                    b.feed_side.properties_out[t].pressure

    def initialize_build(blk,
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
        blk.feed_side.properties_interface.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args['interface'],)
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

    def calculate_scaling_factors(self):
        # setting scaling factors for variables
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t,p,j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, 'area_cross'):
            if iscale.get_scaling_factor(self.area_cross) is None:
                iscale.set_scaling_factor(self.area_cross, 100)

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
