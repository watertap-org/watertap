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
                           check_optimal_termination,
                           Set,
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

    def _process_config(self):
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
        add_object_reference(self, 'length_domain', self.feed_side.length_domain)
        self.first_element = self.length_domain.first()
        self.difference_elements = Set(ordered=True, initialize=(x for x in self.length_domain if x != self.first_element))


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
            self.length_domain,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict)
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict)

        # Membrane interface: indexed state block
        self.feed_side.properties_interface = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.length_domain,
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
            add_object_reference(self, 'dP_dx', feed_side.deltaP)

        self._make_performance()

        self._add_expressions()

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

        # ==========================================================================

        self.width = Var(
            initialize=1,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta('length'),
            doc='Membrane width')

        super()._make_performance()

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, x, p, j):
            return value(self.feed_side.properties[t, x].get_material_flow_terms('Liq', j)
                         * self.recovery_mass_phase_comp[t, 'Liq', j])

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass') * units_meta('time')**-1 * units_meta('length')**-1,
            doc='Mass transfer to permeate')

        if self.config.has_pressure_change:
            self.deltaP = Var(
                self.flowsheet().config.time,
                initialize=-1e5,
                bounds=(-1e6, 0),
                domain=NegativeReals,
                units=units_meta('pressure'),
                doc='Pressure drop across unit')


        # ==========================================================================
        # Mass transfer term equation

        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(b, t, x, p, j):
            return b.mass_transfer_phase_comp[t, x, p, j] == -b.feed_side.mass_transfer_term[t, x, p, j]

        # ==========================================================================
        # Mass flux = feed mass transfer equation

        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_flux_equal_mass_transfer(b, t, x, p, j):
            return b.flux_mass_phase_comp[t, x, p, j] * b.width == -b.feed_side.mass_transfer_term[t, x, p, j]
        # ==========================================================================
        # Mass flux equations (Jw and Js)

        # ==========================================================================
        # Final permeate mass flow rate (of solvent and solute) --> Mp,j, final = sum(Mp,j)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Permeate mass flow rates exiting unit")
        def eq_permeate_production(b, t, p, j):
            return (b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == sum(b.permeate_side[t, x].get_material_flow_terms(p, j)
                           for x in b.difference_elements))
        # ==========================================================================
        # Feed and permeate-side mass transfer connection --> Mp,j = Mf,transfer = Jj * W * L/n

        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, x, p, j):
            return (b.permeate_side[t, x].get_material_flow_terms(p, j)
                    == -b.feed_side.mass_transfer_term[t, x, p, j] * b.length / b.nfe)

        ## ==========================================================================
        # Pressure drop
        if (self.config.pressure_change_type == PressureChangeType.fixed_per_unit_length
                or self.config.pressure_change_type == PressureChangeType.calculated):
            @self.Constraint(self.flowsheet().config.time,
                             doc='Pressure drop across unit')
            def eq_pressure_drop(b, t):
                return (b.deltaP[t] ==
                        sum(b.dP_dx[t, x] * b.length / b.nfe
                            for x in b.difference_elements))

        if (self.config.pressure_change_type == PressureChangeType.fixed_per_stage
                and self.config.has_pressure_change):
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc='Fixed pressure drop across unit')
            def eq_pressure_drop(b, t, x):
                return b.deltaP[t] == b.length * b.dP_dx[t, x]

        ## ==========================================================================
        # Feed-side isothermal conditions
        # NOTE: this could go on the feed_side block, but that seems to hurt initialization
        #       in the tests for this unit
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         doc="Isothermal assumption for feed channel")
        def eq_feed_isothermal(b, t, x):
            return b.feed_side.properties[t, b.first_element].temperature == \
                   b.feed_side.properties[t, x].temperature


    def initialize_build(blk,
                   initialize_guess=None,
                   state_args=None,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None,
                   fail_on_warning=False,
                   ignore_dof=False):
        """
        Initialization routine for 1D-RO unit.

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
        opt = get_solver(solver, optarg)

        source = blk.feed_side.properties[blk.flowsheet().config.time.first(), blk.first_element]
        state_args = blk._get_state_args(source, blk.mixed_permeate[0], initialize_guess, state_args)

        # ---------------------------------------------------------------------
        # Step 1: Initialize feed_side, permeate_side, and mixed_permeate blocks
        flags_feed_side = blk.feed_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['feed_side'],
            hold_state=True)

        init_log.info("Initialization Step 1 Complete")
        if not ignore_dof:
            check_dof(blk, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # base properties on inlet state block

        flag_feed_side_properties_interface = blk.feed_side.properties_interface.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args['interface'])
        flags_permeate_side = blk.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['permeate'])
        flags_mixed_permeate = blk.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args['permeate'])
        init_log.info("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
            # occasionally it might be worth retrying a solve
            if not check_optimal_termination(res):
                init_log.warn("Trouble solving ReverseOsmosis1D unit model, trying one more time")
                res = opt.solve(blk, tee=slc.tee)
        check_solve(res, logger=init_log, fail_flag=fail_on_warning, checkpoint='Initialization Step 3')
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_side.release_state(flags_feed_side, outlvl)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    def calculate_scaling_factors(self):
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties[0, 0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()

        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.width) is None:
            sf = iscale.get_scaling_factor(self.width, default=1, warning=True)
            iscale.set_scaling_factor(self.width, sf)

        if iscale.get_scaling_factor(self.length) is None:
            sf = iscale.get_scaling_factor(self.length, default=10, warning=True)
            iscale.set_scaling_factor(self.length, sf)

        # setting scaling factors for variables

        # will not override if the user provides the scaling factor
        ## default of 1 set by ControlVolume1D
        if iscale.get_scaling_factor(self.area_cross) == 1:
            iscale.set_scaling_factor(self.area_cross, 100)

        for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
            sf = (iscale.get_scaling_factor(self.feed_side.properties[t, x].get_material_flow_terms(p, j)) /
                  iscale.get_scaling_factor(self.feed_side.length)) * value(self.nfe)
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.feed_side.mass_transfer_term[t,x,p,j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, 'deltaP'):
            for v in self.deltaP.values():
                if iscale.get_scaling_factor(v) is None:
                     iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, 'dP_dx'):
            for v in self.feed_side.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e-5)
        else:
            for v in self.feed_side.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e5)
