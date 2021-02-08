##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

# Import Pyomo libraries
from pyomo.environ import Var, Constraint, Param, Expression, SolverFactory, \
    TerminationCondition, Suffix, NonNegativeReals, units as pyunits, Reference
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("RO_0D")
class ROData(UnitModelBlockData):
    """
    Standard RO Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """
    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. RO units do not support dynamic
    behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. RO units do not have defined volume, thus
    this must be False."""))
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
        doc="""Indicates what type of energy balance should be constructed,
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
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}"""))


    def build(self):
        # Call UnitModel.build to setup dynamics
        super(ROData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit parameters
        self.A = Var(
            self.flowsheet().config.time,
            initialize=1e-12,
            bounds=(1e-18, 1e-6),
            domain=NonNegativeReals,
            units=units_meta('length')*units_meta('pressure')**-1*units_meta('time')**-1,
            doc='Water permeability coeff.')
        self.B = Var(
            self.flowsheet().config.time,
            initialize=1e-8,
            bounds=(1e-11, 1e-5),
            domain=NonNegativeReals,
            units=units_meta('length')*units_meta('time')**-1,
            doc='Salt permeability coeff.')
        self.dens_H2O = Param(
            initialize=1000,
            units=units_meta('mass')*units_meta('length')**-3,
            doc='Pure water density')

        # Add unit variables
        self.flux_mass_comp_in = Var(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            initialize=1e-3,
            bounds=(1e-8, 1e6),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Flux at feed inlet')
        self.flux_mass_comp_out = Var(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            initialize=1e-3,
            # bounds=(1e-8, 1e6), # original
            bounds=(2.778e-4, 1e6),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Flux at feed outlet')
        self.area = Var(
            initialize=1,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('length')**2,
            doc='Membrane area')

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

        # Add permeate block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # permeate block is not an inlet
        self.properties_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of permeate",
            default=tmp_dict)

        # Add Ports
        self.add_inlet_port(name='inlet', block=self.feed_side)
        self.add_outlet_port(name='retentate', block=self.feed_side)
        self.add_port(name='permeate', block=self.properties_permeate)

        # References for control volume
        # pressure change
        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != 'none'):
            self.deltaP = Reference(self.feed_side.deltaP)

        # mass transfer
        # redundant, but useful to have on RO block, reference doesn't work with units at the moment
        # self.mass_transfer_comp = Reference(self.feed_side.mass_transfer_term[:, 'Liq', :])
        self.mass_transfer_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            initialize=1,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass') * units_meta('time')**-1,
            doc='Mass transfer to permeate')

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(b, t, j):
            return b.mass_transfer_comp[t, j] == -b.feed_side.mass_transfer_term[t, 'Liq', j]

        # RO performance equations
        @self.Expression(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Average flux expression")
        def flux_mass_comp_avg(b, t, j):
            return 0.5 * (b.flux_mass_comp_in[t, j] + b.flux_mass_comp_out[t, j])

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Permeate production")
        def eq_permeate_production(b, t, j):
            return b.properties_permeate[t].flow_mass_comp[j] == b.area * b.flux_mass_comp_avg[t, j]

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Inlet water and salt flux")
        def eq_flux_in(b, t, j):
            prop_feed = b.feed_side.properties_in[t]
            prop_perm = b.properties_permeate[t]
            if j == 'H2O':
                return (b.flux_mass_comp_in[t, j] == b.A[t] * b.dens_H2O
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - (prop_feed.pressure_osm - prop_perm.pressure_osm)))
            elif j == 'NaCl':
                return (b.flux_mass_comp_in[t, j] == b.B[t]
                        * (prop_feed.conc_mass_comp[j] - prop_perm.conc_mass_comp[j]))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Outlet water and salt flux")
        def eq_flux_out(b, t, j):
            prop_feed = b.feed_side.properties_out[t]
            prop_perm = b.properties_permeate[t]
            if j == 'H2O':
                return (b.flux_mass_comp_out[t, j] == b.A[t] * b.dens_H2O
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - (prop_feed.pressure_osm - prop_perm.pressure_osm)))
            elif j == 'NaCl':
                return (b.flux_mass_comp_out[t, j] == b.B[t]
                        * (prop_feed.conc_mass_comp[j] - prop_perm.conc_mass_comp[j]))

        # Feed and permeate-side connection
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, j):
            return b.properties_permeate[t].flow_mass_comp[j] == -b.feed_side.mass_transfer_term[t, 'Liq', j]

        @self.Constraint(self.flowsheet().config.time,
                         doc="Enthalpy transfer from feed to permeate")
        def eq_connect_enthalpy_transfer(b, t):
            return (sum(b.properties_permeate[t].flow_mass_comp[j]
                        for j in b.config.property_package.component_list)
                    * b.properties_permeate[t].enth_mass
                    == -b.feed_side.enthalpy_transfer[t])

        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t):
            return b.feed_side.properties_out[t].temperature == \
                   b.properties_permeate[t].temperature


    def initialize(
            blk,
            state_args=None,
            routine=None,
            outlvl=idaeslog.NOTSET,
            solver="ipopt",
            optarg={"tol": 1e-6}):
        """
        General wrapper for pressure changer initialization routines
        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - currently no specialized routine for RO unit
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = blk.feed_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize permeate
        # Set state_args from inlet state
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

        blk.properties_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_side.release_state(flags, outlvl + 1)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    def _get_performance_contents(self, time_point=0):
        # TODO: make a unit specific stream table
        var_dict = {}
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # TODO: require users to set scaling factor for area or calculate it based on mass transfer and flux
        iscale.set_scaling_factor(self.area, 1e-1)

        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.A) is None:
            iscale.set_scaling_factor(self.A, 1e12)

        if iscale.get_scaling_factor(self.B) is None:
            iscale.set_scaling_factor(self.B, 1e8)

        if iscale.get_scaling_factor(self.dens_H2O) is None:
            iscale.set_scaling_factor(self.dens_H2O, self.feed_side.properties_in[0].scaling_factor[
                        self.feed_side.properties_in[0].dens_mass])

        for vobj in [self.flux_mass_comp_in, self.flux_mass_comp_out]:
            for ind, v in vobj.items():
                t = ind[0]  # time
                j = ind[1]  # component
                if iscale.get_scaling_factor(v) is None:
                    if j == 'H2O':  # scaling based on water flux equation
                        sf = (iscale.get_scaling_factor(self.A[t])
                              * iscale.get_scaling_factor(self.dens_H2O)
                              * iscale.get_scaling_factor(self.feed_side.properties_in[t].pressure))
                        iscale.set_scaling_factor(v, sf)
                    elif j == 'NaCl':  # scaling based on salt flux equation
                        sf = (iscale.get_scaling_factor(self.B[t])
                              * iscale.get_scaling_factor(self.feed_side.properties_in[t].conc_mass_comp[j]))
                        iscale.set_scaling_factor(v, sf)

        for (t, p, j), v in self.feed_side.mass_transfer_term.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].flow_mass_comp[j])
                if j == 'NaCl':
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        for (t, j), v in self.mass_transfer_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].flow_mass_comp[j])
                if j == 'NaCl':
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        # TODO: update control volume to automatically provide scaling factors for mass_transfer and enthalpy_transfer
        for ind, v in self.feed_side.mass_transfer_term.items():
            (t, p, j) = ind
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.mass_transfer_term[t, p, j])
                iscale.constraint_scaling_transform(self.feed_side.material_balances[t, j], sf)

        for t, v in self.feed_side.enthalpy_transfer.items():
            if iscale.get_scaling_factor(v) is None:
                sf = (iscale.get_scaling_factor(self.feed_side.properties_in[t].flow_mass_comp['H2O'])
                      * iscale.get_scaling_factor(self.feed_side.properties_in[t].enth_mass))
                iscale.set_scaling_factor(v, sf)
                iscale.constraint_scaling_transform(self.feed_side.enthalpy_balances[t], sf)

        # transforming constraints
        for ind, c in self.eq_mass_transfer_term.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.mass_transfer_comp[ind]))

        for ind, c in self.eq_permeate_production.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.mass_transfer_comp[ind]))

        for ind, c in self.eq_flux_in.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.flux_mass_comp_in[ind]))

        for ind, c in self.eq_flux_out.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.flux_mass_comp_out[ind]))

        for ind, c in self.eq_connect_mass_transfer.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.mass_transfer_comp[ind]))

        for ind, c in self.eq_connect_enthalpy_transfer.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                self.feed_side.enthalpy_transfer[ind]))

        for t, c in self.eq_permeate_isothermal.items():
            iscale.constraint_scaling_transform(c,iscale.get_scaling_factor(
                self.feed_side.properties_in[t].temperature))

