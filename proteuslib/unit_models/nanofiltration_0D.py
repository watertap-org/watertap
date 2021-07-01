###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

# Import Pyomo libraries
from pyomo.environ import (Block,
                           Set,
                           Var,
                           Param,
                           Suffix,
                           NonNegativeReals,
                           Reference,
                           units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


@declare_process_block_class("NanoFiltration0D")
class NanoFiltrationData(UnitModelBlockData):
    """
    Standard NF Unit Model Class:
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
    **default** = False. NF units do not support dynamic
    behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. NF units do not have defined volume, thus
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
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        if (len(self.config.property_package.phase_list) > 1
                or 'Liq' not in [p for p in self.config.property_package.phase_list]):
            raise ConfigurationError(
                "NF model only supports one liquid phase ['Liq'],"
                "the property package has specified the following phases {}"
                    .format([p for p in self.config.property_package.phase_list]))

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # TODO: update IDAES such that solvent and solute lists are automatically created on the parameter block
        self.solvent_list = Set()
        self.solute_list = Set()
        for c in self.config.property_package.component_list:
            comp = self.config.property_package.get_component(c)
            try:
                if comp.is_solvent():
                    self.solvent_list.add(c)
                if comp.is_solute():
                    self.solute_list.add(c)
            except TypeError:
                raise ConfigurationError("NF model only supports one solvent and one or more solutes,"
                                         "the provided property package has specified a component '{}' "
                                         "that is not a solvent or solute".format(c))
        if len(self.solvent_list) > 1:
            raise ConfigurationError("NF model only supports one solvent component,"
                                     "the provided property package has specified {} solvent components"
                                     .format(len(self.solvent_list)))

        # Add unit parameters
        self.A_comp = Var(
            self.flowsheet().config.time,
            self.solvent_list,
            initialize=1e-12,
            bounds=(1e-18, 1e-6),
            domain=NonNegativeReals,
            units=units_meta('length')*units_meta('pressure')**-1*units_meta('time')**-1,
            doc='Solvent permeability coeff.')
        self.B_comp = Var(
            self.flowsheet().config.time,
            self.solute_list,
            initialize=1e-8,
            bounds=(1e-11, 1e-5),
            domain=NonNegativeReals,
            units=units_meta('length')*units_meta('time')**-1,
            doc='Solute permeability coeff.')
        self.sigma = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Reflection coefficient')
        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta('mass')*units_meta('length')**-3,
            doc='Pure water density')

        # Add unit variables
        self.flux_mass_phase_comp_in = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1e-3,
            bounds=(1e-12, 1e6),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Flux at feed inlet')
        self.flux_mass_phase_comp_out = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1e-3,
            bounds=(1e-12, 1e6),
            units=units_meta('mass')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Flux at feed outlet')
        self.avg_conc_mass_phase_comp_in = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.solute_list,
            initialize=1e-3,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass')*units_meta('length')**-3,
            doc='Average solute concentration at feed inlet')
        self.avg_conc_mass_phase_comp_out = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.solute_list,
            initialize=1e-3,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass')*units_meta('length')**-3,
            doc='Average solute concentration at feed outlet')
        self.area = Var(
            initialize=1,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('length') ** 2,
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
        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            bounds=(1e-8, 1e6),
            domain=NonNegativeReals,
            units=units_meta('mass')*units_meta('time')**-1,
            doc='Mass transfer to permeate')

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(self, t, p, j):
            return self.mass_transfer_phase_comp[t, p, j] == -self.feed_side.mass_transfer_term[t, p, j]

        # NF performance equations
        @self.Expression(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Average flux expression")
        def flux_mass_phase_comp_avg(b, t, p, j):
            return 0.5 * (b.flux_mass_phase_comp_in[t, p, j] + b.flux_mass_phase_comp_out[t, p, j])

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Permeate production")
        def eq_permeate_production(b, t, p, j):
            return (b.properties_permeate[t].get_material_flow_terms(p, j)
                    == b.area * b.flux_mass_phase_comp_avg[t, p, j])

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Inlet water and salt flux")
        def eq_flux_in(b, t, p, j):
            prop_feed = b.feed_side.properties_in[t]
            prop_perm = b.properties_permeate[t]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return (b.flux_mass_phase_comp_in[t, p, j] == b.A_comp[t, j] * b.dens_solvent
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - b.sigma[t] * (prop_feed.pressure_osm - prop_perm.pressure_osm)))
            elif comp.is_solute():
                return (b.flux_mass_phase_comp_in[t, p, j] == b.B_comp[t, j]
                        * (prop_feed.conc_mass_phase_comp[p, j] - prop_perm.conc_mass_phase_comp[p, j])
                        + ((1 - b.sigma[t]) * b.flux_mass_phase_comp_in[t, p, j]
                           * 1 / b.dens_solvent * b.avg_conc_mass_phase_comp_in[t, p, j])
                        )

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Outlet water and salt flux")
        def eq_flux_out(b, t, p, j):
            prop_feed = b.feed_side.properties_out[t]
            prop_perm = b.properties_permeate[t]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return (b.flux_mass_phase_comp_out[t, p, j] == b.A_comp[t, j] * b.dens_solvent
                        * ((prop_feed.pressure - prop_perm.pressure)
                           - b.sigma[t] * (prop_feed.pressure_osm - prop_perm.pressure_osm)))
            elif comp.is_solute():
                return (b.flux_mass_phase_comp_out[t, p, j] == b.B_comp[t, j]
                        * (prop_feed.conc_mass_phase_comp[p, j] - prop_perm.conc_mass_phase_comp[p, j])
                        + ((1 - b.sigma[t]) * b.flux_mass_phase_comp_out[t, p, j]
                           * 1 / b.dens_solvent * b.avg_conc_mass_phase_comp_out[t, p, j])
                        )

        # Average concentration
        # COMMENT: Chen approximation of logarithmic average implemented
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.solute_list,
                         doc="Average inlet concentration")
        def eq_avg_conc_in(b, t, p, j):
            prop_feed = b.feed_side.properties_in[t]
            prop_perm = b.properties_permeate[t]
            return (b.avg_conc_mass_phase_comp_in[t, p, j] ==
                    (prop_feed.conc_mass_phase_comp[p, j] * prop_perm.conc_mass_phase_comp[p, j] *
                     (prop_feed.conc_mass_phase_comp[p, j] + prop_perm.conc_mass_phase_comp[p, j])/2)**(1/3))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.solute_list,
                         doc="Average inlet concentration")
        def eq_avg_conc_out(b, t, p, j):
            prop_feed = b.feed_side.properties_out[t]
            prop_perm = b.properties_permeate[t]
            return (b.avg_conc_mass_phase_comp_out[t, p, j] ==
                    (prop_feed.conc_mass_phase_comp[p, j] * prop_perm.conc_mass_phase_comp[p, j] *
                     (prop_feed.conc_mass_phase_comp[p, j] + prop_perm.conc_mass_phase_comp[p, j])/2)**(1/3))

        # Feed and permeate-side connection
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer from feed to permeate")
        def eq_connect_mass_transfer(b, t, p, j):
            return (b.properties_permeate[t].get_material_flow_terms(p, j)
                    == -b.feed_side.mass_transfer_term[t, p, j])

        @self.Constraint(self.flowsheet().config.time,
                         doc="Enthalpy transfer from feed to permeate")
        def eq_connect_enthalpy_transfer(b, t):
            return (b.properties_permeate[t].get_enthalpy_flow_terms('Liq')
                    == -b.feed_side.enthalpy_transfer[t])

        @self.Constraint(self.flowsheet().config.time,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t):
            return b.feed_side.properties_out[t].temperature == \
                   b.properties_permeate[t].temperature

    def initialize(
            blk,
            state_args=None,
            outlvl=idaeslog.NOTSET,
            solver=None,
            optarg=None):
        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

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

    def get_costing(self, module=None, **kwargs):
        self.costing = Block()
        module.NanoFiltration_costing(self.costing, **kwargs)

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
        if iscale.get_scaling_factor(self.A_comp) is None:
            iscale.set_scaling_factor(self.A_comp, 1e11)

        if iscale.get_scaling_factor(self.B_comp) is None:
            iscale.set_scaling_factor(self.B_comp, 1e5)

        if iscale.get_scaling_factor(self.sigma) is None:
            iscale.set_scaling_factor(self.sigma, 1)

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].dens_mass_phase['Liq'])
            iscale.set_scaling_factor(self.dens_solvent, sf)

        for vobj in [self.flux_mass_phase_comp_in, self.flux_mass_phase_comp_out]:
            for (t, p, j), v in vobj.items():
                if iscale.get_scaling_factor(v) is None:
                    comp = self.config.property_package.get_component(j)
                    if comp.is_solvent():  # scaling based on solvent flux equation
                        sf = (iscale.get_scaling_factor(self.A_comp[t, j])
                              * iscale.get_scaling_factor(self.dens_solvent)
                              * iscale.get_scaling_factor(self.feed_side.properties_in[t].pressure))
                        iscale.set_scaling_factor(v, sf)
                    elif comp.is_solute():  # scaling based on solute flux equation
                        sf = (iscale.get_scaling_factor(self.B_comp[t, j])
                              * iscale.get_scaling_factor(self.feed_side.properties_in[t].conc_mass_phase_comp[p, j]))
                        iscale.set_scaling_factor(v, sf)

        for vobj in [self.avg_conc_mass_phase_comp_in, self.avg_conc_mass_phase_comp_out]:
            for (t, p, j), v in vobj.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].conc_mass_phase_comp[p, j])
                    iscale.set_scaling_factor(v, sf)

        for (t, p, j), v in self.feed_side.mass_transfer_term.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
                comp = self.config.property_package.get_component(j)
                if comp.is_solute:
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
                comp = self.config.property_package.get_component(j)
                if comp.is_solute:
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        # TODO: update IDAES control volume to scale mass_transfer and enthalpy_transfer
        for ind, v in self.feed_side.mass_transfer_term.items():
            (t, p, j) = ind
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.feed_side.mass_transfer_term[t, p, j])
                iscale.constraint_scaling_transform(self.feed_side.material_balances[t, j], sf)

        for t, v in self.feed_side.enthalpy_transfer.items():
            if iscale.get_scaling_factor(v) is None:
                sf = (iscale.get_scaling_factor(self.feed_side.properties_in[t].enth_flow))
                iscale.set_scaling_factor(v, sf)
                iscale.constraint_scaling_transform(self.feed_side.enthalpy_balances[t], sf)

        # transforming constraints
        for ind, c in self.eq_mass_transfer_term.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_permeate_production.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_flux_in.items():
            sf = iscale.get_scaling_factor(self.flux_mass_phase_comp_in[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_flux_out.items():
            sf = iscale.get_scaling_factor(self.flux_mass_phase_comp_out[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_avg_conc_in.items():
            sf = iscale.get_scaling_factor(self.avg_conc_mass_phase_comp_in[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_avg_conc_out.items():
            sf = iscale.get_scaling_factor(self.avg_conc_mass_phase_comp_out[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_connect_mass_transfer.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_connect_enthalpy_transfer.items():
            sf = iscale.get_scaling_factor(self.feed_side.enthalpy_transfer[ind])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_permeate_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].temperature)
            iscale.constraint_scaling_transform(c, sf)
