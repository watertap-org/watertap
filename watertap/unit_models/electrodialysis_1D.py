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
from pyomo.environ import (Block,
                           Set,
                           Var,
                           Param,
                           Expression,
                           Suffix,
                           NonNegativeReals,
                           Reference,
                           value,
                           units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume1DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault,
                        MaterialFlowBasis)
from idaes.core.control_volume1d import DistributedVars
from idaes.core.util.misc import add_object_reference
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Austin Ladshaw"

_log = idaeslog.getLogger(__name__)

# Name of the unit model
@declare_process_block_class("Electrodialysis1D")
class Electrodialysis1DData(UnitModelBlockData):
    """
    1D Electrodialysis Model
    """
    # CONFIG are options for the unit model
    CONFIG = ConfigBlock()

    # These config args are common to any control volume
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False."""))

    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
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

    # # TODO: For now, Adam's prop pack does not support and energy balance
    #           so we are making this none for now.
    # # TODO: Temporarily disabling energy balances
    '''
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.none,
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
    '''

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

    #These config args are specifically for 1D control volumes
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
            default="dae.finite_difference",
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
    documentation for supported transformations."""))

    CONFIG.declare("transformation_scheme", ConfigValue(
            default="BACKWARD",
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See
    Pyomo documentation for supported schemes."""))

    CONFIG.declare("finite_elements", ConfigValue(
            default=10,
            domain=int,
            description="Number of finite elements in length domain",
            doc="""Number of finite elements to use when discretizing length
            domain (default=10)"""))

    CONFIG.declare("collocation_points", ConfigValue(
            default=2,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=2)"""))


    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit variables and parameters
        # # TODO: Add material props for membranes and such here
        #   membrane area


        # Build control volume for dilute side
        self.dilute_side = ControlVolume1DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args,
            "area_definition": self.config.area_definition,
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points})

        self.dilute_side.add_geometry()
        self.dilute_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.dilute_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True)

        # # TODO: Temporarily disabling energy balances
        if hasattr(self.config, "energy_balance_type"):
            self.dilute_side.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False)

        self.dilute_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False)

        # Apply transformation to dilute_side
        self.dilute_side.apply_transformation()
        self.first_element = self.dilute_side.length_domain.first()
        self.difference_elements = Set(ordered=True, initialize=(x for x in self.dilute_side.length_domain if x != self.first_element))

        # Build control volume for concentrate side
        self.concentrate_side = ControlVolume1DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args,
            "area_definition": self.config.area_definition,
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points})

        self.concentrate_side.add_geometry()
        self.concentrate_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.concentrate_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True)

        # # TODO: Temporarily disabling energy balances
        if hasattr(self.config, "energy_balance_type"):
            self.concentrate_side.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False)

        self.concentrate_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False)

        # Apply transformation to concentrate_side
        self.concentrate_side.apply_transformation()

        # Add ports (creates inlets and outlets for each channel)
        self.add_inlet_port(name='inlet_dilute', block=self.dilute_side)
        self.add_outlet_port(name='outlet_dilute', block=self.dilute_side)

        self.add_inlet_port(name='inlet_concentrate', block=self.concentrate_side)
        self.add_outlet_port(name='outlet_concentrate', block=self.concentrate_side)

        """ Add references to the shared control volume length."""
        add_object_reference(self, 'length', self.dilute_side.length)

        # -------- Add constraints ---------
        # # TODO: Add vars and associated constraints for all flux terms
        #           There will be 1 flux var for water and 1 flux var for all ions
        #           (vars can be indexed by species, so we only write it once)
        #
        #           Those vars will be coupled into the mass_transfer_term below
        #           and will be of opposite sign for each channel




        # # TODO: Summate the flux terms for each mass transfer term in each domain

        # Adds isothermal constraint if no energy balance present
        if not hasattr(self.config, "energy_balance_type"):
            @self.Constraint(self.flowsheet().config.time,
                             self.difference_elements,
                             doc="Isothermal condition for dilute side")
            def eq_isothermal_dilute(self, t, x):
                return (self.dilute_side.properties[t, self.first_element].temperature == \
                        self.dilute_side.properties[t, x].temperature)

        if not hasattr(self.config, "energy_balance_type"):
            @self.Constraint(self.flowsheet().config.time,
                             self.difference_elements,
                             doc="Isothermal condition for concentrate side")
            def eq_isothermal_concentrate(self, t, x):
                return (self.concentrate_side.properties[t, self.first_element].temperature == \
                        self.concentrate_side.properties[t, x].temperature)

        # Add constraint for equal length of each channel
        @self.Constraint(doc="Constraint to ensure each channel has same length")
        def eq_equal_length(self):
            return self.dilute_side.length == self.concentrate_side.length

        # Add constraints for mass transfer terms (dilute_side)
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term on dilute side")
        def eq_mass_transfer_term_dilute(self, t, x, p, j):
            if j == 'H2O':
                return self.dilute_side.mass_transfer_term[t, x, p, j] == 0.0
            else:
                return self.dilute_side.mass_transfer_term[t, x, p, j] == 0.0

        # Add constraints for mass transfer terms (concentrate_side)
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term concentrate side")
        def eq_mass_transfer_term_concentrate(self, t, x, p, j):
            if j == 'H2O':
                return self.concentrate_side.mass_transfer_term[t, x, p, j] == 0.0
            else:
                return self.concentrate_side.mass_transfer_term[t, x, p, j] == 0.0


    # initialize method
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
        # Initialize dilute_side block
        flags = blk.dilute_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # Initialize concentrate_side block
        flags = blk.concentrate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release state
        blk.dilute_side.release_state(flags, outlvl + 1)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )
        blk.concentrate_side.release_state(flags, outlvl + 1)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # # TODO: Add scaling factors

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe({"Dilute Side Inlet": self.inlet_dilute,
                                              "Concentrate Side Inlet": self.inlet_concentrate,
                                              "Dilute Side Outlet": self.outlet_dilute,
                                              "Concentrate Side Outlet": self.outlet_concentrate},
                                                time_point=time_point)
