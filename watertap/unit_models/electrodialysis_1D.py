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
                           Reals,
                           Reference,
                           value,
                           units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import Watertap cores
from watertap.core.util.initialization import check_solve, check_dof

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
from idaes.core.util.constants import Constants
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

        # Refer to the solute_set of the property package for a set of ions
        full_set = self.config.property_package.component_list
        if 'H2O' not in full_set:
            raise ConfigurationError("Property Package MUST constain 'H2O' as a component")

        # cation set reference
        cation_set = self.config.property_package.cation_set

        # anion set reference
        anion_set = self.config.property_package.anion_set

        # Create an ion_charge parameter to reference in model equations
        #   This is for convenience since the charge info is stored in
        #   different places between Adam's prop pack and the generic
        #   prop pack
        self.ion_charge = Param(
            anion_set | cation_set,
            initialize = 1,
            mutable = True,
            units=pyunits.dimensionless,
            doc="Ion charge",
        )

        # Loop through full set and try to assign charge
        for j in full_set:
            if j in anion_set or j in cation_set:
                self.ion_charge[j] = self.config.property_package.get_component(j).config.charge


        # Each dilute_side and concentrate_side pair has 2 membranes
        #   cem = Cation-Exchange Membrane
        #   aem = Anion-Exchange Membrane
        self.membrane_set = Set(initialize = ['cem','aem'])

        self.cell_width = Var(
            initialize = 0.1,
            bounds = (1e-3, 1e2),
            units = pyunits.m,
            doc = "The width of the electrodialysis cell, denoted as b in the model description"
        )

        self.membrane_thickness = Var(
            self.membrane_set,
            initialize = 0.0001,
            bounds = (1e-6, 1e-1),
            units = pyunits.m,
            doc = "Membrane thickness for each membrane"
        )

        self.ion_diffusivity_membrane = Var(
            self.membrane_set,
            anion_set | cation_set,
            initialize = 1e-8,
            bounds = (1e-20, 1e-3),
            units = pyunits.m**2 / pyunits.s,
            doc = "Ion diffusivity through each membrane"
        )

        self.water_permeability_membrane = Var(
            self.membrane_set,
            initialize = 5,
            bounds = (1e-20, 100),
            units = pyunits.m / pyunits.s / pyunits.Pa,
            doc = "Permeability for water through each membrane"
        )


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
        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

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

        """
            Add references to the shared control volume length.
            This is the same length for the full electrodialysis cell,
            and is denoted as 'l' in the model description.
        """
        add_object_reference(self, 'cell_length', self.dilute_side.length)

        # -------- Add constraints ---------
        # # TODO: Add vars and associated constraints for all flux terms
        #           There will be 1 flux var for water and 1 flux var for all ions
        #           (vars can be indexed by species, so we only write it once)
        #
        #           Those vars will be coupled into the mass_transfer_term below
        #           and will be of opposite sign for each channel

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

        # # TODO: Summate the flux terms for each mass transfer term in each domain

        # Non-electrical flux terms
        self.nonelec_flux = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0,
            domain=Reals,
            units=units_meta("amount") * units_meta("length") ** -2 * units_meta("time") ** -1,
            doc="Rate of flux of water and/or ions caused by osmotic pressure and/or diffusion",
        )
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc=" Constraint for rate of flux of water and/or ions caused by osmotic pressure and/or diffusion")
        def eq_nonelec_flux(self, t, x, p, j):
            if j == 'H2O':
                Posm_D = self.dilute_side.properties[t, x].pressure_osm_phase[p]
                Posm_C = self.concentrate_side.properties[t, x].pressure_osm_phase[p]
                L_cem = pyunits.convert( self.water_permeability_membrane['cem'],
                                to_units=units_meta("length") * units_meta("time") ** -1 * units_meta("pressure") ** -1 )
                L_aem = pyunits.convert( self.water_permeability_membrane['aem'],
                                to_units=units_meta("length") * units_meta("time") ** -1 * units_meta("pressure") ** -1 )
                return self.nonelec_flux[t, x, p, j] == (L_cem + L_aem) * \
                                (Posm_C*self.concentrate_side.properties[t,x].conc_mol_phase_comp[p,j] - \
                                Posm_D*self.dilute_side.properties[t,x].conc_mol_phase_comp[p,j])

            elif j in anion_set or j in cation_set:
                cem_rat = pyunits.convert( (self.ion_diffusivity_membrane['cem', j] / self.membrane_thickness['cem']),
                            to_units=units_meta("length") * units_meta("time") ** -1)
                aem_rat = pyunits.convert( (self.ion_diffusivity_membrane['aem', j] / self.membrane_thickness['aem']),
                            to_units=units_meta("length") * units_meta("time") ** -1)
                return self.nonelec_flux[t, x, p, j] == -(cem_rat + aem_rat) * \
                                                    (self.concentrate_side.properties[t, x].conc_mol_phase_comp[p, j] - \
                                                    self.dilute_side.properties[t, x].conc_mol_phase_comp[p, j])

            else:
                return self.nonelec_flux[t, x, p, j] == 0.0


        # Electrical flux terms
        self.elec_flux = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0,
            domain=Reals,
            units=units_meta("amount") * units_meta("length") ** -2 * units_meta("time") ** -1,
            doc="Rate of flux of water and/or ions caused by electrical current",
        )
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc=" Constraint for rate of flux of water and/or ions caused by electrical current")
        def eq_elec_flux(self, t, x, p, j):
            if j == 'H2O':
                return self.elec_flux[t, x, p, j] == 0.0
            elif j in anion_set or j in cation_set:
                return self.elec_flux[t, x, p, j] == 0.0
            else:
                return self.elec_flux[t, x, p, j] == 0.0

        # Add constraints for mass transfer terms (dilute_side)
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term on dilute side")
        def eq_mass_transfer_term_dilute(self, t, x, p, j):
            width = pyunits.convert( self.cell_width, to_units=units_meta("length"))
            return self.dilute_side.mass_transfer_term[t, x, p, j] == -self.nonelec_flux[t, x, p, j]*width \
                                                                     -self.elec_flux[t, x, p, j]*width

        # Add constraints for mass transfer terms (concentrate_side)
        @self.Constraint(self.flowsheet().config.time,
                         self.difference_elements,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term concentrate side")
        def eq_mass_transfer_term_concentrate(self, t, x, p, j):
            width = pyunits.convert( self.cell_width, to_units=units_meta("length"))
            return self.concentrate_side.mass_transfer_term[t, x, p, j] == self.nonelec_flux[t, x, p, j]*width \
                                                                        + self.elec_flux[t, x, p, j]*width


    # initialize method
    def initialize_build(
            blk,
            state_args=None,
            outlvl=idaeslog.NOTSET,
            solver=None,
            optarg=None,
            fail_on_warning=False,
            ignore_dof=False
    ):
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
            fail_on_warning : boolean argument to fail or only produce  warning upon unsuccessful solve (default=False)
            ignore_dof : boolean argument to ignore when DOF != 0 (default=False)

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
        if not ignore_dof:
            check_dof(blk, fail_flag=fail_on_warning, logger=init_log)

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
        check_solve(
            res,
            logger=init_log,
            fail_flag=fail_on_warning,
            checkpoint="Initialization Step 3",
        )

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

    def propogate_initial_state(self):

        i = 0
        for set in self.dilute_side.properties:

            # Should check 'define_state_vars' to see if user has provided
            #   state vars that are outside of the checks in this function
            if ('flow_mol_phase_comp' not in self.dilute_side.properties[set].define_state_vars()
                and 'flow_mass_phase_comp' not in self.dilute_side.properties[set].define_state_vars()):
                raise ConfigurationError("Electrodialysis1D unit model requires "
                                    "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                                    "state variable basis to apply the 'propogate_initial_state' method")

            t = set[0]
            x = set[1]
            if i == 0:
                intial_set = (t,x)

            if 'pressure' in self.dilute_side.properties[set].define_state_vars():
                self.dilute_side.properties[set].pressure = \
                    value(self.dilute_side.properties[intial_set].pressure)
                self.concentrate_side.properties[set].pressure = \
                    value(self.concentrate_side.properties[intial_set].pressure)

            if 'temperature' in self.dilute_side.properties[set].define_state_vars():
                self.dilute_side.properties[set].temperature = \
                    value(self.dilute_side.properties[intial_set].temperature)
                self.concentrate_side.properties[set].temperature = \
                    value(self.concentrate_side.properties[intial_set].temperature)

            if 'flow_mol_phase_comp' in self.dilute_side.properties[set].define_state_vars():
                for ind in self.dilute_side.properties[set].flow_mol_phase_comp:
                    self.dilute_side.properties[set].flow_mol_phase_comp[ind] = \
                        value(self.dilute_side.properties[intial_set].flow_mol_phase_comp[ind])
                    self.concentrate_side.properties[set].flow_mol_phase_comp[ind] = \
                        value(self.concentrate_side.properties[intial_set].flow_mol_phase_comp[ind])

                # Check to see if 'flow_mass_phase_comp' is constructed
                if self.dilute_side.properties[set].is_property_constructed("flow_mass_phase_comp"):
                    for ind in self.dilute_side.properties[set].flow_mass_phase_comp:
                        self.dilute_side.properties[set].flow_mass_phase_comp[ind] = \
                            value(self.dilute_side.properties[intial_set].flow_mol_phase_comp[ind])
                        self.concentrate_side.properties[set].flow_mass_phase_comp[ind] = \
                            value(self.concentrate_side.properties[intial_set].flow_mol_phase_comp[ind])

            if 'flow_mass_phase_comp' in self.dilute_side.properties[set].define_state_vars():
                for ind in self.dilute_side.properties[set].flow_mass_phase_comp:
                    self.dilute_side.properties[set].flow_mass_phase_comp[ind] = \
                        value(self.dilute_side.properties[intial_set].flow_mass_phase_comp[ind])
                    self.concentrate_side.properties[set].flow_mass_phase_comp[ind] = \
                        value(self.concentrate_side.properties[intial_set].flow_mass_phase_comp[ind])

                # Check to see if 'flow_mol_phase_comp' is constructed
                if self.dilute_side.properties[set].is_property_constructed("flow_mol_phase_comp"):
                    for ind in self.dilute_side.properties[set].flow_mol_phase_comp:
                        self.dilute_side.properties[set].flow_mol_phase_comp[ind] = \
                            value(self.dilute_side.properties[intial_set].flow_mass_phase_comp[ind])
                        self.concentrate_side.properties[set].flow_mol_phase_comp[ind] = \
                            value(self.concentrate_side.properties[intial_set].flow_mass_phase_comp[ind])

            i += 1

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # First, should set some values for all state vars based on inlet state
        #   We do this here because good scaling factors will depend on the
        #   starting values for the state vars
        self.propogate_initial_state()

        # Scaling factors that user may setup
        if iscale.get_scaling_factor(self.cell_width) is None:
            sf = iscale.get_scaling_factor(self.cell_width, default=1, warning=False)
            iscale.set_scaling_factor(self.cell_width, sf)

        if iscale.get_scaling_factor(self.membrane_thickness) is None:
            sf = iscale.get_scaling_factor(self.membrane_thickness, default=1e3, warning=False)
            iscale.set_scaling_factor(self.membrane_thickness, sf)

        if iscale.get_scaling_factor(self.ion_diffusivity_membrane) is None:
            sf = iscale.get_scaling_factor(self.ion_diffusivity_membrane, default=1e7, warning=False)
            iscale.set_scaling_factor(self.ion_diffusivity_membrane, sf)

        if iscale.get_scaling_factor(self.water_permeability_membrane) is None:
            sf = iscale.get_scaling_factor(self.water_permeability_membrane, default=1, warning=False)
            iscale.set_scaling_factor(self.water_permeability_membrane, sf)

        # Transform length constraint
        if ((iscale.get_scaling_factor(self.dilute_side.length) is None) and
            (iscale.get_scaling_factor(self.concentrate_side.length) is None)):
            sf = iscale.get_scaling_factor(self.dilute_side.length, default=1, warning=False)
        elif ((iscale.get_scaling_factor(self.dilute_side.length) is None) and
            (iscale.get_scaling_factor(self.concentrate_side.length) is not None)):
            sf = iscale.get_scaling_factor(self.concentrate_side.length, default=1, warning=False)
        else:
            sf = iscale.get_scaling_factor(self.dilute_side.length, default=1, warning=False)
        iscale.constraint_scaling_transform(self.eq_equal_length, sf)

        # Scaling factors for isothermal temperature
        if iscale.get_scaling_factor(self.dilute_side.properties[0, self.first_element].temperature) is None:
            sf = iscale.get_scaling_factor(self.dilute_side.properties[0, self.first_element].temperature, default=1e-2, warning=True)
        if iscale.get_scaling_factor(self.concentrate_side.properties[0, self.first_element].temperature) is None:
            sf = iscale.get_scaling_factor(self.concentrate_side.properties[0, self.first_element].temperature, default=1e-2, warning=True)

        sf = iscale.get_scaling_factor(self.concentrate_side.properties[0, self.first_element].temperature)
        for t in self.flowsheet().config.time:
            for x in self.difference_elements:
                iscale.constraint_scaling_transform(self.eq_isothermal_dilute[t,x], sf)
                iscale.constraint_scaling_transform(self.eq_isothermal_concentrate[t,x], sf)

        # # TODO: Add scaling factors


    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe({"Dilute Side Inlet": self.inlet_dilute,
                                              "Concentrate Side Inlet": self.inlet_concentrate,
                                              "Dilute Side Outlet": self.outlet_dilute,
                                              "Concentrate Side Outlet": self.outlet_concentrate},
                                                time_point=time_point)
