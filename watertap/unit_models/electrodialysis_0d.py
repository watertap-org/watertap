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
from cmath import inf
from tkinter.messagebox import NO
from xmlrpc.client import Boolean
from attr import mutable
from numpy import integer
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Expression,
    Suffix,
    NonNegativeReals,
    Reference,
    value,
    Constraint,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
    components,
)
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from sympy import Domain, Integer, Integers
from idaes.core.util.constants import Constants

__author__ = "Austin Ladshaw, Xiangyu Bi"

_log = idaeslog.getLogger(__name__)

# Name of the unit model
@declare_process_block_class("Electrodialysis0D")
class Electrodialysis0DData(UnitModelBlockData):
    """
    0D Electrodialysis Model
    """

    # CONFIG are options for the unit model
    CONFIG = ConfigBlock()  #

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
        ),
    )

    CONFIG.declare(
        "operation_mode",
        ConfigValue(
            default="Constant Current",
            domain=In(["Constant Current", "Constant Voltage"]),
            description="The electrical operation mode. To be selected between Constant Current and Constant Voltage",
        ),
    )

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
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
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )

    # # TODO: For now, Adam's prop pack does not support and energy balance
    #           so we are making this none for now.
    # # TODO: Temporarily disabling energy balances
    # Should the commented part below be delected? -XB
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

    CONFIG.declare(
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

    CONFIG.declare(
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

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()
        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition

        # Create essential sets.
        ion_set = (
            self.config.property_package.ion_set
        )  # here it is assumed solute_set contains only ionic species.
        self.membrane_set = Set(initialize=["cem", "aem"])
        # Note: solute_set is used for now. Ion_set maybe used in the property package later.
        # Add unit variables and parameters

        # Create unit model parameters and vars
        self.water_density = Param(
            initialize=1000,
            mutable=False,
            units=pyunits.kg * pyunits.m**-3,
            doc="density of water",
        )
        self.water_MW = Param(
            initialize=18.015e-3,
            mutable=False,
            units=pyunits.kg * pyunits.mole**-1,
            doc="molecular weight of water",
        )

        # electrodialysis cell dimensianl properties
        self.cell_width = Var(
            initialize=0.1,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="The width of the electrodialysis cell, denoted as b in the model description",
        )
        self.cell_length = Var(
            initialize=0.5,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="The length of the electrodialysis cell, denoted as l in the model description",
        )
        self.spacer_thickness = Var(
            initialize=0.0001,
            units=pyunits.meter,
            doc="The distance between the concecutive aem and cem",
        )

        # Material and Operational properties
        self.membrane_thickness = Var(
            self.membrane_set,
            initialize=0.0001,
            bounds=(1e-6, 1e-1),
            units=pyunits.meter,
            doc="Membrane thickness",
        )
        self.ion_diffusivity_membrane = Var(
            self.membrane_set,
            ion_set,
            initialize=1e-8,
            bounds=(0, 1e-3),
            units=pyunits.meter**2 * pyunits.second**-1,
            doc="Ion diffusivity in the membrane phase",
        )
        self.ion_trans_number_membrane = Var(
            self.membrane_set,
            ion_set,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Ion transference number in the membrane phase",
        )
        self.water_trans_number_membrane = Var(
            self.membrane_set,
            initialize=5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Transference number of water in membranes",
        )
        self.water_permeability_membrane = Var(
            self.membrane_set,
            initialize=5,
            units=pyunits.meter * pyunits.second**-1 * pyunits.pascal**-1,
            doc="Water permeability coefficient",
        )
        self.membrane_surface_resistence = Var(
            self.membrane_set,
            initialize=2e-4,
            bounds=(1e-6, 1),
            units=pyunits.ohm * pyunits.meter**2,
            doc="Surface resistence of membrane",
        )
        self.current = Var(
            initialize=1,
            bounds=(0, 100),
            units=pyunits.amp,
            doc="Current across a cell-pair",
        )
        self.voltage = Var(
            initialize=10,
            units=pyunits.volt,
            doc="Voltage across a cell-pair, declared under the 'Constant Voltage' mode only",
        )
        self.current_utilization = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The current utilization efficiency",
        )
        # Fluxes Vars for constructing mass transfer terms
        self.elec_migration_flux_in = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by electrical migration",
        )
        self.elec_migration_flux_out = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_out of a component across the membrane driven by electrical migration",
        )
        self.nonelec_flux_in = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by non-electrical forces",
        )
        self.nonelec_flux_out = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_out of a component across the membrane driven by non-electrical forces",
        )

        # Build control volume for the dilute channel
        self.diluate_channel = ControlVolume0DBlock(
            default={
                "dynamic": False,
                "has_holdup": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )
        self.diluate_channel.add_state_blocks(has_phase_equilibrium=False)
        self.diluate_channel.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.diluate_channel.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )
        # # TODO: Consider adding energy balances

        # Build control volume for the concentrate channel
        self.concentrate_channel = ControlVolume0DBlock(
            default={
                "dynamic": False,
                "has_holdup": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )
        self.concentrate_channel.add_state_blocks(has_phase_equilibrium=False)
        self.concentrate_channel.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.concentrate_channel.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )
        # # TODO: Consider adding energy balances

        # Add ports (creates inlets and outlets for each channel)
        self.add_inlet_port(name="inlet_diluate", block=self.diluate_channel)
        self.add_outlet_port(name="outlet_diluate", block=self.diluate_channel)
        self.add_inlet_port(name="inlet_concentrate", block=self.concentrate_channel)
        self.add_outlet_port(name="outlet_concentrate", block=self.concentrate_channel)

        # Build Constraints
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            doc="Current-Voltage relationship",
        )
        def eq_current_voltage_relation(self, t, p):
            surface_resistance = (
                self.membrane_surface_resistence["aem"]
                + self.membrane_surface_resistence["cem"]
                + self.spacer_thickness
                / (
                    0.5
                    * (
                        self.concentrate_channel.properties_in[
                            t
                        ].electrical_conductivity_phase[p]
                        + self.concentrate_channel.properties_out[
                            t
                        ].electrical_conductivity_phase[p]
                        + self.diluate_channel.properties_in[
                            t
                        ].electrical_conductivity_phase[p]
                        + self.diluate_channel.properties_out[
                            t
                        ].electrical_conductivity_phase[p]
                    )
                )
            )
            return (
                self.current * surface_resistance
                == self.voltage * self.cell_width * self.cell_length
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for electrical migration flux_in",
        )
        def eq_elec_migration_flux_in(self, t, p, j):
            if j == "H2O":
                return self.elec_migration_flux_in[t, p, j] == (
                    self.water_trans_number_membrane["cem"]
                    + self.water_trans_number_membrane["aem"]
                ) * (
                    self.current
                    / (self.cell_width * self.cell_length)
                    / Constants.faraday_constant
                )
            elif self.config.property_package.ion_set:
                return self.elec_migration_flux_in[t, p, j] == (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.current
                    / (self.cell_width * self.cell_length)
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                return self.elec_migration_flux_out[t, p, j] == 0

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for electrical migration flux_out",
        )
        def eq_elec_migration_flux_out(self, t, p, j):
            if j == "H2O":
                return self.elec_migration_flux_out[t, p, j] == (
                    self.water_trans_number_membrane["cem"]
                    + self.water_trans_number_membrane["aem"]
                ) * (
                    self.current
                    / (self.cell_width * self.cell_length)
                    / Constants.faraday_constant
                )
            elif j in self.config.property_package.ion_set:
                return self.elec_migration_flux_out[t, p, j] == (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.current
                    / (self.cell_width * self.cell_length)
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                return self.elec_migration_flux_out[t, p, j] == 0

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux_in",
        )
        def eq_nonelec_flux_in(self, t, p, j):
            if j == "H2O":
                return self.nonelec_flux_in[
                    t, p, j
                ] == self.water_density / self.water_MW * (
                    self.water_permeability_membrane["cem"]
                    + self.water_permeability_membrane["aem"]
                ) * (
                    self.concentrate_channel.properties_in[t].pressure_osm_phase[p]
                    - self.diluate_channel.properties_in[t].pressure_osm_phase[p]
                )
            else:
                return self.nonelec_flux_in[t, p, j] == -(
                    self.ion_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.ion_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate_channel.properties_in[t].conc_mol_phase_comp[p, j]
                    - self.diluate_channel.properties_in[t].conc_mol_phase_comp[p, j]
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux_out",
        )
        def eq_nonelec_flux_out(self, t, p, j):
            if j == "H2O":
                return self.nonelec_flux_out[
                    t, p, j
                ] == self.water_density / self.water_MW * (
                    self.water_permeability_membrane["cem"]
                    + self.water_permeability_membrane["aem"]
                ) * (
                    self.concentrate_channel.properties_out[t].pressure_osm_phase[p]
                    - self.diluate_channel.properties_out[t].pressure_osm_phase[p]
                )
            else:
                return self.nonelec_flux_out[t, p, j] == -(
                    self.ion_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.ion_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate_channel.properties_out[t].conc_mol_phase_comp[p, j]
                    - self.diluate_channel.properties_out[t].conc_mol_phase_comp[p, j]
                )

        # Add constraints for mass transfer terms (diluate_channel)
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the diluate channel",
        )
        def eq_mass_transfer_term_diluate(self, t, p, j):
            return self.diluate_channel.mass_transfer_term[t, p, j] == -0.5 * (
                self.elec_migration_flux_in[t, p, j]
                + self.elec_migration_flux_out[t, p, j]
                + self.nonelec_flux_in[t, p, j]
                + self.nonelec_flux_out[t, p, j]
            ) * (self.cell_width * self.cell_length)

        # Add constraints for mass transfer terms (concentrate_channel)
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the concentrate channel",
        )
        def eq_mass_transfer_term_concentrate(self, t, p, j):
            return self.concentrate_channel.mass_transfer_term[t, p, j] == 0.5 * (
                self.elec_migration_flux_in[t, p, j]
                + self.elec_migration_flux_out[t, p, j]
                + self.nonelec_flux_in[t, p, j]
                + self.nonelec_flux_out[t, p, j]
            ) * (self.cell_width * self.cell_length)

        # Add isothermal condition
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal condition for the diluate channel",
        )
        def eq_isothermal_diluate(self, t):
            return (
                self.diluate_channel.properties_in[t].temperature
                == self.diluate_channel.properties_out[t].temperature
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal condition for the concentrate channel",
        )
        def eq_isothermal_concentrate(self, t):
            return (
                self.concentrate_channel.properties_in[t].temperature
                == self.concentrate_channel.properties_out[t].temperature
            )

    # initialize method
    def initialize(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize diluate_channel block
        flags_diluate = blk.diluate_channel.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize concentrate_side block
        flags_concentrate = blk.concentrate_channel.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,  # inlet var
            hold_state=True,
        )
        init_log.info_high("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release state
        blk.diluate_channel.release_state(flags_diluate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        blk.concentrate_channel.release_state(flags_concentrate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Diluate Channel Inlet": self.inlet_diluate,
                "Concentrate Channel Inlet": self.inlet_concentrate,
                "Diluate Channel Outlet": self.outlet_diluate,
                "Concentrate Channel Outlet": self.outlet_concentrate,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Var scaling
        # The following Vars' sf are allowed to be provided by users or set by default.
        if (
            iscale.get_scaling_factor(self.ion_diffusivity_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.ion_diffusivity_membrane, 1e9)
        if iscale.get_scaling_factor(self.membrane_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e8)
        if (
            iscale.get_scaling_factor(self.water_permeability_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.water_permeability_membrane, 1e14)
        if iscale.get_scaling_factor(self.cell_length, warning=True) is None:
            iscale.set_scaling_factor(self.cell_length, 1e1)
        if iscale.get_scaling_factor(self.cell_width, warning=True) is None:
            iscale.set_scaling_factor(self.cell_width, 1e1)
        if iscale.get_scaling_factor(self.spacer_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.spacer_thickness, 1e4)
        if (
            iscale.get_scaling_factor(self.membrane_surface_resistence, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.membrane_surface_resistence, 1e4)
        if iscale.get_scaling_factor(self.current, warning=True) is None:
            iscale.set_scaling_factor(self.current, 1)
        if iscale.get_scaling_factor(self.voltage, warning=True) is None:
            iscale.set_scaling_factor(self.voltage, 1)
        # The folloing Vars are built for constructing constraints and their sf are computed from other Vars.
        iscale.set_scaling_factor(
            self.elec_migration_flux_in,
            iscale.get_scaling_factor(self.current)
            * iscale.get_scaling_factor(self.cell_length) ** -1
            * iscale.get_scaling_factor(self.cell_width) ** -1
            * 1e5,
        )
        iscale.set_scaling_factor(
            self.elec_migration_flux_out,
            iscale.get_scaling_factor(self.current)
            * iscale.get_scaling_factor(self.cell_length) ** -1
            * iscale.get_scaling_factor(self.cell_width) ** -1
            * 1e5,
        )

        # Constraint scaling
        for ind, c in self.eq_current_voltage_relation.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.membrane_surface_resistence)
            )
        for ind, c in self.eq_elec_migration_flux_in.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.elec_migration_flux_in)
            )
        for ind, c in self.eq_elec_migration_flux_out.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.elec_migration_flux_out)
            )
        for ind, c in self.eq_nonelec_flux_in.items():
            if ind[2] == "H2O":
                sf = iscale.get_scaling_factor(
                    self.water_permeability_membrane
                ) * iscale.get_scaling_factor(
                    self.concentrate_channel.properties_in[ind[0]].pressure_osm_phase[
                        ind[1]
                    ]
                )
            sf = (
                iscale.get_scaling_factor(self.ion_diffusivity_membrane)
                / iscale.get_scaling_factor(self.membrane_thickness)
                * iscale.get_scaling_factor(
                    self.concentrate_channel.properties_in[ind[0]].conc_mol_phase_comp[
                        ind[1], ind[2]
                    ]
                )
            )
            iscale.set_scaling_factor(self.nonelec_flux_in[ind], sf)
            iscale.constraint_scaling_transform(c, sf)
        for ind, c in self.eq_nonelec_flux_out.items():
            if ind[2] == "H2O":
                sf = iscale.get_scaling_factor(
                    self.water_permeability_membrane
                ) * iscale.get_scaling_factor(
                    self.concentrate_channel.properties_out[ind[0]].pressure_osm_phase[
                        ind[1]
                    ]
                )
            sf = (
                iscale.get_scaling_factor(self.ion_diffusivity_membrane)
                / iscale.get_scaling_factor(self.membrane_thickness)
                * iscale.get_scaling_factor(
                    self.concentrate_channel.properties_out[ind[0]].conc_mol_phase_comp[
                        ind[1], ind[2]
                    ]
                )
            )
            iscale.set_scaling_factor(self.nonelec_flux_out[ind], sf)
            iscale.constraint_scaling_transform(c, sf)
        for ind, c in self.eq_mass_transfer_term_diluate.items():
            iscale.constraint_scaling_transform(
                c,
                min(
                    iscale.get_scaling_factor(self.elec_migration_flux_in[ind]),
                    iscale.get_scaling_factor(
                        self.nonelec_flux_in[ind], self.elec_migration_flux_out[ind]
                    ),
                    iscale.get_scaling_factor(self.nonelec_flux_out[ind]),
                ),
            )
        for ind, c in self.eq_mass_transfer_term_concentrate.items():
            iscale.constraint_scaling_transform(
                c,
                min(
                    iscale.get_scaling_factor(self.elec_migration_flux_in[ind]),
                    iscale.get_scaling_factor(
                        self.nonelec_flux_in[ind], self.elec_migration_flux_out[ind]
                    ),
                    iscale.get_scaling_factor(self.nonelec_flux_out[ind]),
                ),
            )
        for ind, c in self.eq_isothermal_diluate.items():
            iscale.constraint_scaling_transform(
                c, self.diluate_channel.properties_in[ind].temperature
            )
        for ind, c in self.eq_isothermal_concentrate.items():
            iscale.constraint_scaling_transform(
                c, self.concentrate_channel.properties_in[ind].temperature
            )
