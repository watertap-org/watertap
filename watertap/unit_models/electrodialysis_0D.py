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
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Expression,
    Suffix,
    NonNegativeReals,
    NonNegativeIntegers,
    Reference,
    value,
    log,
    Constraint,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

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
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.constants import Constants
from enum import Enum

__author__ = " Xiangyu Bi, Austin Ladshaw,"

_log = idaeslog.getLogger(__name__)


class LimitingCurrentDensityMethod(Enum):
    InitialValue = 0
    # Empirical = 1
    # Theoretical = 2 TODO: 1 and 2


class ElectricalOperationMode(Enum):
    Constant_Current = 0
    Constant_Voltage = 1


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
            default=ElectricalOperationMode.Constant_Current,
            domain=In(ElectricalOperationMode),
            description="The electrical operation mode. To be selected between Constant Current and Constant Voltage",
        ),
    )

    CONFIG.declare(
        "limiting_current_density_method",
        ConfigValue(
            default=LimitingCurrentDensityMethod.InitialValue,
            domain=In(LimitingCurrentDensityMethod),
            description="Configuration for method to compute the limiting current density",
            doc="""
           **default** - ``LimitingCurrentDensityMethod.InitialValue``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``LimitingCurrentDensityMethod.InitialValue``", "Limiting current is calculated from a single initial value of the feed solution tested by the user."
           "``LimitingCurrentDensityMethod.Empirical``", "Limiting current density is caculated from the empirical equation: TODO"
           "``LimitingCurrentDensityMethod.Theoretical``", "Limiting current density is calculated from a theoretical equation: TODO"
       """,
        ),
    )

    CONFIG.declare(
        "limiting_current_density_data",
        ConfigValue(
            default=500,
            description="Limiting current density data input",
        ),
    )

    CONFIG.declare(
        "has_nonohmic_potential_membrane",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Configuration for whether to model the nonohmic potential across ion exchange membranes",
        ),
    )

    CONFIG.declare(
        "has_Nernst_diffusion_layer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Configuration for whether to simulate the concentration-polarized diffusion layers",
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

    # # TODO: Consider adding the EnergyBalanceType config using the following code
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

        # Create essential sets.
        self.membrane_set = Set(initialize=["cem", "aem"])
        self.electrode_side = Set(initialize=["cathode_left", "anode_right"])
        add_object_reference(self, "ion_set", self.config.property_package.ion_set)

        add_object_reference(
            self, "cation_set", self.config.property_package.cation_set
        )
        add_object_reference(self, "anion_set", self.config.property_package.anion_set)
        # Create unit model parameters and vars
        self.water_density = Param(
            initialize=1000,
            mutable=False,
            units=pyunits.kg * pyunits.m**-3,
            doc="density of water",
        )

        self.cell_pair_num = Var(
            initialize=1,
            domain=NonNegativeReals,
            bounds=(1, 10000),
            units=pyunits.dimensionless,
            doc="cell pair number in a stack",
        )

        # electrodialysis cell dimensional properties
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
            doc="The distance between the consecutive aem and cem",
        )

        # Material and Operational properties
        self.membrane_thickness = Var(
            self.membrane_set,
            initialize=0.0001,
            bounds=(1e-6, 1e-1),
            units=pyunits.meter,
            doc="Membrane thickness",
        )
        self.solute_diffusivity_membrane = Var(
            self.membrane_set,
            self.ion_set | self.config.property_package.solute_set,
            initialize=1e-10,
            bounds=(0.0, 1e-6),
            units=pyunits.meter**2 * pyunits.second**-1,
            doc="Solute (ionic and neutral) diffusivity in the membrane phase",
        )
        self.ion_trans_number_membrane = Var(
            self.membrane_set,
            self.ion_set,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Ion transference number in the membrane phase",
        )
        self.water_trans_number_membrane = Var(
            self.membrane_set,
            initialize=5,
            bounds=(0, 50),
            units=pyunits.dimensionless,
            doc="Transference number of water in membranes",
        )
        self.water_permeability_membrane = Var(
            self.membrane_set,
            initialize=1e-14,
            units=pyunits.meter * pyunits.second**-1 * pyunits.pascal**-1,
            doc="Water permeability coefficient",
        )
        self.membrane_areal_resistance = Var(
            self.membrane_set,
            initialize=2e-4,
            bounds=(1e-6, 1),
            units=pyunits.ohm * pyunits.meter**2,
            doc="Surface resistance of membrane",
        )
        self.electrodes_resistance = Var(
            initialize=0,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.ohm * pyunits.meter**2,
            doc="areal resistance of TWO electrode compartments of a stack",
        )
        self.current = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp,
            doc="Current across a cell-pair or stack",
        )
        self.voltage = Var(
            self.flowsheet().time,
            initialize=100,
            bounds=(0, 1000),
            units=pyunits.volt,
            doc="Voltage across a stack, declared under the 'Constant Voltage' mode only",
        )
        self.current_utilization = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The current utilization including water electro-osmosis and ion diffusion",
        )

        # Performance metrics
        self.current_efficiency = Var(
            self.flowsheet().time,
            initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The overall current efficiency for deionizaiton",
        )
        self.power_electrical = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, 12100),
            domain=NonNegativeReals,
            units=pyunits.watt,
            doc="Electrical power consumption of a stack",
        )
        self.specific_power_electrical = Var(
            self.flowsheet().time,
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.kW * pyunits.hour * pyunits.meter**-3,
            doc="Diluate-volume-flow-rate-specific electrical power consumption",
        )
        self.recovery_mass_H2O = Var(
            self.flowsheet().time,
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="water recovery ratio calculated by mass",
        )
        # TODO: consider adding more performance as needed.

        # Fluxes Vars for constructing mass transfer terms
        self.elec_migration_flux_in = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by electrical migration",
        )
        self.elec_migration_flux_out = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_out of a component across the membrane driven by electrical migration",
        )
        self.nonelec_flux_in = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by non-electrical forces",
        )
        self.nonelec_flux_out = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_out of a component across the membrane driven by non-electrical forces",
        )

        if (
            self.config.has_nonohmic_potential_membrane
            or self.config.has_Nernst_diffusion_layer
        ):
            self.conc_mem_surf_mol_ioa = Var(
                self.membrane_set,
                self.electrode_side,
                self.flowsheet().time,
                self.ion_set,
                initialize=500,
                bounds=(0, 1e5),
                units=pyunits.mol * pyunits.meter**-3,
                doc="Membane surface concentration of components",
            )

        # Build control volume for the dilute channel
        self.diluate = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.diluate.add_state_blocks(has_phase_equilibrium=False)
        self.diluate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.diluate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )
        # # TODO: Consider adding energy balances

        # Build control volume for the concentrate channel
        self.concentrate = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.concentrate.add_state_blocks(has_phase_equilibrium=False)
        self.concentrate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.concentrate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )
        # # TODO: Consider adding energy balances

        # Add ports (creates inlets and outlets for each channel)
        self.add_inlet_port(name="inlet_diluate", block=self.diluate)
        self.add_outlet_port(name="outlet_diluate", block=self.diluate)
        self.add_inlet_port(name="inlet_concentrate", block=self.concentrate)
        self.add_outlet_port(name="outlet_concentrate", block=self.concentrate)

        # extension options
        if self.config.has_nonohmic_potential_membrane:
            self._make_performance_nonohm_mem()
        if self.config.has_Nernst_diffusion_layer:
            self.current_dens_lim_ioa = Var(
                self.flowsheet().time,
                initialize=500,
                bounds=(0, 1000),
                units=pyunits.amp * pyunits.meter**-2,
                doc="Limiting Current Density accross the membrane as a function of the normalized length",
            )
            self._make_performance_dl_polarization()
        # Build Constraints
        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            doc="Current-Voltage relationship",
        )
        def eq_current_voltage_relation(self, t, p):
            if self.config.has_Nernst_diffusion_layer:
                total_areal_resistance = (
                    self.membrane_areal_resistance["aem"]
                    + self.membrane_areal_resistance["cem"]
                    + (
                        self.spacer_thickness
                        - self.dl_thickness_ioa["cem", "cathode_left", t]
                        - self.dl_thickness_ioa["aem", "anode_right", t]
                    )
                    * 0.5**-1
                    * (
                        self.concentrate.properties_in[t].elec_cond_phase["Liq"]
                        + self.concentrate.properties_out[t].elec_cond_phase["Liq"]
                    )
                    ** -1
                    + (
                        self.spacer_thickness
                        - self.dl_thickness_ioa["cem", "anode_right", t]
                        - self.dl_thickness_ioa["aem", "cathode_left", t]
                    )
                    * 0.5**-1
                    * (
                        self.diluate.properties_in[t].elec_cond_phase["Liq"]
                        + self.diluate.properties_out[t].elec_cond_phase["Liq"]
                    )
                    ** -1
                ) * self.cell_pair_num + self.electrodes_resistance
            else:
                total_areal_resistance = (
                    self.membrane_areal_resistance["aem"]
                    + self.membrane_areal_resistance["cem"]
                    + self.spacer_thickness
                    * (
                        0.5**-1
                        * (
                            self.concentrate.properties_in[t].elec_cond_phase["Liq"]
                            + self.concentrate.properties_out[t].elec_cond_phase["Liq"]
                        )
                        ** -1
                        + 0.5**-1
                        * (
                            self.diluate.properties_in[t].elec_cond_phase["Liq"]
                            + self.diluate.properties_out[t].elec_cond_phase["Liq"]
                        )
                        ** -1
                    )
                ) * self.cell_pair_num + self.electrodes_resistance  # the average conductivity of each channel's inlet and outlet is taken to represent that of the entire channel
            if self.config.has_nonohmic_potential_membrane:
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.current[t]
                        * (self.cell_width * self.cell_length) ** -1
                        * total_areal_resistance
                        + (
                            self.potential_ohm_dl_ioa["cem", t]
                            + self.potential_ohm_dl_ioa["aem", t]
                            + self.potential_nonohm_dl_ioa["cem", t]
                            + self.potential_nonohm_dl_ioa["aem", t]
                            + self.potential_nonohm_membrane_ioa["cem", t]
                            + self.potential_nonohm_membrane_ioa["aem", t]
                        )
                        * self.cell_pair_num
                        == self.voltage[t]
                    )
                else:
                    return (
                        self.current[t]
                        * (self.cell_width * self.cell_length) ** -1
                        * total_areal_resistance
                        + (
                            +self.potential_nonohm_membrane_ioa["cem", t]
                            + self.potential_nonohm_membrane_ioa["aem", t]
                        )
                        * self.cell_pair_num
                        == self.voltage[t]
                    )
            else:
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.current[t]
                        * (self.cell_width * self.cell_length) ** -1
                        * total_areal_resistance
                        + (
                            self.potential_ohm_dl_ioa["cem", t]
                            + self.potential_ohm_dl_ioa["aem", t]
                            + self.potential_nonohm_dl_ioa["cem", t]
                            + self.potential_nonohm_dl_ioa["aem", t]
                        )
                        * self.cell_pair_num
                        == self.voltage[t]
                    )
                else:
                    return (
                        self.current[t]
                        * (self.cell_width * self.cell_length) ** -1
                        * total_areal_resistance
                        == self.voltage[t]
                    )

        @self.Constraint(
            self.flowsheet().time,
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
                    self.current[t]
                    / (self.cell_width * self.cell_length)
                    / Constants.faraday_constant
                )
            elif j in self.ion_set:
                return self.elec_migration_flux_in[t, p, j] == (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.current[t]
                    / (self.cell_width * self.cell_length)
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                self.elec_migration_flux_in[t, p, j].fix(0)
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
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
                    self.current[t]
                    / (self.cell_width * self.cell_length)
                    / Constants.faraday_constant
                )
            elif j in self.ion_set:
                return self.elec_migration_flux_out[t, p, j] == (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.current[t]
                    / (self.cell_width * self.cell_length)
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                self.elec_migration_flux_out[t, p, j].fix(0)
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux_in",
        )
        def eq_nonelec_flux_in(self, t, p, j):
            if j == "H2O":
                if self.config.has_Nernst_diffusion_layer:
                    return self.nonelec_flux_in[
                        t, p, j
                    ] == self.water_density / self.config.property_package.mw_comp[
                        j
                    ] * (
                        self.water_permeability_membrane["cem"]
                        + self.water_permeability_membrane["aem"]
                    ) * (
                        self.concentrate.properties_in[t].pressure_osm_phase[p]
                        * (
                            1
                            + self.current[t]
                            * (self.cell_width * self.cell_length) ** -1
                            / self.current_dens_lim_ioa[t]
                        )
                        - self.diluate.properties_in[t].pressure_osm_phase[p]
                        * (
                            1
                            - self.current[t]
                            * (self.cell_width * self.cell_length) ** -1
                            / self.current_dens_lim_ioa[t]
                        )
                    )
                else:
                    return self.nonelec_flux_in[
                        t, p, j
                    ] == self.water_density / self.config.property_package.mw_comp[
                        j
                    ] * (
                        self.water_permeability_membrane["cem"]
                        + self.water_permeability_membrane["aem"]
                    ) * (
                        self.concentrate.properties_in[t].pressure_osm_phase[p]
                        - self.diluate.properties_in[t].pressure_osm_phase[p]
                    )

            else:
                if self.config.has_Nernst_diffusion_layer:
                    return self.nonelec_flux_in[t, p, j] == -(
                        self.solute_diffusivity_membrane["cem", j]
                        * self.membrane_thickness["cem"] ** -1
                        * (
                            self.conc_mem_surf_mol_ioa["cem", "cathode_left", t, j]
                            - self.conc_mem_surf_mol_ioa["cem", "anode_right", t, j]
                        )
                        + self.solute_diffusivity_membrane["aem", j]
                        * self.membrane_thickness["aem"] ** -1
                        * (
                            self.conc_mem_surf_mol_ioa["aem", "anode_right", t, j]
                            - self.conc_mem_surf_mol_ioa["aem", "cathode_left", t, j]
                        )
                    )

                else:
                    return self.nonelec_flux_in[t, p, j] == -(
                        self.solute_diffusivity_membrane["cem", j]
                        / self.membrane_thickness["cem"]
                        + self.solute_diffusivity_membrane["aem", j]
                        / self.membrane_thickness["aem"]
                    ) * (
                        self.concentrate.properties_in[t].conc_mol_phase_comp[p, j]
                        - self.diluate.properties_in[t].conc_mol_phase_comp[p, j]
                    )

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux_out",
        )
        def eq_nonelec_flux_out(self, t, p, j):
            if j == "H2O":
                if self.config.has_Nernst_diffusion_layer:
                    return self.nonelec_flux_out[
                        t, p, j
                    ] == self.water_density / self.config.property_package.mw_comp[
                        j
                    ] * (
                        self.water_permeability_membrane["cem"]
                        + self.water_permeability_membrane["aem"]
                    ) * (
                        self.concentrate.properties_out[t].pressure_osm_phase[p]
                        * (
                            1
                            + self.current[t]
                            * (self.cell_width * self.cell_length) ** -1
                            / self.current_dens_lim_ioa[t]
                        )
                        - self.diluate.properties_out[t].pressure_osm_phase[p]
                        * (
                            1
                            - self.current[t]
                            * (self.cell_width * self.cell_length) ** -1
                            / self.current_dens_lim_ioa[t]
                        )
                    )
                else:
                    return self.nonelec_flux_out[
                        t, p, j
                    ] == self.water_density / self.config.property_package.mw_comp[
                        j
                    ] * (
                        self.water_permeability_membrane["cem"]
                        + self.water_permeability_membrane["aem"]
                    ) * (
                        self.concentrate.properties_out[t].pressure_osm_phase[p]
                        - self.diluate.properties_out[t].pressure_osm_phase[p]
                    )

            else:
                if self.config.has_Nernst_diffusion_layer:
                    return self.nonelec_flux_out[t, p, j] == -(
                        self.solute_diffusivity_membrane["cem", j]
                        * self.membrane_thickness["cem"] ** -1
                        * (
                            self.conc_mem_surf_mol_ioa["cem", "cathode_left", t, j]
                            - self.conc_mem_surf_mol_ioa["cem", "anode_right", t, j]
                        )
                        + self.solute_diffusivity_membrane["aem", j]
                        * self.membrane_thickness["aem"] ** -1
                        * (
                            self.conc_mem_surf_mol_ioa["aem", "anode_right", t, j]
                            - self.conc_mem_surf_mol_ioa["aem", "cathode_left", t, j]
                        )
                    )

                else:
                    return self.nonelec_flux_out[t, p, j] == -(
                        self.solute_diffusivity_membrane["cem", j]
                        / self.membrane_thickness["cem"]
                        + self.solute_diffusivity_membrane["aem", j]
                        / self.membrane_thickness["aem"]
                    ) * (
                        self.concentrate.properties_out[t].conc_mol_phase_comp[p, j]
                        - self.diluate.properties_out[t].conc_mol_phase_comp[p, j]
                    )

        # Add constraints for mass transfer terms (diluate)
        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the diluate channel",
        )
        def eq_mass_transfer_term_diluate(self, t, p, j):
            return (
                self.diluate.mass_transfer_term[t, p, j]
                == -0.5
                * (
                    self.elec_migration_flux_in[t, p, j]
                    + self.elec_migration_flux_out[t, p, j]
                    + self.nonelec_flux_in[t, p, j]
                    + self.nonelec_flux_out[t, p, j]
                )
                * (self.cell_width * self.cell_length)
                * self.cell_pair_num
            )

        # Add constraints for mass transfer terms (concentrate)
        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the concentrate channel",
        )
        def eq_mass_transfer_term_concentrate(self, t, p, j):
            return (
                self.concentrate.mass_transfer_term[t, p, j]
                == 0.5
                * (
                    self.elec_migration_flux_in[t, p, j]
                    + self.elec_migration_flux_out[t, p, j]
                    + self.nonelec_flux_in[t, p, j]
                    + self.nonelec_flux_out[t, p, j]
                )
                * (self.cell_width * self.cell_length)
                * self.cell_pair_num
            )

        # Add isothermal condition
        @self.Constraint(
            self.flowsheet().time,
            doc="Isothermal condition for the diluate channel",
        )
        def eq_isothermal_diluate(self, t):
            return (
                self.diluate.properties_in[t].temperature
                == self.diluate.properties_out[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Isothermal condition for the concentrate channel",
        )
        def eq_isothermal_concentrate(self, t):
            return (
                self.concentrate.properties_in[t].temperature
                == self.concentrate.properties_out[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Electrical power consumption of a stack",
        )
        def eq_power_electrical(self, t):
            return self.power_electrical[t] == self.current[t] * self.voltage[t]

        @self.Constraint(
            self.flowsheet().time,
            doc="Diluate_volume_flow_rate_specific electrical power consumption of a stack",
        )
        def eq_specific_power_electrical(self, t):
            return (
                pyunits.convert(
                    self.specific_power_electrical[t],
                    pyunits.watt * pyunits.second * pyunits.meter**-3,
                )
                * self.diluate.properties_out[t].flow_vol_phase["Liq"]
                == self.current[t] * self.voltage[t]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Overall current efficiency evaluation",
        )
        def eq_current_efficiency(self, t):
            return (
                self.current_efficiency[t] * self.current[t] * self.cell_pair_num
                == sum(
                    self.diluate.properties_in[t].flow_mol_phase_comp["Liq", j]
                    * self.config.property_package.charge_comp[j]
                    - self.diluate.properties_out[t].flow_mol_phase_comp["Liq", j]
                    * self.config.property_package.charge_comp[j]
                    for j in self.config.property_package.cation_set
                )
                * Constants.faraday_constant
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Water recovery by mass",
        )
        def eq_recovery_mass_H2O(self, t):
            return (
                self.recovery_mass_H2O[t]
                * (
                    self.diluate.properties_in[t].flow_mass_phase_comp["Liq", "H2O"]
                    + self.concentrate.properties_in[t].flow_mass_phase_comp[
                        "Liq", "H2O"
                    ]
                )
                == self.diluate.properties_out[t].flow_mass_phase_comp["Liq", "H2O"]
            )

    def _make_performance_nonohm_mem(self):

        self.potential_nonohm_membrane_ioa = Var(
            self.membrane_set,
            self.flowsheet().time,
            initialize=0.01,  # to reinspect
            bounds=(-50, 50),
            units=pyunits.volt,
            doc="Nonohmic potential across a membane",
        )
        # ioa = in-out average

        @self.Constraint(
            self.membrane_set,
            self.electrode_side,
            self.flowsheet().time,
            self.ion_set,
            doc="calcualte current density from the electrical input",
        )
        def eq_set_surface_conc_ioa(self, mem, side, t, j):
            if not self.config.has_Nernst_diffusion_layer:
                if (mem == "cem" and side == "cathode_left") or (
                    mem == "aem" and side == "anode_right"
                ):
                    return self.conc_mem_surf_mol_ioa[mem, side, t, j] == 0.5 * (
                        self.concentrate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        + self.concentrate.properties_out[t].conc_mol_phase_comp[
                            "Liq", j
                        ]
                    )
                else:
                    return self.conc_mem_surf_mol_ioa[mem, side, t, j] == 0.5 * (
                        self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        + self.diluate.properties_out[t].conc_mol_phase_comp["Liq", j]
                    )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            doc="Calculate the total non-ohmic potential across an iem; this takes account of diffusion and Donnan Potentials",
        )
        def eq_potential_nonohm_membrane_ioa(self, mem, t):

            return self.potential_nonohm_membrane_ioa[mem, t] == (
                Constants.gas_constant
                * self.diluate.properties_in[t].temperature
                / Constants.faraday_constant
                * (
                    sum(
                        self.ion_trans_number_membrane[mem, j]
                        / self.config.property_package.charge_comp[j]
                        for j in self.cation_set
                    )
                    + sum(
                        self.ion_trans_number_membrane[mem, j]
                        / self.config.property_package.charge_comp[j]
                        for j in self.anion_set
                    )
                )
                * log(
                    sum(
                        self.conc_mem_surf_mol_ioa[mem, "cathode_left", t, j]
                        for j in self.ion_set
                    )
                    / sum(
                        self.conc_mem_surf_mol_ioa[mem, "anode_right", t, j]
                        for j in self.ion_set
                    )
                )
            )

    def _make_performance_dl_polarization(self):

        self.potential_nonohm_dl_ioa = Var(
            self.membrane_set,
            self.flowsheet().time,
            initialize=0,
            bounds=(-50, 50),
            units=pyunits.volt,
            doc="Nonohmic potential in two diffusion layers on the two sides of a membrane",
        )

        self.potential_ohm_dl_ioa = Var(
            self.membrane_set,
            self.flowsheet().time,
            initialize=0.001,
            bounds=(0, 50),
            units=pyunits.volt,
            doc="Ohmic potential in two diffusion layers on the two sides of a membrane",
        )

        self.dl_thickness_ioa = Var(
            self.membrane_set,
            self.electrode_side,
            self.flowsheet().time,
            initialize=0.0005,
            bounds=(0, 1e-2),
            units=pyunits.m,
            doc="Thickness of the diffusion layer",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Calculate length-indexed limitting current density",
        )
        def eq_current_dens_lim_ioa(self, t):
            if (
                self.config.limiting_current_density_method
                == LimitingCurrentDensityMethod.InitialValue
            ):
                return self.current_dens_lim_ioa[t] == (
                    self.config.limiting_current_density_data
                    * pyunits.amp
                    * pyunits.meter**-2
                    / sum(
                        self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * 0.5
                    * sum(
                        self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        + self.diluate.properties_out[t].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.electrode_side,
            self.flowsheet().time,
            self.ion_set,
            doc="Establish relationship between interfacial concentration polarization ratio and current density",
        )
        def eq_conc_polarization_ratio(self, mem, side, t, j):
            if (mem == "cem" and side == "cathode_left") or (
                mem == "aem" and side == "anode_right"
            ):
                return self.conc_mem_surf_mol_ioa[mem, side, t, j] / (
                    0.5
                    * (
                        self.concentrate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        + self.concentrate.properties_out[t].conc_mol_phase_comp[
                            "Liq", j
                        ]
                    )
                ) == (
                    1
                    + self.current[t]
                    * (self.cell_width * self.cell_length) ** -1
                    / self.current_dens_lim_ioa[t]
                )
            else:
                return self.conc_mem_surf_mol_ioa[mem, side, t, j] / (
                    0.5
                    * (
                        self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                        + self.diluate.properties_out[t].conc_mol_phase_comp["Liq", j]
                    )
                ) == (
                    1
                    - self.current[t]
                    * (self.cell_width * self.cell_length) ** -1
                    / self.current_dens_lim_ioa[t]
                )

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            doc="Calculate the total non-ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_potential_nonohm_dl_ioa(self, mem, t):
            if mem == "cem":
                return self.potential_nonohm_dl_ioa[mem, t] == (
                    Constants.gas_constant
                    * self.diluate.properties_in[t].temperature
                    / Constants.faraday_constant
                    * (
                        sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.anion_set
                        )
                    )
                    * log(
                        sum(
                            self.conc_mem_surf_mol_ioa[mem, "anode_right", t, j]
                            for j in self.ion_set
                        )
                        * 0.5
                        * (
                            sum(
                                self.concentrate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.concentrate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                        * sum(
                            self.conc_mem_surf_mol_ioa[mem, "cathode_left", t, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * (
                            0.5
                            * (
                                sum(
                                    self.diluate.properties_in[t].conc_mol_phase_comp[
                                        "Liq", j
                                    ]
                                    for j in self.ion_set
                                )
                                + sum(
                                    self.diluate.properties_out[t].conc_mol_phase_comp[
                                        "Liq", j
                                    ]
                                    for j in self.ion_set
                                )
                            )
                        )
                        ** -1
                    )
                )
            else:
                return self.potential_nonohm_dl_ioa[mem, t] == (
                    Constants.gas_constant
                    * self.diluate.properties_in[t].temperature
                    / Constants.faraday_constant
                    * (
                        sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.anion_set
                        )
                    )
                    * log(
                        sum(
                            self.conc_mem_surf_mol_ioa[mem, "anode_right", t, j]
                            for j in self.ion_set
                        )
                        * 0.5
                        * (
                            sum(
                                self.diluate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.diluate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                        * sum(
                            self.conc_mem_surf_mol_ioa[mem, "cathode_left", t, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * (
                            0.5
                            * (
                                sum(
                                    self.concentrate.properties_in[
                                        t
                                    ].conc_mol_phase_comp["Liq", j]
                                    for j in self.ion_set
                                )
                                + sum(
                                    self.concentrate.properties_out[
                                        t
                                    ].conc_mol_phase_comp["Liq", j]
                                    for j in self.ion_set
                                )
                            )
                        )
                        ** -1
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            doc="Calculate the total ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_potential_ohm_dl_ioa(self, mem, t):
            if mem == "cem":
                return self.potential_ohm_dl_ioa[mem, t] == (
                    Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * self.diluate.properties_in[t].equiv_conductivity_phase["Liq"]
                    ** -1
                    * log(
                        sum(
                            self.conc_mem_surf_mol_ioa[mem, "anode_right", t, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * 0.5**-1
                        * (
                            sum(
                                self.concentrate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.concentrate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                        ** -1
                        * sum(
                            self.conc_mem_surf_mol_ioa[mem, "cathode_left", t, j]
                            for j in self.ion_set
                        )
                        * 0.5
                        * (
                            sum(
                                self.diluate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.diluate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                    )
                )
            else:
                return self.potential_ohm_dl_ioa[mem, t] == (
                    -Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * self.diluate.properties_in[t].equiv_conductivity_phase["Liq"]
                    ** -1
                    * log(
                        sum(
                            self.conc_mem_surf_mol_ioa[mem, "anode_right", t, j]
                            for j in self.ion_set
                        )
                        * 0.5
                        * (
                            sum(
                                self.diluate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.diluate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                        * sum(
                            self.conc_mem_surf_mol_ioa[mem, "cathode_left", t, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * 0.5**-1
                        * (
                            sum(
                                self.concentrate.properties_in[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                            + sum(
                                self.concentrate.properties_out[t].conc_mol_phase_comp[
                                    "Liq", j
                                ]
                                for j in self.ion_set
                            )
                        )
                        ** -1
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.electrode_side,
            self.flowsheet().time,
            doc="Calculate the total non-ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_dl_thickness_ioa(self, mem, side, t):
            if mem == "cem" and side == "cathode_left":
                return self.dl_thickness_ioa[mem, side, t] == (
                    Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * 0.5
                    * (
                        sum(
                            self.concentrate.properties_in[t].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.cation_set
                        )
                        + sum(
                            self.concentrate.properties_out[t].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.cation_set
                        )
                    )
                    * self.current_dens_lim_ioa[t] ** -1
                )
            elif mem == "cem" and side == "anode_right":
                return self.dl_thickness_ioa[mem, side, t] == (
                    Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * 0.5
                    * (
                        sum(
                            self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties_out[t].conc_mol_phase_comp["Liq", j]
                            for j in self.cation_set
                        )
                    )
                    * self.current_dens_lim_ioa[t] ** -1
                )
            elif mem == "aem" and side == "cathode_left":
                return self.dl_thickness_ioa[mem, side, t] == (
                    -Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * 0.5
                    * (
                        sum(
                            self.diluate.properties_in[t].conc_mol_phase_comp["Liq", j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties_out[t].conc_mol_phase_comp["Liq", j]
                            for j in self.cation_set
                        )
                    )
                    * self.current_dens_lim_ioa[t] ** -1
                )
            else:
                return self.dl_thickness_ioa[mem, side, t] == (
                    -Constants.faraday_constant
                    * (
                        sum(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                            ** -1
                            for j in self.ion_set
                        )
                        ** -1
                        * len(self.ion_set)
                    )
                    * (
                        sum(
                            self.ion_trans_number_membrane[mem, j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        - sum(
                            self.diluate.properties_in[t].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * 0.5
                    * (
                        sum(
                            self.concentrate.properties_in[t].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.cation_set
                        )
                        + sum(
                            self.concentrate.properties_out[t].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.cation_set
                        )
                    )
                    * self.current_dens_lim_ioa[t] ** -1
                )

    # initialize method
    def initialize_build(
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

        # Set the outlet has the same intial condition of the inlet.
        for k in blk.keys():
            for j in blk[k].config.property_package.component_list:
                blk[k].diluate.properties_out[0].flow_mol_phase_comp["Liq", j] = value(
                    blk[k].diluate.properties_in[0].flow_mol_phase_comp["Liq", j]
                )
                blk[k].concentrate.properties_out[0].flow_mol_phase_comp[
                    "Liq", j
                ] = value(
                    blk[k].concentrate.properties_in[0].flow_mol_phase_comp["Liq", j]
                )
        if hasattr(blk[k], "conc_mem_surf_mol_ioa"):
            for mem in blk[k].membrane_set:
                for side in blk[k].electrode_side:
                    for j in blk[k].ion_set:
                        blk[k].conc_mem_surf_mol_ioa[mem, side, 0, j].set_value(
                            blk[k]
                            .concentrate.properties_in[0]
                            .conc_mol_phase_comp["Liq", j]
                        )
        # Initialize diluate block
        flags_diluate = blk.diluate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize concentrate_side block
        flags_concentrate = blk.concentrate.initialize(
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
        blk.diluate.release_state(flags_diluate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        blk.concentrate.release_state(flags_concentrate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Scaling factors that user may setup
        # The users are highly encouraged to provide scaling factors for assessable vars below.
        # Not providing these vars will give a warning.
        if (
            iscale.get_scaling_factor(self.solute_diffusivity_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.solute_diffusivity_membrane, 1e10)
        if iscale.get_scaling_factor(self.membrane_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e4)
        if (
            iscale.get_scaling_factor(self.water_permeability_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.water_permeability_membrane, 1e14)
        if iscale.get_scaling_factor(self.cell_pair_num, warning=True) is None:
            iscale.set_scaling_factor(self.cell_pair_num, 0.1)
        if iscale.get_scaling_factor(self.cell_length, warning=True) is None:
            iscale.set_scaling_factor(self.cell_length, 1e1)
        if iscale.get_scaling_factor(self.cell_width, warning=True) is None:
            iscale.set_scaling_factor(self.cell_width, 1e1)
        if iscale.get_scaling_factor(self.spacer_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.spacer_thickness, 1e4)
        if (
            iscale.get_scaling_factor(self.membrane_areal_resistance, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.membrane_areal_resistance, 1e4)
        if iscale.get_scaling_factor(self.electrodes_resistance, warning=True) is None:
            iscale.set_scaling_factor(self.electrodes_resistance, 1e4)
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

        for ind in self.nonelec_flux_in:
            if ind[2] == "H2O":
                sf = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * iscale.get_scaling_factor(
                        self.concentrate.properties_in[ind[0]].pressure_osm_phase[
                            ind[1]
                        ]
                    )
                )
            else:
                sf = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * iscale.get_scaling_factor(
                        self.concentrate.properties_in[ind[0]].conc_mol_phase_comp[
                            ind[1], ind[2]
                        ]
                    )
                )
            iscale.set_scaling_factor(self.nonelec_flux_in[ind], sf)
        for ind in self.nonelec_flux_out:
            if ind[2] == "H2O":
                sf = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * iscale.get_scaling_factor(
                        self.concentrate.properties_out[ind[0]].pressure_osm_phase[
                            ind[1]
                        ]
                    )
                )
            else:
                sf = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * iscale.get_scaling_factor(
                        self.concentrate.properties_out[ind[0]].conc_mol_phase_comp[
                            ind[1], ind[2]
                        ]
                    )
                )

            iscale.set_scaling_factor(self.nonelec_flux_out[ind], sf)
        iscale.set_scaling_factor(
            self.power_electrical,
            iscale.get_scaling_factor(self.current)
            * iscale.get_scaling_factor(self.voltage),
        )
        for ind, c in self.specific_power_electrical.items():
            iscale.set_scaling_factor(
                self.specific_power_electrical[ind],
                3.6e6
                * iscale.get_scaling_factor(self.current[ind])
                * iscale.get_scaling_factor(self.voltage[ind])
                * iscale.get_scaling_factor(
                    self.diluate.properties_out[ind].flow_vol_phase["Liq"]
                )
                ** -1,
            )
        if hasattr(self, "conc_mem_surf_mol_ioa"):
            for ind in self.conc_mem_surf_mol_ioa:
                if iscale.get_scaling_factor(self.conc_mem_surf_mol_ioa[ind]) is None:
                    if (ind[0] == "cem" and ind[1] == "cathode_left") or (
                        ind[0] == "aem" and ind[1] == "anode_right"
                    ):
                        iscale.set_scaling_factor(
                            self.conc_mem_surf_mol_ioa[ind],
                            iscale.get_scaling_factor(
                                self.concentrate.properties_in[
                                    ind[2]
                                ].conc_mol_phase_comp["Liq", ind[3]]
                            ),
                        )
                    else:
                        iscale.set_scaling_factor(
                            self.conc_mem_surf_mol_ioa[ind],
                            iscale.get_scaling_factor(
                                self.diluate.properties_in[ind[2]].conc_mol_phase_comp[
                                    "Liq", ind[3]
                                ]
                            ),
                        )

        if hasattr(self, "current_dens_lim_ioa"):
            if iscale.get_scaling_factor(self.current_dens_lim_ioa) is None:
                if (
                    self.config.limiting_current_density_method
                    == LimitingCurrentDensityMethod.InitialValue
                ):
                    sf = self.config.limiting_current_density_data**-1
                    iscale.set_scaling_factor(self.current_dens_lim_ioa, sf)

        if hasattr(self, "potential_nonohm_membrane_ioa"):
            if iscale.get_scaling_factor(self.potential_nonohm_membrane_ioa) is None:
                sf = (
                    value(Constants.faraday_constant)
                    * value(Constants.gas_constant) ** -1
                    * 298.15**-1
                )
                iscale.set_scaling_factor(self.potential_nonohm_membrane_ioa, sf)
        if hasattr(self, "potential_nonohm_dl_ioa"):
            if iscale.get_scaling_factor(self.potential_nonohm_dl_ioa) is None:
                sf = (
                    value(Constants.faraday_constant)
                    * value(Constants.gas_constant) ** -1
                    * 298.15**-1
                )
                iscale.set_scaling_factor(self.potential_nonohm_dl_ioa, sf)
        if hasattr(self, "potential_ohm_dl_ioa"):
            if iscale.get_scaling_factor(self.potential_ohm_dl_ioa) is None:
                sf = (
                    96485**-1
                    * sum(
                        iscale.get_scaling_factor(
                            self.config.property_package.diffus_phase_comp["Liq", j]
                        )
                        ** -2
                        for j in self.ion_set
                    )
                    ** -0.5
                    * float(len(self.ion_set)) ** -1
                )
                iscale.set_scaling_factor(self.potential_ohm_dl_ioa, sf)
        if hasattr(self, "dl_thickness_ioa"):
            if iscale.get_scaling_factor(self.dl_thickness_ioa) is None:
                for ind in self.dl_thickness_ioa:
                    sf = (
                        96485**-1
                        * sum(
                            iscale.get_scaling_factor(
                                self.config.property_package.diffus_phase_comp["Liq", j]
                            )
                            ** -2
                            for j in self.ion_set
                        )
                        ** -0.5
                        * len(self.ion_set) ** -1
                        * sum(
                            iscale.get_scaling_factor(
                                self.conc_mem_surf_mol_ioa[ind, j]
                            )
                            ** 2
                            for j in self.cation_set
                        )
                        ** 0.5
                        * iscale.get_scaling_factor(self.current) ** -1
                        * iscale.get_scaling_factor(self.cell_width)
                        * iscale.get_scaling_factor(self.cell_length)
                    )
                    iscale.set_scaling_factor(self.dl_thickness_ioa[ind], sf)

        # Constraint scaling
        for ind, c in self.eq_current_voltage_relation.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.membrane_areal_resistance)
            )
        for ind, c in self.eq_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.power_electrical)
            )
        for ind, c in self.eq_specific_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.specific_power_electrical[ind])
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
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.nonelec_flux_in[ind])
            )

        for ind, c in self.eq_nonelec_flux_out.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.nonelec_flux_out[ind])
            )

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
        if hasattr(self, "eq_potential_nonohm_membrane_ioa"):
            for ind, c in self.eq_potential_nonohm_membrane_ioa.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.potential_nonohm_membrane_ioa[ind]),
                )
        if hasattr(self, "eq_current_dens_lim_ioa"):
            for ind, c in self.eq_current_dens_lim_ioa.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.current_dens_lim_ioa[ind])
                )
        if hasattr(self, "eq_potential_nonohm_dl_ioa"):
            for ind, c in self.eq_potential_nonohm_dl_ioa.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.potential_nonohm_dl_ioa[ind])
                )
        if hasattr(self, "eq_potential_ohm_dl_ioa"):
            for ind, c in self.eq_potential_ohm_dl_ioa.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.potential_ohm_dl_ioa[ind])
                )
        if hasattr(self, "eq_dl_thickness_ioa"):
            for ind, c in self.eq_dl_thickness_ioa.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.dl_thickness_ioa[ind])
                )

        for ind, c in self.eq_recovery_mass_H2O.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.diluate.properties_out[ind].flow_mass_phase_comp["Liq", "H2O"]
                ),
            )

        for ind, c in self.eq_power_electrical.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(self.power_electrical[ind]),
            )

        for ind, c in self.eq_specific_power_electrical.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(self.specific_power_electrical[ind])
                * iscale.get_scaling_factor(
                    self.diluate.properties_out[ind].flow_vol_phase["Liq"]
                ),
            )
        for ind, c in self.eq_current_efficiency.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.current[ind])
            )

        for ind, c in self.eq_isothermal_diluate.items():
            iscale.constraint_scaling_transform(
                c, self.diluate.properties_in[ind].temperature
            )
        for ind, c in self.eq_isothermal_concentrate.items():
            iscale.constraint_scaling_transform(
                c, self.concentrate.properties_in[ind].temperature
            )

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

    def _get_performance_contents(self, time_point=0):
        return {
            "vars": {
                "Total electrical power consumption(Watt)": self.power_electrical[
                    time_point
                ],
                "Specific electrical power consumption (kW*h/m**3)": self.specific_power_electrical[
                    time_point
                ],
                "Current efficiency for deionzation": self.current_efficiency[
                    time_point
                ],
                "Water recovery by mass": self.recovery_mass_H2O[time_point],
            },
            "exprs": {},
            "params": {},
        }

    def get_power_electrical(self, time_point=0):
        return self.power_electrical[time_point]
