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

from pyomo.environ import (
    Set,
    Var,
    Suffix,
    Constraint,
    NonNegativeReals,
    NonNegativeIntegers,
    value,
    units as pyunits,
    log,
)
from pyomo.dae import (
    DerivativeVar,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import Watertap cores
from watertap.core.util.initialization import check_solve, check_dof

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.constants import Constants
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from enum import Enum

__author__ = "Xiangyu Bi, Austin Ladshaw"

_log = idaeslog.getLogger(__name__)


class LimitingCurrentDensityMethod(Enum):
    InitialValue = 0
    # Empirical = 1
    # Theoretical = 2 TODO: 1 and 2


class ElectricalOperationMode(Enum):
    Constant_Current = 0
    Constant_Voltage = 1


# Name of the unit model
@declare_process_block_class("Electrodialysis1D")
class Electrodialysis1DData(UnitModelBlockData):
    """
    1D Electrodialysis Model
    """

    # CONFIG are options for the unit model
    CONFIG = ConfigBlock()

    # These config args are common to any control volume
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
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "operation_mode",
        ConfigValue(
            default=ElectricalOperationMode.Constant_Voltage,
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
           "``LimitingCurrentDensityMethod.Empirical``", "Limiting current density is calculated from the empirical equation: TODO "
           "``LimitingCurrentDensityMethod.Theoretical``", "Limiting current density is calculated from a theoretical equation: TODO "
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

    # # Consider adding the EnergyBalanceType config using the following code
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

    # These config args are specifically for 1D control volumes
    '''
    CONFIG.declare(
        "area_definition",
        ConfigValue(
            default=DistributedVars.uniform,
            domain=In(DistributedVars),
            description="Argument for defining form of area variable",
            doc="""Argument defining whether area variable should be spatially
    variant or not. **default** - DistributedVars.uniform.
    **Valid values:** {
    DistributedVars.uniform - area does not vary across spatial domain,
    DistributedVars.variant - area can vary over the domain and is indexed
    by time and space.}""",
        ),
    )
    '''

    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.finite_difference",
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
    documentation for supported transformations.""",
        ),
    )

    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default="BACKWARD",
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See
    Pyomo documentation for supported schemes.""",
        ),
    )

    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=10,
            domain=int,
            description="Number of finite elements in length domain",
            doc="""Number of finite elements to use when discretizing length
            domain (default=10)""",
        ),
    )

    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=2,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=2)""",
        ),
    )

    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()
        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # Get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units
        # Create essential sets
        add_object_reference(
            self, "component_set", self.config.property_package.component_list
        )
        add_object_reference(
            self, "solute_set", self.config.property_package.solute_set
        )
        add_object_reference(self, "ion_set", self.config.property_package.ion_set)

        add_object_reference(
            self, "cation_set", self.config.property_package.cation_set
        )
        add_object_reference(self, "anion_set", self.config.property_package.anion_set)
        self.membrane_set = Set(
            initialize=["cem", "aem"]
        )  #   cem = Cation-Exchange Membrane aem = Anion-Exchange Membrane
        self.electrode_side_set = Set(initialize=["cathode_left", "anode_right"])
        # Length var for building 1D control volume
        self.cell_length = Var(
            initialize=0.5,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="Length of the electrodialysis cell/stack, in parallel to the convective flow",
        )

        # Control Volume for the Diluate channel:
        self.diluate = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )
        self.diluate.add_geometry(length_var=self.cell_length)
        self.diluate.add_state_blocks(has_phase_equilibrium=False)
        self.diluate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        # # Note: the energy balance is disabled currently
        if hasattr(self.config, "energy_balance_type"):
            self.diluate.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )
        self.diluate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Below is declared the electrical power var and its derivative var,
        # which is a performance metric of the entire electrodialysis stack.
        # This var takes the "diluate" as the parent to utilize the discretization (as in Pyomo DAE) of this block for solving.
        self.diluate.power_electrical_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0,
            bounds=(0, 12100),
            domain=NonNegativeReals,
            units=pyunits.watt,
            doc="Electrical power consumption of a stack",
        )
        self.diluate.Dpower_electrical_Dx = DerivativeVar(
            self.diluate.power_electrical_x,
            wrt=self.diluate.length_domain,
            units=pyunits.watt,
        )

        # Apply the discretization transformation (Pyomo DAE) to the diluate block
        self.diluate.apply_transformation()

        # Control volume for the Concentrate channel
        self.concentrate = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )
        self.concentrate.add_geometry(length_var=self.cell_length)
        self.concentrate.add_state_blocks(has_phase_equilibrium=False)
        self.concentrate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        # # Note: the energy balance is disabled currently
        if hasattr(self.config, "energy_balance_type"):
            self.concentrate.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )
        self.concentrate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )
        self.concentrate.apply_transformation()

        # Add inlet and outlet as ports for Diluate and Concentrate
        self.add_inlet_port(name="inlet_diluate", block=self.diluate)
        self.add_outlet_port(name="outlet_diluate", block=self.diluate)
        self.add_inlet_port(name="inlet_concentrate", block=self.concentrate)
        self.add_outlet_port(name="outlet_concentrate", block=self.concentrate)

        # Create unit model parameters and vars
        self.cell_pair_num = Var(
            initialize=1,
            domain=NonNegativeReals,
            bounds=(1, 10000),
            units=pyunits.dimensionless,
            doc="cell pair number in a stack",
        )
        # Note: domain is set as Reals for optimizaiton convenience,
        # although this var is technically integer.

        # electrodialysis cell dimensional properties
        self.cell_width = Var(
            initialize=0.1,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="The width of the electrodialysis cell, parallel to the current/ion flux direction",
        )
        self.spacer_thickness = Var(
            initialize=0.0001,
            units=pyunits.meter,
            doc="The distance between consecutive aem and cem",
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
            self.config.property_package.ion_set
            | self.config.property_package.solute_set,
            initialize=1e-10,
            bounds=(1e-16, 1e-6),
            units=pyunits.meter**2 * pyunits.second**-1,
            doc="Solute (ionic and neutral) diffusivity in the membrane phase",
        )
        self.ion_trans_number_membrane = Var(
            self.membrane_set,
            self.config.property_package.ion_set,
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
        self.total_areal_resistance_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1e-2,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.ohm * pyunits.meter**2,
            doc="Total areal resistance of a stack ",
        )
        if self.config.operation_mode == ElectricalOperationMode.Constant_Current:
            self.current_applied = Var(
                self.flowsheet().time,
                initialize=1,
                bounds=(0, 1000),
                units=pyunits.amp,
                doc="Current across a cell-pair or stack, declared under the 'Constant Current' mode only",
            )

        self.current_density_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp * pyunits.meter**-2,
            doc="Current density accross the membrane as a function of the normalized length",
        )

        if self.config.operation_mode == ElectricalOperationMode.Constant_Voltage:
            self.voltage_applied = Var(
                self.flowsheet().time,
                initialize=100,
                bounds=(0, 1000),
                units=pyunits.volt,
                doc="Voltage across a stack, declared under the 'Constant Voltage' mode only",
            )

        self.voltage_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=100,
            bounds=(0, 1000),
            units=pyunits.volt,
            doc="Voltage across a stack",
        )
        self.current_utilization = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The current utilization including water electro-osmosis and ion diffusion",
        )
        # Performance metrics
        # Note that the power_electrical_x was declared under the diluate block
        self.specific_power_electrical = Var(
            self.flowsheet().time,
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.kW * pyunits.hour * pyunits.meter**-3,
            doc="Diluate-volume-flow-rate-specific electrical power consumption",
        )
        self.current_efficiency_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The overall current efficiency for deionizaiton",
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
        self._make_performance()

    def _make_performance(self):
        if (
            self.config.has_nonohmic_potential_membrane
            or self.config.has_Nernst_diffusion_layer
        ):
            self.conc_mem_surf_mol_x = Var(
                self.membrane_set,
                self.electrode_side_set,
                self.flowsheet().time,
                self.diluate.length_domain,
                self.config.property_package.ion_set,
                initialize=500,
                bounds=(0, 1e5),
                units=pyunits.mol * pyunits.meter**-3,
                doc="Membane surface concentration of components",
            )
        if self.config.has_nonohmic_potential_membrane:
            self._make_performance_nonohm_mem()
        if self.config.has_Nernst_diffusion_layer:
            self.current_dens_lim_x = Var(
                self.flowsheet().time,
                self.diluate.length_domain,
                initialize=500,
                bounds=(0, 10000),
                units=pyunits.amp * pyunits.meter**-2,
                doc="Limiting Current Density accross the membrane as a function of the normalized length",
            )
            self._make_performance_dl_polarization()
        # To require H2O must be in the component
        if "H2O" not in self.component_set:
            raise ConfigurationError(
                "Property Package MUST constain 'H2O' as a component"
            )

        # -------- Add constraints ---------

        # Adds isothermal constraint if no energy balance present
        if not hasattr(self.config, "energy_balance_type"):

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Isothermal condition for Diluate",
            )
            def eq_isothermal_diluate(self, t, x):
                if x == self.diluate.length_domain.first():
                    return Constraint.Skip
                return (
                    self.diluate.properties[
                        t, self.diluate.length_domain.first()
                    ].temperature
                    == self.diluate.properties[t, x].temperature
                )

        if not hasattr(self.config, "energy_balance_type"):

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Isothermal condition for Concentrate",
            )
            def eq_isothermal_concentrate(self, t, x):
                if x == self.diluate.length_domain.first():
                    return Constraint.Skip
                return (
                    self.concentrate.properties[
                        t, self.diluate.length_domain.first()
                    ].temperature
                    == self.concentrate.properties[t, x].temperature
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total areal resistance of a stack",
        )
        def eq_get_total_areal_resistance_x(self, t, x):
            if self.config.has_Nernst_diffusion_layer:
                return self.total_areal_resistance_x[t, x] == (
                    (
                        self.membrane_areal_resistance["aem"]
                        + self.membrane_areal_resistance["cem"]
                        + (
                            self.spacer_thickness
                            - self.dl_thickness_x["cem", "cathode_left", t, x]
                            - self.dl_thickness_x["aem", "anode_right", t, x]
                        )
                        * self.concentrate.properties[t, x].elec_cond_phase["Liq"] ** -1
                        + (
                            self.spacer_thickness
                            - self.dl_thickness_x["cem", "anode_right", t, x]
                            - self.dl_thickness_x["aem", "cathode_left", t, x]
                        )
                        * self.diluate.properties[t, x].elec_cond_phase["Liq"] ** -1
                    )
                    * self.cell_pair_num
                    + self.electrodes_resistance
                )
            else:
                return self.total_areal_resistance_x[t, x] == (
                    (
                        self.membrane_areal_resistance["aem"]
                        + self.membrane_areal_resistance["cem"]
                        + self.spacer_thickness
                        * (
                            self.concentrate.properties[t, x].elec_cond_phase["Liq"]
                            ** -1
                            + self.diluate.properties[t, x].elec_cond_phase["Liq"] ** -1
                        )
                    )
                    * self.cell_pair_num
                    + self.electrodes_resistance
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="calcualte current density from the electrical input",
        )
        def eq_get_current_density(self, t, x):
            if self.config.operation_mode == ElectricalOperationMode.Constant_Current:
                return (
                    self.current_density_x[t, x] * self.cell_width * self.diluate.length
                    == self.current_applied[t]
                )
            else:
                if self.config.has_nonohmic_potential_membrane:
                    if self.config.has_Nernst_diffusion_layer:
                        return (
                            self.current_density_x[t, x]
                            * self.total_areal_resistance_x[t, x]
                            + (
                                self.potential_ohm_dl_x["cem", t, x]
                                + self.potential_ohm_dl_x["aem", t, x]
                                + self.potential_nonohm_dl_x["cem", t, x]
                                + self.potential_nonohm_dl_x["aem", t, x]
                                + self.potential_nonohm_membrane_x["cem", t, x]
                                + self.potential_nonohm_membrane_x["aem", t, x]
                            )
                            * self.cell_pair_num
                            == self.voltage_applied[t]
                        )
                    else:
                        return (
                            self.current_density_x[t, x]
                            * self.total_areal_resistance_x[t, x]
                            + (
                                self.potential_nonohm_membrane_x["cem", t, x]
                                + self.potential_nonohm_membrane_x["aem", t, x]
                            )
                            * self.cell_pair_num
                            == self.voltage_applied[t]
                        )
                else:
                    if self.config.has_Nernst_diffusion_layer:
                        return (
                            self.current_density_x[t, x]
                            * self.total_areal_resistance_x[t, x]
                            + (
                                self.potential_ohm_dl_x["cem", t, x]
                                + self.potential_ohm_dl_x["aem", t, x]
                                + self.potential_nonohm_dl_x["cem", t, x]
                                + self.potential_nonohm_dl_x["aem", t, x]
                            )
                            * self.cell_pair_num
                            == self.voltage_applied[t]
                        )
                    else:
                        return (
                            self.current_density_x[t, x]
                            * self.total_areal_resistance_x[t, x]
                            == self.voltage_applied[t]
                        )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="calcualte length_indexed voltage",
        )
        def eq_get_voltage_x(self, t, x):
            if self.config.has_nonohmic_potential_membrane:
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.current_density_x[t, x]
                        * self.total_areal_resistance_x[t, x]
                        + (
                            self.potential_ohm_dl_x["cem", t, x]
                            + self.potential_ohm_dl_x["aem", t, x]
                            + self.potential_nonohm_dl_x["cem", t, x]
                            + self.potential_nonohm_dl_x["aem", t, x]
                            + self.potential_nonohm_membrane_x["cem", t, x]
                            + self.potential_nonohm_membrane_x["aem", t, x]
                        )
                        * self.cell_pair_num
                        == self.voltage_x[t, x]
                    )
                else:
                    return (
                        self.current_density_x[t, x]
                        * self.total_areal_resistance_x[t, x]
                        + (
                            self.potential_nonohm_membrane_x["cem", t, x]
                            + self.potential_nonohm_membrane_x["aem", t, x]
                        )
                        * self.cell_pair_num
                        == self.voltage_x[t, x]
                    )
            else:
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.current_density_x[t, x]
                        * self.total_areal_resistance_x[t, x]
                        + (
                            self.potential_ohm_dl_x["cem", t, x]
                            + self.potential_ohm_dl_x["aem", t, x]
                            + self.potential_nonohm_dl_x["cem", t, x]
                            + self.potential_nonohm_dl_x["aem", t, x]
                        )
                        * self.cell_pair_num
                        == self.voltage_x[t, x]
                    )
                else:
                    return (
                        self.current_density_x[t, x]
                        * self.total_areal_resistance_x[t, x]
                        == self.voltage_x[t, x]
                    )

        # Mass Transfer for the Diluate
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term diluate channel",
        )
        def eq_mass_transfer_term_diluate(self, t, x, p, j):

            if j == "H2O":
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.diluate.mass_transfer_term[t, x, p, j]
                        == (
                            -self.cell_width
                            * (
                                self.water_trans_number_membrane["cem"]
                                + self.water_trans_number_membrane["aem"]
                            )
                            * (
                                self.current_density_x[t, x]
                                / Constants.faraday_constant
                            )
                            - self.cell_width
                            * self.diluate.properties[
                                0, self.diluate.length_domain.first()
                            ].dens_mass_solvent
                            / self.config.property_package.mw_comp[j]
                            * (
                                self.water_permeability_membrane["cem"]
                                + self.water_permeability_membrane["aem"]
                            )
                            * (
                                self.concentrate.properties[t, x].pressure_osm_phase[p]
                                * (
                                    1
                                    + self.current_density_x[t, x]
                                    / self.current_dens_lim_x[t, x]
                                )
                                - self.diluate.properties[t, x].pressure_osm_phase[p]
                                * (
                                    1
                                    - self.current_density_x[t, x]
                                    / self.current_dens_lim_x[t, x]
                                )
                            )
                        )
                        * self.cell_pair_num
                    )
                else:
                    return (
                        self.diluate.mass_transfer_term[t, x, p, j]
                        == (
                            -self.cell_width
                            * (
                                self.water_trans_number_membrane["cem"]
                                + self.water_trans_number_membrane["aem"]
                            )
                            * (
                                self.current_density_x[t, x]
                                / Constants.faraday_constant
                            )
                            - self.cell_width
                            * self.diluate.properties[
                                0, self.diluate.length_domain.first()
                            ].dens_mass_solvent
                            / self.config.property_package.mw_comp[j]
                            * (
                                self.water_permeability_membrane["cem"]
                                + self.water_permeability_membrane["aem"]
                            )
                            * (
                                self.concentrate.properties[t, x].pressure_osm_phase[p]
                                - self.diluate.properties[t, x].pressure_osm_phase[p]
                            )
                        )
                        * self.cell_pair_num
                    )
            elif j in self.config.property_package.ion_set:
                if self.config.has_Nernst_diffusion_layer:
                    return (
                        self.diluate.mass_transfer_term[t, x, p, j]
                        == (
                            -self.cell_width
                            * (
                                self.ion_trans_number_membrane["cem", j]
                                - self.ion_trans_number_membrane["aem", j]
                            )
                            * (self.current_utilization * self.current_density_x[t, x])
                            / (
                                self.config.property_package.charge_comp[j]
                                * Constants.faraday_constant
                            )
                            + self.cell_width
                            * (
                                self.solute_diffusivity_membrane["cem", j]
                                / self.membrane_thickness["cem"]
                                * (
                                    self.conc_mem_surf_mol_x[
                                        "cem", "cathode_left", t, x, j
                                    ]
                                    - self.conc_mem_surf_mol_x[
                                        "cem", "anode_right", t, x, j
                                    ]
                                )
                                + self.solute_diffusivity_membrane["aem", j]
                                / self.membrane_thickness["aem"]
                                * (
                                    self.conc_mem_surf_mol_x[
                                        "aem", "anode_right", t, x, j
                                    ]
                                    - self.conc_mem_surf_mol_x[
                                        "aem", "cathode_left", t, x, j
                                    ]
                                )
                            )
                        )
                        * self.cell_pair_num
                    )
                else:
                    return (
                        self.diluate.mass_transfer_term[t, x, p, j]
                        == (
                            -self.cell_width
                            * (
                                self.ion_trans_number_membrane["cem", j]
                                - self.ion_trans_number_membrane["aem", j]
                            )
                            * (self.current_utilization * self.current_density_x[t, x])
                            / (
                                self.config.property_package.charge_comp[j]
                                * Constants.faraday_constant
                            )
                            + self.cell_width
                            * (
                                self.solute_diffusivity_membrane["cem", j]
                                / self.membrane_thickness["cem"]
                                + self.solute_diffusivity_membrane["aem", j]
                                / self.membrane_thickness["aem"]
                            )
                            * (
                                self.concentrate.properties[t, x].conc_mol_phase_comp[
                                    p, j
                                ]
                                - self.diluate.properties[t, x].conc_mol_phase_comp[
                                    p, j
                                ]
                            )
                        )
                        * self.cell_pair_num
                    )
            else:
                return (
                    self.diluate.mass_transfer_term[t, x, p, j]
                    == (
                        self.cell_width
                        * (
                            self.solute_diffusivity_membrane["cem", j]
                            / self.membrane_thickness["cem"]
                            + self.solute_diffusivity_membrane["aem", j]
                            / self.membrane_thickness["aem"]
                        )
                        * (
                            self.concentrate.properties[t, x].conc_mol_phase_comp[p, j]
                            - self.diluate.properties[t, x].conc_mol_phase_comp[p, j]
                        )
                    )
                    * self.cell_pair_num
                )

        # Mass Transfer for the Concentrate
        @self.Constraint(
            self.flowsheet().config.time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term concentrate channel",
        )
        def eq_mass_transfer_term_concentrate(self, t, x, p, j):

            return (
                self.concentrate.mass_transfer_term[t, x, p, j]
                == -self.diluate.mass_transfer_term[t, x, p, j]
            )

        # Performance Metrics calculation
        @self.Constraint(
            self.flowsheet().config.time,
            self.diluate.length_domain,
            doc="Electrical power consumption of a stack",
        )
        def eq_power_electrical(self, t, x):
            if x == self.diluate.length_domain.first():
                self.diluate.power_electrical_x[t, x].fix(0)
                return Constraint.Skip
            else:
                return (
                    self.diluate.Dpower_electrical_Dx[t, x]
                    == self.voltage_x[t, x]
                    * self.current_density_x[t, x]
                    * self.cell_width
                    * self.diluate.length
                )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Diluate_volume_flow_rate_specific electrical power consumption of a stack",
        )
        def eq_specific_power_electrical(self, t):
            return (
                pyunits.convert(
                    self.specific_power_electrical[t],
                    pyunits.watt * pyunits.second * pyunits.meter**-3,
                )
                * self.diluate.properties[
                    t, self.diluate.length_domain.last()
                ].flow_vol_phase["Liq"]
                == self.diluate.power_electrical_x[t, self.diluate.length_domain.last()]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            doc="Overall current efficiency evaluation",
        )
        def eq_current_efficiency_x(self, t, x, p):

            return (
                self.current_efficiency_x[t, x]
                * self.current_density_x[t, x]
                * self.cell_width
                * self.cell_pair_num
                == -sum(
                    self.diluate.mass_transfer_term[t, x, p, j]
                    * self.config.property_package.charge_comp[j]
                    for j in self.cation_set
                )
                * Constants.faraday_constant
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Water recovery by mass",
        )
        def eq_recovery_mass_H2O(self, t):
            return (
                self.recovery_mass_H2O[t]
                * (
                    self.diluate.properties[
                        t, self.diluate.length_domain.first()
                    ].flow_mass_phase_comp["Liq", "H2O"]
                    + self.concentrate.properties[
                        t, self.diluate.length_domain.first()
                    ].flow_mass_phase_comp["Liq", "H2O"]
                )
                == self.diluate.properties[
                    t, self.diluate.length_domain.last()
                ].flow_mass_phase_comp["Liq", "H2O"]
            )

    def _make_performance_nonohm_mem(self):

        self.potential_nonohm_membrane_x = Var(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,  # to inspect
            bounds=(-50, 50),
            units=pyunits.volt,
            doc="Nonohmic potential across a membane",
        )

        @self.Constraint(
            self.membrane_set,
            self.electrode_side_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            self.ion_set,
            doc="calcualte current density from the electrical input",
        )
        def eq_set_surface_conc(self, mem, side, t, x, j):
            if not self.config.has_Nernst_diffusion_layer:
                if (mem == "cem" and side == "cathode_left") or (
                    mem == "aem" and side == "anode_right"
                ):
                    return (
                        self.conc_mem_surf_mol_x[mem, side, t, x, j]
                        == self.concentrate.properties[t, x].conc_mol_phase_comp[
                            "Liq", j
                        ]
                    )
                else:
                    return (
                        self.conc_mem_surf_mol_x[mem, side, t, x, j]
                        == self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                    )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total non-ohmic potential across an iem; this takes account of diffusion and Donnan Potentials",
        )
        def eq_potential_nonohm_membrane_x(self, mem, t, x):
            if (
                not self.config.has_Nernst_diffusion_layer
            ) and x == self.diluate.length_domain.first():
                self.potential_nonohm_membrane_x[mem, t, x].fix(0)
                return Constraint.Skip

            return self.potential_nonohm_membrane_x[mem, t, x] == (
                Constants.gas_constant
                * self.diluate.properties[t, x].temperature
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
                        self.conc_mem_surf_mol_x[mem, "cathode_left", t, x, j]
                        for j in self.ion_set
                    )
                    / sum(
                        self.conc_mem_surf_mol_x[mem, "anode_right", t, x, j]
                        for j in self.ion_set
                    )
                )
            )  # Note: the log argument is taken as *total ion concentration* while in the original equation as *electrolyte concentration*
            # This is mathematically valid as the ratio eliminate the difference and electro-neutrality stands.

    def _make_performance_dl_polarization(self):

        self.potential_nonohm_dl_x = Var(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,  # to inspect
            bounds=(-50, 50),
            units=pyunits.volt,
            doc="Nonohmic potential in two diffusion layers on the two sides of a membrane",
        )

        self.potential_ohm_dl_x = Var(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,  # to inspect
            bounds=(0, 50),
            units=pyunits.volt,
            doc="Ohmic potential in two diffusion layers on the two sides of a membrane",
        )

        self.dl_thickness_x = Var(
            self.membrane_set,
            self.electrode_side_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.0005,
            bounds=(0, 1e-2),
            units=pyunits.m,
            doc="Thickness of the diffusion layer",
        )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate length-indexed limitting current density",
        )
        def eq_current_dens_lim_x(self, t, x):
            if (
                self.config.limiting_current_density_method
                == LimitingCurrentDensityMethod.InitialValue
            ):
                return self.current_dens_lim_x[t, x] == (
                    self.config.limiting_current_density_data
                    * pyunits.amp
                    * pyunits.meter**-2
                    / sum(
                        self.diluate.properties[t, 0].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * sum(
                        self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.electrode_side_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.ion_set,
            doc="Establish relationship between interfacial concentration polarization ratio and current density",
        )
        def eq_conc_polarization_ratio(self, mem, side, t, x, j):
            if (mem == "cem" and side == "cathode_left") or (
                mem == "aem" and side == "anode_right"
            ):
                return self.conc_mem_surf_mol_x[
                    mem, side, t, x, j
                ] / self.concentrate.properties[t, x].conc_mol_phase_comp["Liq", j] == (
                    1 + self.current_density_x[t, x] / self.current_dens_lim_x[t, x]
                )
            else:
                return self.conc_mem_surf_mol_x[
                    mem, side, t, x, j
                ] / self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j] == (
                    1 - self.current_density_x[t, x] / self.current_dens_lim_x[t, x]
                )

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total non-ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_potential_nonohm_dl(self, mem, t, x):
            if mem == "cem":
                return self.potential_nonohm_dl_x[mem, t, x] == (
                    Constants.gas_constant
                    * self.diluate.properties[t, x].temperature
                    / Constants.faraday_constant
                    * (
                        sum(
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.anion_set
                        )
                    )
                    * log(
                        sum(
                            self.conc_mem_surf_mol_x[mem, "anode_right", t, x, j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.concentrate.properties[t, x].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.ion_set
                        )
                        * sum(
                            self.conc_mem_surf_mol_x[mem, "cathode_left", t, x, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * sum(
                            self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                            for j in self.ion_set
                        )
                        ** -1
                    )
                )
            else:
                return self.potential_nonohm_dl_x[mem, t, x] == (
                    Constants.gas_constant
                    * self.diluate.properties[t, x].temperature
                    / Constants.faraday_constant
                    * (
                        sum(
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                        + sum(
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.anion_set
                        )
                    )
                    * log(
                        sum(
                            self.conc_mem_surf_mol_x[mem, "anode_right", t, x, j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.conc_mem_surf_mol_x[mem, "cathode_left", t, x, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * sum(
                            self.concentrate.properties[t, x].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.ion_set
                        )
                        ** -1
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_potential_ohm_dl_x(self, mem, t, x):
            if mem == "cem":
                return self.potential_ohm_dl_x[mem, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * self.diluate.properties[t, x].equiv_conductivity_phase["Liq"]
                    ** -1
                    * log(
                        sum(
                            self.conc_mem_surf_mol_x[mem, "anode_right", t, x, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * sum(
                            self.concentrate.properties[t, x].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.ion_set
                        )
                        ** -1
                        * sum(
                            self.conc_mem_surf_mol_x[mem, "cathode_left", t, x, j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                            for j in self.ion_set
                        )
                    )
                )
            else:
                return self.potential_ohm_dl_x[mem, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * self.diluate.properties[t, x].equiv_conductivity_phase["Liq"]
                    ** -1
                    * log(
                        sum(
                            self.conc_mem_surf_mol_x[mem, "anode_right", t, x, j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                            for j in self.ion_set
                        )
                        * sum(
                            self.conc_mem_surf_mol_x[mem, "cathode_left", t, x, j]
                            for j in self.ion_set
                        )
                        ** -1
                        * sum(
                            self.concentrate.properties[t, x].conc_mol_phase_comp[
                                "Liq", j
                            ]
                            for j in self.ion_set
                        )
                        ** -1
                    )
                )

        @self.Constraint(
            self.membrane_set,
            self.electrode_side_set,
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total non-ohmic potential across the two diffusion layers of an iem.",
        )
        def eq_dl_thickness(self, mem, side, t, x):
            if mem == "cem" and side == "cathode_left":
                return self.dl_thickness_x[mem, side, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * sum(
                        self.concentrate.properties[t, x].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * self.current_dens_lim_x[t, x] ** -1
                )
            elif mem == "cem" and side == "anode_right":
                return self.dl_thickness_x[mem, side, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * sum(
                        self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * self.current_dens_lim_x[t, x] ** -1
                )
            elif mem == "aem" and side == "cathode_left":
                return self.dl_thickness_x[mem, side, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * sum(
                        self.diluate.properties[t, x].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * self.current_dens_lim_x[t, x] ** -1
                )
            else:
                return self.dl_thickness_x[mem, side, t, x] == (
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
                            self.diluate.properties[t, x].trans_num_phase_comp["Liq", j]
                            / self.config.property_package.charge_comp[j]
                            for j in self.cation_set
                        )
                    )
                    ** -1
                    * sum(
                        self.concentrate.properties[t, x].conc_mol_phase_comp["Liq", j]
                        for j in self.cation_set
                    )
                    * self.current_dens_lim_x[t, x] ** -1
                )

    # Intialization routines
    def initialize_build(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        fail_on_warning=False,
        ignore_dof=False,
    ):
        """
        General wrapper for electrodialysis_1D initialization routines

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
        # Set the intial conditions over the 1D length from the state vars
        for k in blk.keys():
            for set in blk[k].diluate.properties:
                if ("flow_mol_phase_comp" or "flow_mass_phase_comp") not in blk[
                    k
                ].diluate.properties[set].define_state_vars():
                    raise ConfigurationError(
                        "Electrodialysis1D unit model requires "
                        "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                        "state variable basis to apply the 'propogate_initial_state' method"
                    )
                if "temperature" in blk[k].diluate.properties[set].define_state_vars():
                    blk[k].diluate.properties[set].temperature = value(
                        blk[k].diluate.properties[(0.0, 0.0)].temperature
                    )
                if "pressure" in blk[k].diluate.properties[set].define_state_vars():
                    blk[k].diluate.properties[set].pressure = value(
                        blk[k].diluate.properties[(0.0, 0.0)].pressure
                    )
                if (
                    "flow_mol_phase_comp"
                    in blk[k].diluate.properties[set].define_state_vars()
                ):
                    for ind in blk[k].diluate.properties[set].flow_mol_phase_comp:
                        blk[k].diluate.properties[set].flow_mol_phase_comp[ind] = value(
                            blk[k]
                            .diluate.properties[(0.0, 0.0)]
                            .flow_mol_phase_comp[ind]
                        )
                if (
                    "flow_mass_phase_comp"
                    in blk[k].diluate.properties[set].define_state_vars()
                ):
                    for ind in blk[k].diluate.properties[set].flow_mass_phase_comp:
                        blk[k].diluate.properties[set].flow_mass_phase_comp[
                            ind
                        ] = value(
                            blk[k]
                            .diluate.properties[(0.0, 0.0)]
                            .flow_mass_phase_comp[ind]
                        )
                if hasattr(blk[k], "conc_mem_surf_mol_x"):
                    for mem in blk[k].membrane_set:
                        for side in blk[k].electrode_side_set:
                            for j in blk[k].ion_set:
                                blk[k].conc_mem_surf_mol_x[mem, side, set, j].set_value(
                                    blk[k]
                                    .diluate.properties[set]
                                    .conc_mol_phase_comp["Liq", j]
                                )

            for set in blk[k].concentrate.properties:
                if ("flow_mol_phase_comp" or "flow_mass_phase_comp") not in blk[
                    k
                ].concentrate.properties[set].define_state_vars():
                    raise ConfigurationError(
                        "Electrodialysis1D unit model requires "
                        "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                        "state variable basis to apply the 'propogate_initial_state' method"
                    )
                if (
                    "temperature"
                    in blk[k].concentrate.properties[set].define_state_vars()
                ):
                    blk[k].concentrate.properties[set].temperature = value(
                        blk[k].concentrate.properties[(0.0, 0.0)].temperature
                    )
                if "pressure" in blk[k].concentrate.properties[set].define_state_vars():
                    blk[k].concentrate.properties[set].pressure = value(
                        blk[k].concentrate.properties[(0.0, 0.0)].pressure
                    )
                if (
                    "flow_mol_phase_comp"
                    in blk[k].concentrate.properties[set].define_state_vars()
                ):
                    for ind in blk[k].concentrate.properties[set].flow_mol_phase_comp:
                        blk[k].concentrate.properties[set].flow_mol_phase_comp[
                            ind
                        ] = value(
                            blk[k]
                            .concentrate.properties[(0.0, 0.0)]
                            .flow_mol_phase_comp[ind]
                        )
                if (
                    "flow_mass_phase_comp"
                    in blk[k].concentrate.properties[set].define_state_vars()
                ):
                    for ind in blk[k].concentrate.properties[set].flow_mass_phase_comp:
                        blk[k].concentrate.properties[set].flow_mass_phase_comp[
                            ind
                        ] = value(
                            blk[k]
                            .concentrate.properties[(0.0, 0.0)]
                            .flow_mass_phase_comp[ind]
                        )
                if hasattr(blk[k], "conc_mem_surf_mol_x"):
                    for mem in blk[k].membrane_set:
                        for side in blk[k].electrode_side_set:
                            for j in blk[k].ion_set:
                                blk[k].conc_mem_surf_mol_x[mem, side, set, j].set_value(
                                    blk[k]
                                    .concentrate.properties[set]
                                    .conc_mol_phase_comp["Liq", j]
                                )
                blk[k].total_areal_resistance_x[set].set_value(
                    (
                        blk[k].membrane_areal_resistance["aem"]
                        + blk[k].membrane_areal_resistance["cem"]
                        + blk[k].spacer_thickness
                        * (
                            blk[k].concentrate.properties[set].elec_cond_phase["Liq"]
                            ** -1
                            + blk[k].diluate.properties[set].elec_cond_phase["Liq"]
                            ** -1
                        )
                    )
                    * blk[k].cell_pair_num
                    + blk[k].electrodes_resistance
                )

        # ---------------------------------------------------------------------
        # Initialize diluate block
        flags_diluate = blk.diluate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        if not ignore_dof:
            check_dof(blk, fail_flag=fail_on_warning, logger=init_log)

        # Initialize concentrate block
        flags_concentrate = blk.concentrate.initialize(
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
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        check_solve(
            res,
            logger=init_log,
            fail_flag=fail_on_warning,
            checkpoint="Initialization Step 3",
        )
        # ---------------------------------------------------------------------
        # Release state
        blk.diluate.release_state(flags_diluate, outlvl)
        blk.concentrate.release_state(flags_concentrate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Scaling factors that user may setup
        # The users are highly encouraged to provide scaling factors for assessable vars below.
        # Not providing these vars will give a warning.
        if iscale.get_scaling_factor(self.cell_width, warning=True) is None:
            iscale.set_scaling_factor(self.cell_width, 1)
        if iscale.get_scaling_factor(self.cell_length, warning=True) is None:
            iscale.set_scaling_factor(self.cell_length, 1)
        if iscale.get_scaling_factor(self.cell_pair_num, warning=True) is None:
            iscale.set_scaling_factor(self.cell_pair_num, 0.1)
        if iscale.get_scaling_factor(self.membrane_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e4)

        if (
            iscale.get_scaling_factor(self.membrane_areal_resistance, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.membrane_areal_resistance, 1e4)

        if (
            iscale.get_scaling_factor(self.solute_diffusivity_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.solute_diffusivity_membrane, 1e10)

        if (
            iscale.get_scaling_factor(self.water_permeability_membrane, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.water_permeability_membrane, 1e4)

        if iscale.get_scaling_factor(self.spacer_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.spacer_thickness, 1e4)
        if iscale.get_scaling_factor(self.electrodes_resistance, warning=True) is None:
            iscale.set_scaling_factor(self.electrodes_resistance, 1e4)
        if hasattr(self, "voltage_applied") and (
            iscale.get_scaling_factor(self.voltage_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.voltage_applied, 1)
        if hasattr(self, "current_applied") and (
            iscale.get_scaling_factor(self.current_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.current_applied, 1)

        # For vars below, the users can choose but not required to provide scaling factors.
        # No warnings if no providing.

        for ind in self.total_areal_resistance_x:
            if (
                iscale.get_scaling_factor(
                    self.total_areal_resistance_x[ind], warning=False
                )
                is None
            ):
                sf = (
                    iscale.get_scaling_factor(self.membrane_areal_resistance) ** 2
                    + iscale.get_scaling_factor(self.spacer_thickness) ** 2
                    * (
                        iscale.get_scaling_factor(
                            self.diluate.properties[ind].elec_cond_phase["Liq"]
                        )
                        ** -2
                        + iscale.get_scaling_factor(
                            self.concentrate.properties[ind].elec_cond_phase["Liq"]
                        )
                        ** -2
                    )
                ) ** 0.5 * iscale.get_scaling_factor(self.cell_pair_num)
                iscale.set_scaling_factor(self.total_areal_resistance_x[ind], sf)

        for ind in self.current_density_x:
            if (
                iscale.get_scaling_factor(self.current_density_x[ind], warning=False)
                is None
            ):
                if (
                    self.config.operation_mode
                    == ElectricalOperationMode.Constant_Current
                ):
                    sf = (
                        iscale.get_scaling_factor(self.current_applied)
                        / iscale.get_scaling_factor(self.cell_width)
                        / iscale.get_scaling_factor(self.cell_length)
                    )
                    iscale.set_scaling_factor(self.current_density_x[ind], sf)
                else:
                    sf = iscale.get_scaling_factor(
                        self.voltage_applied
                    ) / iscale.get_scaling_factor(self.total_areal_resistance_x[ind])
                    iscale.set_scaling_factor(self.current_density_x[ind], sf)

        for ind in self.voltage_x:
            if iscale.get_scaling_factor(self.voltage_x[ind], warning=False) is None:
                sf = iscale.get_scaling_factor(
                    self.current_density_x[ind]
                ) * iscale.get_scaling_factor(self.total_areal_resistance_x[ind])
                iscale.set_scaling_factor(self.voltage_x[ind], sf)

        if hasattr(self, "conc_mem_surf_mol_x"):
            for ind in self.conc_mem_surf_mol_x:
                if iscale.get_scaling_factor(self.conc_mem_surf_mol_x[ind]) is None:
                    if (ind[0] == "cem" and ind[1] == "cathode_left") or (
                        ind[0] == "aem" and ind[1] == "anode_right"
                    ):
                        iscale.set_scaling_factor(
                            self.conc_mem_surf_mol_x[ind],
                            iscale.get_scaling_factor(
                                self.concentrate.properties[
                                    ind[2], ind[3]
                                ].conc_mol_phase_comp["Liq", ind[4]]
                            ),
                        )
                    else:
                        iscale.set_scaling_factor(
                            self.conc_mem_surf_mol_x[ind],
                            iscale.get_scaling_factor(
                                self.diluate.properties[
                                    ind[2], ind[3]
                                ].conc_mol_phase_comp["Liq", ind[4]]
                            ),
                        )

        if hasattr(self, "current_dens_lim_x"):
            if iscale.get_scaling_factor(self.current_dens_lim_x) is None:
                if (
                    self.config.limiting_current_density_method
                    == LimitingCurrentDensityMethod.InitialValue
                ):
                    for ind in self.current_dens_lim_x:
                        sf = (
                            self.config.limiting_current_density_data**-1
                            * sum(
                                iscale.get_scaling_factor(
                                    self.diluate.properties[
                                        ind[0], 0
                                    ].conc_mol_phase_comp["Liq", j]
                                )
                                ** 2
                                for j in self.cation_set
                            )
                            ** -0.5
                            * sum(
                                iscale.get_scaling_factor(
                                    self.diluate.properties[ind].conc_mol_phase_comp[
                                        "Liq", j
                                    ]
                                )
                                ** 2
                                for j in self.cation_set
                            )
                            ** 0.5
                        )
                        iscale.set_scaling_factor(self.current_dens_lim_x[ind], sf)

        if hasattr(self, "potential_nonohm_membrane_x"):
            if iscale.get_scaling_factor(self.potential_nonohm_membrane_x) is None:
                for ind in self.potential_nonohm_membrane_x:
                    sf = (
                        value(Constants.faraday_constant)
                        * value(Constants.gas_constant) ** -1
                        * 298.15**-1
                    )
                    iscale.set_scaling_factor(self.potential_nonohm_membrane_x[ind], sf)
        if hasattr(self, "potential_nonohm_dl_x"):
            if iscale.get_scaling_factor(self.potential_nonohm_dl_x) is None:
                for ind in self.potential_nonohm_dl_x:
                    sf = (
                        value(Constants.faraday_constant)
                        * value(Constants.gas_constant) ** -1
                        * 298.15**-1
                    )
                    iscale.set_scaling_factor(self.potential_nonohm_dl_x[ind], sf)
        if hasattr(self, "potential_ohm_dl_x"):
            if iscale.get_scaling_factor(self.potential_ohm_dl_x) is None:
                for ind in self.potential_nonohm_dl_x:
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
                    iscale.set_scaling_factor(self.potential_ohm_dl_x[ind], sf)
        if hasattr(self, "dl_thickness_x"):
            if iscale.get_scaling_factor(self.dl_thickness_x) is None:
                for ind in self.dl_thickness_x:
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
                            iscale.get_scaling_factor(self.conc_mem_surf_mol_x[ind, j])
                            ** 2
                            for j in self.cation_set
                        )
                        ** 0.5
                        * iscale.get_scaling_factor(
                            self.current_density_x[ind[2], ind[3]]
                        )
                        ** -1
                    )
                    iscale.set_scaling_factor(self.dl_thickness_x[ind], sf)

        for ind in self.diluate.power_electrical_x:
            if (
                iscale.get_scaling_factor(
                    self.diluate.power_electrical_x[ind], warning=False
                )
                is None
            ):
                iscale.set_scaling_factor(
                    self.diluate.power_electrical_x[ind],
                    iscale.get_scaling_factor(self.voltage_x[ind])
                    * iscale.get_scaling_factor(self.current_density_x[ind])
                    * iscale.get_scaling_factor(self.cell_width)
                    * iscale.get_scaling_factor(self.cell_length),
                )
        if (
            iscale.get_scaling_factor(self.specific_power_electrical, warning=False)
            is None
        ):
            iscale.set_scaling_factor(
                self.specific_power_electrical,
                3.6e6
                * iscale.get_scaling_factor(
                    self.diluate.power_electrical_x[
                        0, self.diluate.length_domain.last()
                    ]
                )
                * (
                    iscale.get_scaling_factor(
                        self.diluate.properties[
                            0, self.diluate.length_domain.last()
                        ].flow_vol_phase["Liq"]
                    )
                    * iscale.get_scaling_factor(self.cell_pair_num)
                )
                ** -1,
            )

        # Set up constraint scaling
        for ind, c in self.eq_get_total_areal_resistance_x.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.total_areal_resistance_x[ind])
            )

        for ind, c in self.eq_get_current_density.items():
            if self.config.operation_mode == ElectricalOperationMode.Constant_Current:
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.current_applied[ind[0]])
                )
            else:
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.voltage_applied[ind[0]])
                )
        for ind, c in self.eq_get_voltage_x.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.voltage_x[ind])
            )
        if hasattr(self, "eq_potential_nonohm_membrane_x"):
            for ind, c in self.eq_potential_nonohm_membrane_x.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.potential_nonohm_membrane_x[ind])
                )
        if hasattr(self, "eq_current_dens_lim_x"):
            for ind, c in self.eq_current_dens_lim_x.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.current_dens_lim_x[ind])
                )
        if hasattr(self, "eq_conc_polarization_ratio"):
            for ind, c in self.eq_conc_polarization_ratio.items():
                iscale.constraint_scaling_transform(c, 1)
        if hasattr(self, "eq_potential_nonohm_dl_x"):
            for ind, c in self.eq_potential_nonohm_dl_x.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.potential_nonohm_dl_x[ind])
                )
        if hasattr(self, "eq_potential_ohm_dl_x"):
            for ind, c in self.eq_potential_ohm_dl_x.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.potential_ohm_dl_x[ind])
                )
        if hasattr(self, "eq_dl_thickness_x"):
            for ind, c in self.eq_dl_thickness_x.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.dl_thickness_x[ind])
                )

        for ind, c in self.eq_mass_transfer_term_diluate.items():
            if ind[3] == "H2O":
                sf_osm = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * (
                        (
                            iscale.get_scaling_factor(
                                self.diluate.properties[
                                    (ind[0], ind[1])
                                ].pressure_osm_phase[ind[2]]
                            )
                        )
                        ** 2
                        + (
                            iscale.get_scaling_factor(
                                self.diluate.properties[
                                    (ind[0], ind[1])
                                ].pressure_osm_phase[ind[2]]
                            )
                        )
                        ** 2
                    )
                    ** 0.5
                )
                sf_eleosm = value(Constants.faraday_constant)
                iscale.constraint_scaling_transform(
                    c,
                    (sf_osm**2 + sf_eleosm**2) ** 0.5
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
            elif ind[3] in self.config.property_package.ion_set:
                sf_diff = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * (
                        (
                            (
                                iscale.get_scaling_factor(
                                    self.concentrate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                            + (
                                iscale.get_scaling_factor(
                                    self.diluate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                        )
                        ** 0.5
                    )
                )
                sf_elemig = value(Constants.faraday_constant)
                iscale.constraint_scaling_transform(
                    c,
                    (sf_diff**2 + sf_elemig**2) ** 0.5
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
            else:
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * (
                        (
                            (
                                iscale.get_scaling_factor(
                                    self.concentrate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                            + (
                                iscale.get_scaling_factor(
                                    self.diluate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                        )
                        ** 0.5
                    )
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
        for ind, c in self.eq_mass_transfer_term_concentrate.items():
            if ind[3] == "H2O":
                sf_osm = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * (
                        (
                            iscale.get_scaling_factor(
                                self.diluate.properties[
                                    (ind[0], ind[1])
                                ].pressure_osm_phase[ind[2]]
                            )
                        )
                        ** 2
                        + (
                            iscale.get_scaling_factor(
                                self.diluate.properties[
                                    (ind[0], ind[1])
                                ].pressure_osm_phase[ind[2]]
                            )
                        )
                        ** 2
                    )
                    ** 0.5
                )
                sf_eleosm = value(Constants.faraday_constant)
                iscale.constraint_scaling_transform(
                    c,
                    (sf_osm**2 + sf_eleosm**2) ** 0.5
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
            elif ind[3] in self.config.property_package.ion_set:
                sf_diff = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * (
                        (
                            (
                                iscale.get_scaling_factor(
                                    self.concentrate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                            + (
                                iscale.get_scaling_factor(
                                    self.diluate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                        )
                        ** 0.5
                    )
                )
                sf_elemig = value(Constants.faraday_constant)
                iscale.constraint_scaling_transform(
                    c,
                    (sf_diff**2 + sf_elemig**2) ** 0.5
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
            else:
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * (
                        (
                            (
                                iscale.get_scaling_factor(
                                    self.concentrate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                            + (
                                iscale.get_scaling_factor(
                                    self.diluate.properties[
                                        (ind[0], ind[1])
                                    ].conc_mol_phase_comp[ind[2], ind[3]]
                                )
                            )
                            ** 2
                        )
                        ** 0.5
                    )
                    * iscale.get_scaling_factor(self.cell_pair_num),
                )
        for ind, c in self.eq_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.diluate.power_electrical_x[ind])
            )

        for ind, c in self.eq_specific_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.specific_power_electrical)
            )

        for ind, c in self.eq_current_efficiency_x.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(self.current_density_x[ind[0], ind[1]])
                * iscale.get_scaling_factor(self.cell_width),
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
                "Total electrical power consumption(Watt)": self.diluate.power_electrical_x[
                    time_point, self.diluate.length_domain.last()
                ],
                "Specific electrical power consumption (kW*h/m**3)": self.specific_power_electrical[
                    time_point
                ],
                "Water recovery by mass": self.recovery_mass_H2O[time_point],
            },
            "exprs": {},
            "params": {},
        }

    def get_power_electrical(self, time_point=0):
        return self.diluate.power_electrical_x[
            time_point, self.diluate.length_domain.last()
        ]
