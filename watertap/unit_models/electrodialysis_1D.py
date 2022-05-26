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
    Constraint,
    NonNegativeReals,
    NonNegativeIntegers,
    Reals,
    Reference,
    value,
    units as pyunits,
)
from pyomo.dae import(
    ContinuousSet,
    DerivativeVar,
    Integral,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from sympy import Derivative

# Import Watertap cores
from watertap.core.util.initialization import check_solve, check_dof

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.control_volume1d import DistributedVars
from idaes.core.util.constants import Constants
from idaes.core.util.misc import add_object_reference
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Xiangyu Bi, Austin Ladshaw"

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
            default="Constant_Current",
            domain=In(["Constant_Current", "Constant_Voltage"]),
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

    # These config args are specifically for 1D control volumes
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

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit variables and parameters
        # # TODO: Add material props for membranes and such here
        # Create essential sets
        component_set = self.config.property_package.component_list
        ion_set=self.config.property_package.ion_set
        cation_set = self.config.property_package.cation_set
        anion_set = self.config.property_package.anion_set
        self.membrane_set = Set(initialize=["cem", "aem"]) #   cem = Cation-Exchange Membrane aem = Anion-Exchange Membrane
        length_set = ContinuousSet(bounds=(0, 1))
        
        # To require H2O must be in the component
        if "H2O" not in component_set:
            raise ConfigurationError(
                "Property Package MUST constain 'H2O' as a component"
            )

        # Create an ion_charge parameter for easier reference in model equations
        self.ion_charge = Param(
            ion_set,
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ion charge",
        )
        for j in component_set:
            if j in ion_set:
                self.ion_charge[j] = self.config.property_package.get_component(
                    j
                ).config.charge

        # Create unit model parameters and vars
        self.water_density = Param(
            initialize=1000,
            mutable=False,
            units=pyunits.kg * pyunits.m**-3,
            doc="density of water",
        )

        self.cell_pair_num = Var(
            initialize=1,
            domain=NonNegativeIntegers,
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
        self.membrane_surface_resistance = Var(
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
        self.current_applied = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp,
            doc="Current across a cell-pair or stack",
        )
       

        #self.i_x =DerivativeVar(self.current, wrt = length_set)

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
        '''
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
        # TODO: consider adding more performance as needed.
        '''
        # Build control volume for the diluate channel
        self.diluate = ControlVolume1DBlock(
            default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
                "area_definition": self.config.area_definition,
                "transformation_method": self.config.transformation_method,
                "transformation_scheme": self.config.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

        self.diluate.add_geometry()#length_domain = length_set 
        self.diluate.i_x_D = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp * pyunits.meter ** -2,
            doc="Current density at length = x",
        )
        #print("?????", self.diluate.length_domain.parent_block(),len(length_set))
        
        self.diluate.add_state_blocks(has_phase_equilibrium=False)

        self.diluate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        # # TODO: Temporarily disabling energy balances
        if hasattr(self.config, "energy_balance_type"):
            self.diluate.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )

        self.diluate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )

        # Apply transformation to the diluate channel
        self.diluate.apply_transformation()
        print("?????", len(length_set))
        self.first_element = self.diluate.length_domain.first()
        self.difference_elements = Set(
            ordered=True,
            initialize=(
                x for x in self.diluate.length_domain if x != self.first_element
            ),
        )
        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

        # Build control volume for concentrate side
        self.concentrate = ControlVolume1DBlock(
            default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
                "area_definition": self.config.area_definition,
                "transformation_method": self.config.transformation_method,
                "transformation_scheme": self.config.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

        self.concentrate.add_geometry() #length_domain_set=length_set
        self.concentrate.i_x_C = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp * pyunits.meter ** -2,
            doc="Current density at length = x",
        )
        self.concentrate.add_state_blocks(has_phase_equilibrium=False)

        self.concentrate.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        # # TODO: Temporarily disabling energy balances
        if hasattr(self.config, "energy_balance_type"):
            self.concentrate.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )

        self.concentrate.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )

        # Apply transformation to concentrate
        self.concentrate.apply_transformation()
        print("?????????????", self.diluate.length_domain.parent_block(),len(length_set),type(length_set),len(self.concentrate.length_domain))

        # Add ports (creates inlets and outlets for each channel)
        self.add_inlet_port(name="inlet_diluate", block=self.diluate)
        self.add_outlet_port(name="outlet_diluate", block=self.diluate)

        self.add_inlet_port(name="inlet_concentrate", block=self.concentrate)
        self.add_outlet_port(name="outlet_concentrate", block=self.concentrate)

        """
            Add references to the shared control volume length.
            This is the same length for the full electrodialysis cell,
            and is denoted as 'l' in the model description.
        """
        add_object_reference(self, "cell_normalized_length", self.diluate.length)

        #To report an Configuration Error if the two channels havae different lengths. 
        if self.diluate.length.value != self.concentrate.length.value:
            raise ConfigurationError("The Concentrate and Diluate channels"
            "must have the same lengths as the stack cell length")

        # -------- Add constraints ---------
        # # TODO: Add vars and associated constraints for all flux terms
        #           There will be 1 flux var for water and 1 flux var for all ions
        #           (vars can be indexed by species, so we only write it once)
        #
        #           Those vars will be coupled into the mass_transfer_term below
        #           and will be of opposite sign for each channel

        # Adds isothermal constraint if no energy balance present
        if not hasattr(self.config, "energy_balance_type"):

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Isothermal condition for diluate side",
            )
            def eq_isothermal_diluate(self, t, x):
                return (
                    self.diluate.properties[t, self.diluate.length_domain.first()].temperature
                        == self.diluate.properties[t, x].temperature
                )

        if not hasattr(self.config, "energy_balance_type"):

            @self.Constraint(
                self.flowsheet().time,
                self.difference_elements,
                doc="Isothermal condition for concentrate side",
            )
            def eq_isothermal_concentrate(self, t, x):

                return (
                    self.concentrate.properties[t, self.diluate.length_domain.first()].temperature
                    == self.concentrate.properties[t, x].temperature
                )

        # Add constraint for equal length of each channel
        
        @self.Constraint(doc="Constraint to ensure each channel has same length")

        def eq_equal_length(self):
            return self.diluate.length == self.concentrate.length
        
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            doc="electricity input condition for diluate channel",
        )
        def eq_elec_input_condition_D(self, t, x, p):
            if self.config.operation_mode == 'Constant_Current':
                print("!!!This is Constant_Current model")
                return self.diluate.i_x_D[t, x] * self.cell_width * self.cell_length == self.current_applied[t]
                
            else:
                print("!!!This is Constant_Voltage model")
                surface_resistance_cp = (
                self.membrane_surface_resistance["aem"]
                + self.membrane_surface_resistance["cem"]
                + self.spacer_thickness
                / (
                        self.concentrate.properties[t, x].electrical_conductivity_phase[p]
                        + self.diluate.properties[t, x].electrical_conductivity_phase[p]
                       
                    )       
                )
                return (
                    self.diluate.i_x_D[t, x]
                    * (
                        surface_resistance_cp * self.cell_pair_num
                        + self.electrodes_resistance
                    )
                    == self.voltage[t] 
                )

        @self.Constraint(
            self.flowsheet().time,
            self.concentrate.length_domain,
            self.config.property_package.phase_list,
            doc="electricity input condition for concentrate channel",
        )
        def eq_elec_input_condition_C(self, t, x, p):
            if self.config.operation_mode == 'Constant_Current':
                print("!!!This is Constant_Current model")
                return self.concentrate.i_x_C[t, x] * self.cell_width * self.cell_length == self.current_applied[t]
                
            else:
                print("!!!This is Constant_Voltage model")
                surface_resistance_cp = (
                self.membrane_surface_resistance["aem"]
                + self.membrane_surface_resistance["cem"]
                + self.spacer_thickness
                / (
                        self.concentrate.properties[t, x].electrical_conductivity_phase[p]
                        + self.diluate.properties[t, x].electrical_conductivity_phase[p]
                       
                    )       
                )
                return (
                    self.concentrate.i_x_C[t, x]
                    * (
                        surface_resistance_cp * self.cell_pair_num
                        + self.electrodes_resistance
                    )
                    == self.voltage[t] 
                )
        # Mass Transfer construct for the Diluate Channel
        
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term diluate channel",
        )
        def eq_mass_transfer_term_diluate(self, t, x, p, j):
            if x == self.diluate.length_domain.first():
                return Constraint.Skip 
            elif j == "H2O":
                return self.diluate.mass_transfer_term[t, x, p, j] == - self.cell_width * (
                    self.water_trans_number_membrane["cem"]
                    + self.water_trans_number_membrane["aem"]
                ) * (
                    self.diluate.i_x_D[t,x]
                    / Constants.faraday_constant
                ) - self.cell_width * self.water_density / self.config.property_package.mw_comp[j] * (
                    self.water_permeability_membrane["cem"]
                    + self.water_permeability_membrane["aem"]
                ) * (
                    self.concentrate.properties[t, x].pressure_osm_phase[p]
                    - self.diluate.properties[t, x].pressure_osm_phase[p]
                )
            elif j in self.config.property_package.ion_set:
                
                return self.diluate.mass_transfer_term[t, x, p, j] == - self.cell_width * (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.diluate.i_x_D[t, x]
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant 
                ) + self.cell_width * (
                    self.solute_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.solute_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x ].conc_mol_phase_comp[p, j]
                )
            else:
                return self.diluate.mass_transfer_term[t, x, p, j] == self.cell_width *  (
                    self.solute_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.solute_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x ].conc_mol_phase_comp[p, j]
                )
        
        # Mass Transfer construct for the Concentrate Channel
        @self.Constraint(
            self.flowsheet().config.time,
            self.concentrate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term concentrate channel",
        )
        def eq_mass_transfer_term_concentrate(self, t, x, p, j):
            if x == self.concentrate.length_domain.first:
                return Constraint.Skip     
            elif j == "H2O":
                return self.concentrate.mass_transfer_term[t, x, p, j] == self.cell_width * (
                    self.water_trans_number_membrane["cem"]
                    + self.water_trans_number_membrane["aem"]
                ) * (
                    self.concentrate.i_x_C[t,x]
                    / Constants.faraday_constant
                ) + self.cell_width * self.water_density / self.config.property_package.mw_comp[j] * (
                    self.water_permeability_membrane["cem"]
                    + self.water_permeability_membrane["aem"]
                ) * (
                    self.concentrate.properties[t, x].pressure_osm_phase[p]
                    - self.diluate.properties[t, x].pressure_osm_phase[p]
                )
            elif j in self.config.property_package.ion_set:
                return self.concentrate.mass_transfer_term[t, x, p, j] ==  self.cell_width * (
                    self.ion_trans_number_membrane["cem", j]
                    - self.ion_trans_number_membrane["aem", j]
                ) * (
                    self.current_utilization
                    * self.concentrate.i_x_C[t, x]
                ) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant 
                ) - self.cell_width * (
                    self.solute_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.solute_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x ].conc_mol_phase_comp[p, j]
                )
            else:
                return self.concentrate.mass_transfer_term[t, x, p, j] ==  -self.cell_width * (
                    self.solute_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                    + self.solute_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.concentrate.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x ].conc_mol_phase_comp[p, j]
                )
            


    # initialize method
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
        # Initialize diluate block
        flags = blk.diluate.initialize(
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

        # Initialize concentrate block
        flags = blk.concentrate.initialize(
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
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        check_solve(
            res,
            logger=init_log,
            fail_flag=fail_on_warning,
            checkpoint="Initialization Step 3",
        )

        # ---------------------------------------------------------------------
        # Release state
        blk.diluate.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        blk.concentrate.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
    '''

    def propogate_initial_state(self):

        i = 0
        for set in self.diluate.properties:

            # Should check 'define_state_vars' to see if user has provided
            #   state vars that are outside of the checks in this function
            if (
                "flow_mol_phase_comp"
                not in self.diluate.properties[set].define_state_vars()
            ):
                raise ConfigurationError(
                    "Electrodialysis1D unit model requires "
                    "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                    "state variable basis to apply the 'propogate_initial_state' method"
                )

            t = set[0]
            x = set[1]
            if i == 0:
                intial_set = (t, x)

            if "pressure" in self.diluate.properties[set].define_state_vars():
                self.diluate.properties[set].pressure = value(
                    self.diluate.properties[intial_set].pressure
                )
                self.concentrate.properties[set].pressure = value(
                    self.concentrate.properties[intial_set].pressure
                )

            if "temperature" in self.diluate.properties[set].define_state_vars():
                self.diluate.properties[set].temperature = value(
                    self.diluate.properties[intial_set].temperature
                )
                self.concentrate.properties[set].temperature = value(
                    self.concentrate.properties[intial_set].temperature
                )

            if (
                "flow_mol_phase_comp"
                in self.diluate.properties[set].define_state_vars()
            ):
                for ind in self.diluate.properties[set].flow_mol_phase_comp:
                    self.diluate.properties[set].flow_mol_phase_comp[ind] = value(
                        self.diluate.properties[intial_set].flow_mol_phase_comp[ind]
                    )
                    self.concentrate.properties[set].flow_mol_phase_comp[
                        ind
                    ] = value(
                        self.concentrate.properties[
                            intial_set
                        ].flow_mol_phase_comp[ind]
                    )

                # Check to see if 'flow_mass_phase_comp' is constructed
                if self.diluate.properties[set].is_property_constructed(
                    "flow_mass_phase_comp"
                ):
                    for ind in self.diluate.properties[set].flow_mass_phase_comp:
                        self.diluate.properties[set].flow_mass_phase_comp[
                            ind
                        ] = value(
                            self.diluate.properties[intial_set].flow_mol_phase_comp[
                                ind
                            ]
                        )
                        self.concentrate.properties[set].flow_mass_phase_comp[
                            ind
                        ] = value(
                            self.concentrate.properties[
                                intial_set
                            ].flow_mol_phase_comp[ind]
                        )

            i += 1
    '''
    '''
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
            sf = iscale.get_scaling_factor(
                self.membrane_thickness, default=1e3, warning=False
            )
            iscale.set_scaling_factor(self.membrane_thickness, sf)

        if iscale.get_scaling_factor(self.ion_diffusivity_membrane) is None:
            sf = iscale.get_scaling_factor(
                self.ion_diffusivity_membrane, default=1e9, warning=False
            )
            iscale.set_scaling_factor(self.ion_diffusivity_membrane, sf)

        if iscale.get_scaling_factor(self.water_permeability_membrane) is None:
            sf = iscale.get_scaling_factor(
                self.water_permeability_membrane, default=1e10, warning=False
            )
            iscale.set_scaling_factor(self.water_permeability_membrane, sf)

        if iscale.get_scaling_factor(self.diluate.area) is None:
            sf = iscale.get_scaling_factor(
                self.diluate.area, default=1, warning=False
            )
            iscale.set_scaling_factor(self.diluate.area, sf)

        if iscale.get_scaling_factor(self.concentrate.area) is None:
            sf = iscale.get_scaling_factor(
                self.concentrate.area, default=1, warning=False
            )
            iscale.set_scaling_factor(self.concentrate.area, sf)

        # Transform length constraint
        if (iscale.get_scaling_factor(self.diluate.length) is None) and (
            iscale.get_scaling_factor(self.concentrate.length) is None
        ):
            sf = iscale.get_scaling_factor(
                self.diluate.length, default=1, warning=False
            )
        elif (iscale.get_scaling_factor(self.diluate.length) is None) and (
            iscale.get_scaling_factor(self.concentrate.length) is not None
        ):
            sf = iscale.get_scaling_factor(
                self.concentrate.length, default=1, warning=False
            )
        else:
            sf = iscale.get_scaling_factor(
                self.diluate.length, default=1, warning=False
            )
        iscale.constraint_scaling_transform(self.eq_equal_length, sf)

        # Scaling factors for isothermal temperature
        if (
            iscale.get_scaling_factor(
                self.diluate.properties[0, self.first_element].temperature
            )
            is None
        ):
            sf = iscale.get_scaling_factor(
                self.diluate.properties[0, self.first_element].temperature,
                default=1e-2,
                warning=True,
            )
        if (
            iscale.get_scaling_factor(
                self.concentrate.properties[0, self.first_element].temperature
            )
            is None
        ):
            sf = iscale.get_scaling_factor(
                self.concentrate.properties[0, self.first_element].temperature,
                default=1e-2,
                warning=True,
            )

        sf = iscale.get_scaling_factor(
            self.concentrate.properties[0, self.first_element].temperature
        )
        for t in self.flowsheet().time:
            for x in self.difference_elements:
                iscale.constraint_scaling_transform(self.eq_isothermal_diluate[t, x], sf)
                iscale.constraint_scaling_transform(
                    self.eq_isothermal_concentrate[t, x], sf
                )

        # Iterate through and evaluate constraints to approximate scaling
        for p in self.config.property_package.phase_list:
            for j in self.config.property_package.component_list:
                for t in self.flowsheet().time:
                    # This index is tracking the first time through the x loop
                    ix = 0
                    for x in self.difference_elements:

                        # Reset scaling factors for _flow_terms
                        sf_flow = iscale.get_scaling_factor(
                            self.diluate._flow_terms[t, x, p, j]
                        )
                        if sf_flow >= 1000:
                            sf_flow = 1000
                        if sf_flow <= 0.001:
                            sf_flow = 0.001
                        iscale.set_scaling_factor(
                            self.diluate._flow_terms[t, x, p, j], sf_flow
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.diluate._flow_terms[t, 0, p, j], sf_flow
                            )

                        sf_flow = iscale.get_scaling_factor(
                            self.concentrate._flow_terms[t, x, p, j]
                        )
                        if sf_flow >= 1000:
                            sf_flow = 1000
                        if sf_flow <= 0.001:
                            sf_flow = 0.001
                        iscale.set_scaling_factor(
                            self.concentrate._flow_terms[t, x, p, j], sf_flow
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.concentrate._flow_terms[t, 0, p, j], sf_flow
                            )

                        ix += 1

                        # Calculations for nonelec_flux
                        sf1 = 0
                        if p == "Liq":
                            if j == "H2O":
                                Posm_D = self.diluate.properties[
                                    t, x
                                ].pressure_osm_phase[p]
                                Posm_C = self.concentrate.properties[
                                    t, x
                                ].pressure_osm_phase[p]
                                L_cem = pyunits.convert(
                                    self.water_permeability_membrane["cem"],
                                    to_units=units_meta("length")
                                    * units_meta("time") ** -1
                                    * units_meta("pressure") ** -1,
                                )
                                L_aem = pyunits.convert(
                                    self.water_permeability_membrane["aem"],
                                    to_units=units_meta("length")
                                    * units_meta("time") ** -1
                                    * units_meta("pressure") ** -1,
                                )
                                sf1 = (L_cem + L_aem) * (
                                    Posm_C
                                    * self.concentrate.properties[
                                        t, x
                                    ].conc_mol_phase_comp[p, j]
                                )
                            elif (
                                j in self.config.property_package.anion_set
                                or j in self.config.property_package.cation_set
                            ):
                                cem_rat = pyunits.convert(
                                    (
                                        self.ion_diffusivity_membrane["cem", j]
                                        / self.membrane_thickness["cem"]
                                    ),
                                    to_units=units_meta("length")
                                    * units_meta("time") ** -1,
                                )
                                aem_rat = pyunits.convert(
                                    (
                                        self.ion_diffusivity_membrane["aem", j]
                                        / self.membrane_thickness["aem"]
                                    ),
                                    to_units=units_meta("length")
                                    * units_meta("time") ** -1,
                                )
                                sf1 = -(cem_rat + aem_rat) * (
                                    self.concentrate.properties[
                                        t, x
                                    ].conc_mol_phase_comp[p, j]
                                )

                            else:
                                sf1 = 0.0
                        else:
                            sf1 = 0.0

                        if abs(value(sf1)) < 0.001:
                            sf1 = 0.001

                        iscale.set_scaling_factor(
                            self.nonelec_flux[t, x, p, j], 1 / abs(value(sf1))
                        )
                        iscale.constraint_scaling_transform(
                            self.eq_nonelec_flux[t, x, p, j], 1 / abs(value(sf1))
                        )

                        # Calculations for elec_flux
                        sf2 = 0
                        if p == "Liq":
                            if j == "H2O":
                                sf2 = 0
                            elif (
                                j in self.config.property_package.anion_set
                                or j in self.config.property_package.cation_set
                            ):
                                sf2 = 0
                            else:
                                sf2 = 0.0
                        else:
                            sf2 = 0

                        if abs(value(sf2)) < 0.001:
                            sf2 = 0.001

                        iscale.set_scaling_factor(
                            self.elec_flux[t, x, p, j], 1 / abs(value(sf2))
                        )
                        iscale.constraint_scaling_transform(
                            self.eq_elec_flux[t, x, p, j], 1 / abs(value(sf2))
                        )

                        # Calculations for mass transfer constraints
                        sf3 = 0
                        width = pyunits.convert(
                            self.cell_width, to_units=units_meta("length")
                        )
                        sf3 = sf1 * width + sf2 * width
                        if abs(value(sf3)) < 0.001:
                            sf3 = 0.001
                        iscale.constraint_scaling_transform(
                            self.eq_mass_transfer_term_diluate[t, x, p, j],
                            1 / abs(value(sf3)),
                        )
                        iscale.constraint_scaling_transform(
                            self.eq_mass_transfer_term_concentrate[t, x, p, j],
                            1 / abs(value(sf3)),
                        )

                    # End x
                # End t
            # End j
        # End p

        # Check for existence of true-to-apparent species info (comes from the Generic Property Package)
        for t in self.flowsheet().time:
            # This index is tracking the first time through the x loop
            ix = 0
            for x in self.difference_elements:

                # Constraints to scale
                if self.diluate.properties[t, x].is_property_constructed(
                    "true_to_appr_species"
                ):
                    for set in self.diluate.properties[t, x].true_to_appr_species:
                        p = set[0]
                        j = set[1]
                        iscale.constraint_scaling_transform(
                            self.diluate.properties[t, x].true_to_appr_species[
                                p, j
                            ],
                            1,
                        )
                        if ix == 0:
                            iscale.constraint_scaling_transform(
                                self.diluate.properties[t, 0].true_to_appr_species[
                                    p, j
                                ],
                                1,
                            )
                if self.concentrate.properties[t, x].is_property_constructed(
                    "true_to_appr_species"
                ):
                    for set in self.concentrate.properties[
                        t, x
                    ].true_to_appr_species:
                        p = set[0]
                        j = set[1]
                        iscale.constraint_scaling_transform(
                            self.concentrate.properties[t, x].true_to_appr_species[
                                p, j
                            ],
                            1,
                        )
                        if ix == 0:
                            iscale.constraint_scaling_transform(
                                self.concentrate.properties[
                                    t, 0
                                ].true_to_appr_species[p, j],
                                1,
                            )

                if self.diluate.properties[t, x].is_property_constructed(
                    "appr_mole_frac_constraint"
                ):
                    for set in self.diluate.properties[
                        t, x
                    ].appr_mole_frac_constraint:
                        p = set[0]
                        j = set[1]
                        iscale.constraint_scaling_transform(
                            self.diluate.properties[t, x].appr_mole_frac_constraint[
                                p, j
                            ],
                            1,
                        )
                        if ix == 0:
                            iscale.constraint_scaling_transform(
                                self.diluate.properties[
                                    t, 0
                                ].appr_mole_frac_constraint[p, j],
                                1,
                            )
                if self.concentrate.properties[t, x].is_property_constructed(
                    "appr_mole_frac_constraint"
                ):
                    for set in self.concentrate.properties[
                        t, x
                    ].appr_mole_frac_constraint:
                        p = set[0]
                        j = set[1]
                        iscale.constraint_scaling_transform(
                            self.concentrate.properties[
                                t, x
                            ].appr_mole_frac_constraint[p, j],
                            1,
                        )
                        if ix == 0:
                            iscale.constraint_scaling_transform(
                                self.concentrate.properties[
                                    t, 0
                                ].appr_mole_frac_constraint[p, j],
                                1,
                            )

                # Variables to scale
                if self.diluate.properties[t, x].is_property_constructed(
                    "flow_mol_phase_comp_apparent"
                ):
                    for set in self.diluate.properties[
                        t, x
                    ].flow_mol_phase_comp_apparent:
                        p = set[0]
                        j = set[1]
                        sf = iscale.get_scaling_factor(
                            self.diluate.properties[t, x].flow_mol_phase_comp[p, j]
                        )
                        iscale.set_scaling_factor(
                            self.diluate.properties[
                                t, x
                            ].flow_mol_phase_comp_apparent[p, j],
                            sf,
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.diluate.properties[
                                    t, 0
                                ].flow_mol_phase_comp_apparent[p, j],
                                sf,
                            )
                if self.concentrate.properties[t, x].is_property_constructed(
                    "flow_mol_phase_comp_apparent"
                ):
                    for set in self.concentrate.properties[
                        t, x
                    ].flow_mol_phase_comp_apparent:
                        p = set[0]
                        j = set[1]
                        sf = iscale.get_scaling_factor(
                            self.concentrate.properties[t, x].flow_mol_phase_comp[
                                p, j
                            ]
                        )
                        iscale.set_scaling_factor(
                            self.concentrate.properties[
                                t, x
                            ].flow_mol_phase_comp_apparent[p, j],
                            sf,
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.concentrate.properties[
                                    t, 0
                                ].flow_mol_phase_comp_apparent[p, j],
                                sf,
                            )

                if self.diluate.properties[t, x].is_property_constructed(
                    "mole_frac_phase_comp_apparent"
                ):
                    for set in self.diluate.properties[
                        t, x
                    ].mole_frac_phase_comp_apparent:
                        p = set[0]
                        j = set[1]
                        sf = iscale.get_scaling_factor(
                            self.diluate.properties[t, x].mole_frac_phase_comp[p, j]
                        )
                        iscale.set_scaling_factor(
                            self.diluate.properties[
                                t, x
                            ].mole_frac_phase_comp_apparent[p, j],
                            sf,
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.diluate.properties[
                                    t, 0
                                ].mole_frac_phase_comp_apparent[p, j],
                                sf,
                            )
                if self.concentrate.properties[t, x].is_property_constructed(
                    "mole_frac_phase_comp_apparent"
                ):
                    for set in self.concentrate.properties[
                        t, x
                    ].mole_frac_phase_comp_apparent:
                        p = set[0]
                        j = set[1]
                        sf = iscale.get_scaling_factor(
                            self.concentrate.properties[t, x].mole_frac_phase_comp[
                                p, j
                            ]
                        )
                        iscale.set_scaling_factor(
                            self.concentrate.properties[
                                t, x
                            ].mole_frac_phase_comp_apparent[p, j],
                            sf,
                        )
                        if ix == 0:
                            iscale.set_scaling_factor(
                                self.concentrate.properties[
                                    t, 0
                                ].mole_frac_phase_comp_apparent[p, j],
                                sf,
                            )

                ix += 1
            # End x
        # End t

        # # TODO: Add scaling factors
'''
    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Diluate Side Inlet": self.inlet_diluate,
                "Concentrate Side Inlet": self.inlet_concentrate,
                "Diluate Side Outlet": self.outlet_diluate,
                "Concentrate Side Outlet": self.outlet_concentrate,
            },
            time_point=time_point,
        )