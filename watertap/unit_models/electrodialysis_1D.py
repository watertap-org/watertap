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
from email.policy import default
from logging import warning
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
            default=4,
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

        self.diluate.add_geometry()
        
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

        '''
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
        '''
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
        #add_object_reference(self, "cell_normalized_length", self.diluate.length)
        self._make_performance()
    
    def _make_performance(self):
        # Create essential sets
        component_set = self.config.property_package.component_list
        ion_set=self.config.property_package.ion_set
        cation_set = self.config.property_package.cation_set
        anion_set = self.config.property_package.anion_set
        self.membrane_set = Set(initialize=["cem", "aem"]) #   cem = Cation-Exchange Membrane aem = Anion-Exchange Membrane
        #length_set = ContinuousSet(bounds=(0, 1))
        
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
            self.config.property_package.phase_list,
            initialize=1e-2,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.ohm * pyunits.meter**2,
            doc="Total areal resistance of a stack ",
        )
        self.current_applied = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp,
            doc="Current across a cell-pair or stack",
        )
        self.current_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp,
            doc="Current across a cell-pair or stack",
        )
        self.current_density_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.amp * pyunits.meter ** -2,
            doc="Current across a cell-pair or stack",
        )

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
            doc="Voltage across a stack, declared under the 'Constant Voltage' mode only",
        )
        self.current_utilization = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The current utilization including water electro-osmosis and ion diffusion",
        )
        # Performance metrics
        '''
        self.power_electrical = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 12100),
            domain=NonNegativeReals,
            units=pyunits.watt,
            doc="Electrical power consumption of a stack",
        )
        self.Dpower_electrical_Dx =DerivativeVar(self.power_electrical, wrt = self.diluate.length_domain, units=pyunits.watt * pyunits.meter **-1)
        
        self.specific_power_electrical = Var(
            self.flowsheet().time,
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.kW * pyunits.hour * pyunits.meter**-3,
            doc="Diluate-volume-flow-rate-specific electrical power consumption",
        )
         self.current_efficiency = Var(
            self.flowsheet().time,
            initialize=0.9,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The overall current efficiency for deionizaiton",
        )
        '''
        # TODO: consider adding more performance as needed.

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
                self.diluate.length_domain,
                doc="Isothermal condition for concentrate side",
            )
            def eq_isothermal_concentrate(self, t, x):

                return (
                    self.concentrate.properties[t, self.diluate.length_domain.first()].temperature
                    == self.concentrate.properties[t, x].temperature
                )

        # Add constraint for equal length of each channel
        '''
        @self.Constraint(doc="Constraint to ensure each channel has same length")
        def eq_equal_length(self):
            return self.diluate.length == self.concentrate.length
        '''
    
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            doc="Calculate the total area resistance of a stack"
        )
        def eq_get_total_areal_resistance_x(self, t, x, p):
            return self.total_areal_resistance_x[t, x, p] ==((
                self.membrane_areal_resistance["aem"]
                + self.membrane_areal_resistance["cem"]
                + self.spacer_thickness
                / (
                    self.concentrate.properties[t, x].electrical_conductivity_phase[p]
                    + self.diluate.properties[t, x].electrical_conductivity_phase[p]
                       
                    )       
                )
                * self.cell_pair_num
                + self.electrodes_resistance)
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            doc="calcualte current density from the electrical input"
        )
        def eq_get_current_density(self, t, x, p):
            if self.config.operation_mode == 'Constant_Current':
                return self.current_density_x[t, x] * self.cell_width * self.diluate.length == self.current_applied[t]
            else:
                return self.current_density_x[t, x] * self.total_areal_resistance_x[t, x, p] == self.voltage_applied[t]
    
        
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            doc="calcualte length_indexed voltage"
        )
        def eq_get_voltage_x(self, t, x, p):
            return self.voltage_x[t, x] == self.current_density_x[t, x] * self.total_areal_resistance_x[t, x, p]
        
        
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
                    self.current_density_x[t,x]
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
                    * self.current_density_x[t, x]
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
                    self.current_density_x[t,x]
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
                    * self.current_density_x[t, x]
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
        
        # Performance Metrics
        '''
        @self.Constraint(
            self.flowsheet().config.time,
            self.diluate.length_domain,
            doc="Electrical power consumption of a stack",
        )
        def eq_power_electrical(self, t, x):
            return self.power_electrical[t, x] == self.current_applied[t] * self.voltage_x[t, x]
        '''
        '''
        def eq_power_electrical(self, t, x):
            if self.config.operation_mode == 'Constant_Current':
                return self.power_electrical[t, x] == self.current_applied[t] * self.voltage_x[t, x]
            else:
                if x == self.diluate.length_domain.first():
                    return Constraint.Skip
                else: 
                    return self.Dpower_electrical_Dx[t, x] == self.voltage_applied[t] * self.current_density_x[t, x] * self.cell_width
        '''
        '''
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
                * self.diluate.properties[t, self.diluate.length_domain.last()].flow_vol_phase["Liq"]
                == self.current[t] * self.voltage[t]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Overall current efficiency evaluation",
        )
        def eq_current_efficiency(self, t):
            return (
                self.current_efficiency[t] * self.current[t]
                == sum(
                    self.diluate.properties[t, self.diluate.length_domain.first()].flow_mol_phase_comp["Liq", j]
                    * self.config.property_package.charge_comp[j]
                    - self.diluate.properties[t, self.diluate.length_domain.last()].flow_mol_phase_comp[
                        "Liq", j
                    ]
                    * self.config.property_package.charge_comp[j]
                    for j in self.config.property_package.cation_set
                )
                * Constants.faraday_constant
            )    
        '''

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
        
        for k in blk.keys():
            for set in blk[k].diluate.properties:
               # print('SET9', blk[k].diluate.properties.items().
                if "flow_mol_phase_comp" not in blk[k].diluate.properties[set].define_state_vars():
                    raise ConfigurationError(
                            "Electrodialysis1D unit model requires "
                            "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                            "state variable basis to apply the 'propogate_initial_state' method"
                        )
                if "teperature" in blk[k].diluate.properties[set].define_state_vars():
                    blk[k].diluate.properties[set].temperature = value(blk[k].diluate.properties[(0.0, 0.0)].temperature)
                if "pressure" in blk[k].diluate.properties[set].define_state_vars():
                    blk[k].diluate.properties[set].pressure = value(blk[k].diluate.properties[(0.0, 0.0)].pressure)
                if "flow_mol_phase_comp" in blk[k].diluate.properties[set].define_state_vars():
                    for ind in blk[k].diluate.properties[set].flow_mol_phase_comp:
                        blk[k].diluate.properties[set].flow_mol_phase_comp[ind] = value(blk[k].diluate.properties[(0.0, 0.0)].flow_mol_phase_comp[ind])
                        print('###INITINIT###',set, ind, value(blk[k].diluate.properties[set].flow_mol_phase_comp[ind]))
                if "flow_mass_phase_comp" in blk[k].diluate.properties[set].define_state_vars():
                    for ind in blk[k].diluate.properties[set].flow_mass_phase_comp:
                        blk[k].diluate.properties[set].flow_mass_phase_comp[ind] = value(blk[k].diluate.properties[(0.0, 0.0)].flow_mass_phase_comp[ind])

            for set in blk[k].concentrate.properties:
               # print('SET9', blk[k].diluate.properties.items().
                if "flow_mol_phase_comp" not in blk[k].concentrate.properties[set].define_state_vars():
                    raise ConfigurationError(
                            "Electrodialysis1D unit model requires "
                            "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                            "state variable basis to apply the 'propogate_initial_state' method"
                        )
                if "teperature" in blk[k].concentrate.properties[set].define_state_vars():
                    blk[k].concentrate.properties[set].temperature = value(blk[k].concentrate.properties[(0.0, 0.0)].temperature)
                if "pressure" in blk[k].diluate.properties[set].define_state_vars():
                    blk[k].concentrate.properties[set].pressure = value(blk[k].concentrate.properties[(0.0, 0.0)].pressure)
                if "flow_mol_phase_comp" in blk[k].concentrate.properties[set].define_state_vars():
                    for ind in blk[k].concentrate.properties[set].flow_mol_phase_comp:
                        blk[k].concentrate.properties[set].flow_mol_phase_comp[ind] = value(blk[k].concentrate.properties[(0.0, 0.0)].flow_mol_phase_comp[ind])
                if "flow_mass_phase_comp" in blk[k].concentrate.properties[set].define_state_vars():
                    for ind in blk[k].concentrate.properties[set].flow_mass_phase_comp:
                        blk[k].concentrate.properties[set].flow_mass_phase_comp[ind] = value(blk[k].concentrate.properties[(0.0, 0.0)].flow_mass_phase_comp[ind])
        
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
        for ind in blk.diluate._flow_terms_index:
            print(f'flow_term value of {ind} after init STEP2 is {blk.diluate._flow_terms[ind].value}')
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
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        blk.concentrate.release_state(flags_concentrate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
    

    
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Scaling factors that user may setup
        if iscale.get_scaling_factor(self.cell_width, warning=False) is None:
            iscale.set_scaling_factor(self.cell_width, 1)

        if iscale.get_scaling_factor(self.membrane_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e4)

        if iscale.get_scaling_factor(self.membrane_areal_resistance, warning=True) is None:
            iscale.set_scaling_factor(self.membrane_areal_resistance, 1e4)

        if iscale.get_scaling_factor(self.solute_diffusivity_membrane, warning=True) is None:
            iscale.set_scaling_factor(self.solute_diffusivity_membrane, 1e10)

        if iscale.get_scaling_factor(self.water_permeability_membrane, warning=True) is None:
            iscale.set_scaling_factor(self.water_permeability_membrane, 1e4)

        if iscale.get_scaling_factor(self.spacer_thickness, warning=True) is None:
            iscale.set_scaling_factor(self.spacer_thickness, 1e4)

        if iscale.get_scaling_factor(self.electrodes_resistance, warning=True) is None:
            iscale.set_scaling_factor(self.electrodes_resistance, 1e4)

        if iscale.get_scaling_factor(self.current_applied, warning=True) is None:
            iscale.set_scaling_factor(self.current_applied, 1)
        if iscale.get_scaling_factor(self.current_density_x, warning=False) is None:
            sf = iscale.get_scaling_factor(
                self.current_applied) / iscale.get_scaling_factor(self.cell_width)
            iscale.set_scaling_factor(self.current_density_x, sf)
        if iscale.get_scaling_factor(self.voltage_applied, warning=True) is None:
            iscale.set_scaling_factor(self.voltage_applied, 1)
              
        if iscale.get_scaling_factor(self.total_areal_resistance_x, warning=False) is None:
            sf = iscale.get_scaling_factor(
                self.membrane_areal_resistance
            ) #/ self.cell_pair_num
            iscale.set_scaling_factor(self.total_areal_resistance_x, sf)
        
        if iscale.get_scaling_factor(self.voltage_x, warning=False) is None:
            sf = iscale.get_scaling_factor(
                self.current_density_x) * iscale.get_scaling_factor(self.total_areal_resistance_x)
            iscale.set_scaling_factor(self.voltage_x, sf)

        # Transform length constraint
        '''
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
        '''
        '''
        '''
        # set the flow rate scaling factor constant throughout all finite elements. 
        for set in self.diluate.properties:
             if "flow_mol_phase_comp" in self.diluate.properties[set].define_state_vars():
                for ind in self.diluate.properties[set].flow_mol_phase_comp:
                    sf_flow = iscale.get_scaling_factor(self.diluate.properties[(0.0,0.0)].flow_mol_phase_comp[ind], default =1, warning = True)
                    iscale.set_scaling_factor(self.diluate.properties[set].flow_mol_phase_comp[ind], sf_flow)

        for set in self.concentrate.properties:
             if "flow_mol_phase_comp" in self.concentrate.properties[set].define_state_vars():
                for ind in self.concentrate.properties[set].flow_mol_phase_comp:
                    sf_flow = iscale.get_scaling_factor(self.concentrate.properties[(0.0,0.0)].flow_mol_phase_comp[ind], default =1, warning = True)
                    iscale.set_scaling_factor(self.concentrate.properties[set].flow_mol_phase_comp[ind], sf_flow)
                    
        #Transform constraint scaling 
        '''
       #for ind, c in self.eq_equal_length.items():
        #        iscale.constraint_scaling_transform(
         #           c, iscale.get_scaling_factor(self.diluate.length)
          #      )
        '''
        for ind, c in self.eq_get_total_areal_resistance_x.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.total_areal_resistance_x)
            )

        for ind, c in self.eq_get_current_density.items():
            if self.config.operation_mode == 'Constant_Current':
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.current_applied[ind[0]])
                )
            else:
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.voltage_applied[ind[0]])
                )
        for ind, c in self.eq_get_voltage_x.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.voltage_x)
            )

        for ind, c in self.eq_mass_transfer_term_diluate.items():
            if ind[3]== 'H2O':
                sf_osm = (
                        1e-3
                        * 0.018
                        * iscale.get_scaling_factor(self.water_permeability_membrane)
                        * ((iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].pressure_osm_phase[ind[2]]))**2 + (iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].pressure_osm_phase[ind[2]]))**2
                        )** 0.5
                    )
                sf_eleosm = Constants.faraday_constant
                iscale.constraint_scaling_transform(
                    c, (sf_osm ** 2 + sf_eleosm ** 2) ** 0.5
                )
            elif ind in self.config.property_package.ion_set:
                sf_diff = (
                        iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * (((iscale.get_scaling_factor(
                            self.concentrate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2 
                            +(iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2) ** 0.5)
                        )
                sf_elemig = Constants.faraday_constant
                iscale.constraint_scaling_transform(
                    c, (sf_diff ** 2 + sf_elemig ** 2) ** 0.5
                )
            else: 
                iscale.constraint_scaling_transform(
                    c,
                        iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * (((iscale.get_scaling_factor(
                            self.concentrate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2 
                            +(iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2) ** 0.5)
                        
                )
        for ind, c in self.eq_mass_transfer_term_concentrate.items():
            if ind[3]== 'H2O':
                sf_osm = (
                        1e-3
                        * 0.018
                        * iscale.get_scaling_factor(self.water_permeability_membrane)
                        * ((iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].pressure_osm_phase[ind[2]]))**2 + (iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].pressure_osm_phase[ind[2]]))**2
                        )** 0.5
                    )
                sf_eleosm = Constants.faraday_constant
                iscale.constraint_scaling_transform(
                    c, (sf_osm ** 2 + sf_eleosm ** 2) ** 0.5
                )
            elif ind in self.config.property_package.ion_set:
                sf_diff = (
                        iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * (((iscale.get_scaling_factor(
                            self.concentrate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2 
                            +(iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2) ** 0.5)
                        )
                sf_elemig = Constants.faraday_constant
                iscale.constraint_scaling_transform(
                    c, (sf_diff ** 2 + sf_elemig ** 2) ** 0.5
                )
            else: 
                iscale.constraint_scaling_transform(
                    c, 
                        iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                        / iscale.get_scaling_factor(self.membrane_thickness)
                        * (((iscale.get_scaling_factor(
                            self.concentrate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2 
                            +(iscale.get_scaling_factor(
                            self.diluate.properties[(ind[0],ind[1])
                            ].conc_mol_phase_comp[ind[2], ind[3]])) ** 2) ** 0.5)
                        
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
    