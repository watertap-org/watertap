#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import math

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    NonNegativeReals,
    value,
    Constraint,
    sqrt,
    units as pyunits,
)
from pyomo.dae import (
    DerivativeVar,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In

# Import Watertap cores
from watertap.core.util.initialization import check_solve, check_dof

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.misc import add_object_reference
from watertap.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.math import smooth_min

from idaes.core.util.exceptions import ConfigurationError, InitializationError

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.constants import Constants
from enum import Enum

from watertap.core import ControlVolume1DBlock, InitializationMixin
from watertap.costing.unit_models.electrodialysis import cost_electrodialysis

__author__ = "Johnson Dhanasekaran"

_log = idaeslog.getLogger(__name__)


class LimitingCurrentDensitybpmMethod(Enum):
    InitialValue = 0
    Empirical = 1


class ElectricalOperationMode(Enum):
    Constant_Current = 0
    Constant_Voltage = 1


class PressureDropMethod(Enum):
    none = 0
    experimental = 1
    Darcy_Weisbach = 2


class FrictionFactorMethod(Enum):
    fixed = 0
    Gurreri = 1
    Kuroda = 2


class HydraulicDiameterMethod(Enum):
    fixed = 0
    spacer_specific_area_known = 1
    conventional = 2


# Name of the unit model
@declare_process_block_class("Electrodialysis_Bipolar_1D")
class Electrodialysis_Bipolar_1DData(InitializationMixin, UnitModelBlockData):
    """
    1D Bipolar and Electrodialysis Model
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
        "pressure_drop_method",
        ConfigValue(
            default=PressureDropMethod.none,
            domain=In(PressureDropMethod),
            description="Method to calculate the frictional pressure drop in electrodialysis channels",
            doc="""
     **default** - ``PressureDropMethod.none``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``PressureDropMethod.none``", "The frictional pressure drop is neglected." 
           "``PressureDropMethod.experimental``", "The pressure drop is calculated by an experimental data as pressure drop per unit length."
           "``PressureDropMethod.Darcy_Weisbach``", "The pressure drop is calculated by the Darcy-Weisbach equation."
       """,
        ),
    )
    CONFIG.declare(
        "friction_factor_method",
        ConfigValue(
            default=FrictionFactorMethod.fixed,
            domain=In(FrictionFactorMethod),
            description="Method to calculate the Darcy's friction factor",
            doc="""
     **default** - ``FrictionFactorMethod.fixed``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``FrictionFactorMethod.fixed``", "Friction factor is fixed by users" 
           "``FrictionFactorMethod.Gurreri``", "Friction factor evaluated based on Gurreri's work"
           "``FrictionFactorMethod.Kuroda``", "Friction factor evaluated based on Kuroda's work"
       """,
        ),
    )

    CONFIG.declare(
        "hydraulic_diameter_method",
        ConfigValue(
            default=HydraulicDiameterMethod.conventional,
            domain=In(HydraulicDiameterMethod),
            description="Method to calculate the hydraulic diameter for a rectangular channel in ED",
            doc="""
     **default** - ``HydraulicDiameterMethod.conventional``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``HydraulicDiameterMethod.fixed``", "Hydraulic diameter is fixed by users" 
           "``HydraulicDiameterMethod.conventional``", "Conventional method for a rectangular channel with spacer porosity considered" 
           "``HydraulicDiameterMethod.spacer_specific_area_known``", "A method for spacer-filled channel requiring the spacer specific area data"
       """,
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
        "limiting_current_density_method_bpm",
        ConfigValue(
            default=LimitingCurrentDensitybpmMethod.InitialValue,
            domain=In(LimitingCurrentDensitybpmMethod),
            description="Configuration for method to compute the limiting current density across the bipolar membrane",
            doc="""
               **default** - ``LimitingCurrentDensitybpmMethod.InitialValue``

           .. csv-table::
               :header: "Configuration Options", "Description"

               "``LimitingCurrentDensitybpmMethod.InitialValue``", "Limiting current is calculated from a single initial value given by the user."
               "``LimitingCurrentDensitybpmMethod.Empirical``", "Limiting current density is calculated from the empirical equation"
           """,
        ),
    )

    CONFIG.declare(
        "salt_calculation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="""Salt calculation,
                    **default** - False.""",
        ),
    )

    CONFIG.declare(
        "limiting_current_density_bpm_data",
        ConfigValue(
            default=0.5,
            description="Limiting current density data input for bipolar membrane",
        ),
    )
    CONFIG.declare(
        "salt_input_cem",
        ConfigValue(
            default=100,
            description="Specified salt concentration on acid channel of the bipolar membrane",
        ),
    )
    CONFIG.declare(
        "salt_input_aem",
        ConfigValue(
            default=100,
            description="Specified salt concentration on base channel of the bipolar membrane",
        ),
    )

    CONFIG.declare(
        "limiting_potential_data",
        ConfigValue(
            default=0.5,
            description="Limiting potential data input",
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

    CONFIG.declare(
        "is_isothermal",
        ConfigValue(
            default=True,
            domain=Bool,
            description="""Assume isothermal conditions for control volume(s); energy_balance_type must be EnergyBalanceType.none,
    **default** - True.""",
        ),
    )

    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.none,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.none.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )

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
        "terms_fE",
        ConfigValue(
            default=40,
            domain=int,
            description="Number of terms in the second Wien electric field expression",
            doc="""Number of terms in the second Wien electric field expression. Highly recommended to use more than 15 (default=40)""",
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

    def _validate_config(self):
        if (
            self.config.is_isothermal
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            raise ConfigurationError(
                "If the isothermal assumption is used then the energy balance type must be none"
            )

    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()
        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Check configs for errors
        self._validate_config()

        # Create essential sets.
        self.membrane_set = Set(initialize=["cem", "aem", "bpm"])
        self.electrode_side = Set(initialize=["cathode_left", "anode_right"])
        add_object_reference(self, "ion_set", self.config.property_package.ion_set)

        add_object_reference(
            self, "cation_set", self.config.property_package.cation_set
        )
        add_object_reference(self, "anion_set", self.config.property_package.anion_set)
        add_object_reference(
            self, "component_set", self.config.property_package.component_list
        )
        # Create unit model parameters and vars

        self.cell_length = Var(
            initialize=0.5,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="The length of the bipolar electrodialysis cell",
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

        self.diluate.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False,
        )

        if self.config.is_isothermal:
            self.diluate.add_isothermal_assumption()

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

        # den_mass and visc_d in diluate and concentrate channels are the same
        add_object_reference(
            self, "dens_mass", self.diluate.properties[0, 0].dens_mass_phase["Liq"]
        )
        add_object_reference(
            self, "visc_d", self.diluate.properties[0, 0].visc_d_phase["Liq"]
        )

        # Apply the discretization transformation (Pyomo DAE) to the diluate block
        self.diluate.apply_transformation()

        # Build control volume for the base channel of the bipolar channel
        self.basic = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )
        self.basic.add_geometry(length_var=self.cell_length)
        self.basic.add_state_blocks(has_phase_equilibrium=False)
        self.basic.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.basic.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False,
        )

        if self.config.is_isothermal:
            self.basic.add_isothermal_assumption()

        self.basic.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )
        self.basic.apply_transformation()

        # Control volume for the acidic channel
        self.acidic = ControlVolume1DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )
        self.acidic.add_geometry(length_var=self.cell_length)
        self.acidic.add_state_blocks(has_phase_equilibrium=False)
        self.acidic.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.acidic.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False,
        )

        if self.config.is_isothermal:
            self.acidic.add_isothermal_assumption()

        self.acidic.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )
        self.acidic.apply_transformation()

        # Add ports (creates inlets and outlets for each channel)
        self.add_inlet_port(name="inlet_diluate", block=self.diluate)
        self.add_outlet_port(name="outlet_diluate", block=self.diluate)
        self.add_inlet_port(name="inlet_basic", block=self.basic)
        self.add_outlet_port(name="outlet_basic", block=self.basic)
        self.add_inlet_port(name="inlet_acidic", block=self.acidic)
        self.add_outlet_port(name="outlet_acidic", block=self.acidic)

        self.water_density = Param(
            initialize=1000,
            units=pyunits.kg * pyunits.m**-3,
            doc="density of water",
        )

        self.cell_triplet_num = Var(
            initialize=1,
            domain=NonNegativeReals,
            bounds=(1, 10000),
            units=pyunits.dimensionless,
            doc="cell triplet number in a stack",
        )

        # bipolar electrodialysis cell dimensional properties
        self.cell_width = Var(
            initialize=0.1,
            bounds=(1e-3, 1e2),
            units=pyunits.meter,
            doc="The width of the electrodialysis cell, denoted as b in the model description",
        )
        self.channel_height = Var(
            initialize=0.0001,
            units=pyunits.meter,
            doc="The distance between the consecutive membranes",
        )
        self.spacer_porosity = Var(
            initialize=0.7,
            bounds=(0.01, 1),
            units=pyunits.dimensionless,
            doc='The porosity of spacer in the BMED channels. This is also referred to elsewhere as "void fraction" or "volume parameters"',
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
            initialize=0.5,
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
        self.total_areal_resistance_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=2e-4,
            bounds=(0, 1e3),
            units=pyunits.ohm * pyunits.meter**2,
            doc="Total areal resistance of a stack ",
        )

        self.membrane_areal_resistance_coef_0 = Var(
            initialize=2e-4,
            bounds=(1e-6, 1),
            units=pyunits.ohm * pyunits.meter**2,
            doc="Constant areal resistance of membrane at infinity-approximated electrolyte concentration",
        )
        self.membrane_areal_resistance_coef_1 = Var(
            initialize=0,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.ohm * pyunits.kg * pyunits.m**-1,
            doc="Coefficient of membrane areal resistance to 1/c, where c is the electrolyte concentration",
        )
        if self.config.operation_mode == ElectricalOperationMode.Constant_Current:
            self.current_applied = Var(
                self.flowsheet().time,
                initialize=1,
                bounds=(0, 1000),
                units=pyunits.amp,
                doc="Current across a cell-triplet or stack, declared under the 'Constant Current' mode only",
            )

        self.current_density_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1e6),
            units=pyunits.amp * pyunits.meter**-2,
            doc="Current density across the membrane as a function of the normalized length",
        )
        self.voltage_membrane_drop = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            bounds=(0, 1000),
            units=pyunits.volt,
            doc="Potential drop across the bipolar membrane",
        )

        if self.config.operation_mode == ElectricalOperationMode.Constant_Voltage:
            self.voltage_applied = Var(
                self.flowsheet().time,
                initialize=100,
                bounds=(0, 2000 * 1e3),
                units=pyunits.volt,
                doc="Voltage across a stack, declared under the 'Constant Voltage' mode only",
            )

        self.voltage_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=100,
            bounds=(0, 2000 * 1e3),
            units=pyunits.volt,
            doc="Voltage across a stack",
        )

        self.electrodes_resistance = Var(
            initialize=0,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.ohm * pyunits.meter**2,
            doc="areal resistance of TWO electrode compartments of a stack",
        )
        self.current_utilization = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The current utilization including water electro-osmosis and ion diffusion",
        )
        self.shadow_factor = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="The reduction in area due to limited area available for flow",
        )

        # Performance metrics
        self.specific_power_electrical = Var(
            self.flowsheet().time,
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.kW * pyunits.hour * pyunits.meter**-3,
            doc="Diluate-volume-flow-rate-specific electrical power consumption",
        )
        self.velocity_diluate = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,
            units=pyunits.meter * pyunits.second**-1,
            doc="Linear velocity of flow in the diluate",
        )
        self.velocity_basic = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,
            units=pyunits.meter * pyunits.second**-1,
            doc="Linear velocity of flow in the base channel of the bipolar membrane",
        )
        self.velocity_acidic = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=0.01,
            units=pyunits.meter * pyunits.second**-1,
            doc="Linear velocity of flow in the acid channel of the bipolar membrane",
        )
        self.elec_field_non_dim = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            units=pyunits.dimensionless,
            doc="Limiting current density across the bipolar membrane as a function of the normalized length",
        )
        self.relative_permittivity = Var(
            initialize=30,
            bounds=(1, 80),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Relative permittivity",
        )
        self.membrane_fixed_charge = Var(
            initialize=1.5e3,
            bounds=(1e-1, 1e5),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Membrane fixed charge",
        )
        self.k2_zero = Var(
            initialize=2 * 10**-5,
            bounds=(1e-10, 1e2),
            units=pyunits.second**-1,
            doc="Dissociation rate constant at no electric field",
        )
        self.salt_conc_ael_ref = Var(
            initialize=1e3,
            bounds=(1e-8, 1e6),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Fixed salt concentration on the base channel of the bipolar membrane",
        )
        self.salt_conc_cel_ref = Var(
            initialize=1e3,
            bounds=(1e-8, 1e6),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Fixed salt concentration on the acid channel of the bipolar membrane",
        )
        self.salt_conc_dilu_ref = Var(
            initialize=1e3,
            bounds=(1e-8, 1e6),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Fixed salt concentration on the diluate channel ",
        )
        self.salt_conc_ael_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1e3,
            bounds=(1e-8, 1e6),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Salt concentration on the base channel of the bipolar membrane",
        )
        self.salt_conc_cel_x = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1e3,
            bounds=(1e-6, 1e4),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Salt concentration on the acid channel of the bipolar membrane",
        )
        self.current_dens_lim_bpm = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1e2,
            bounds=(0, 1e5),
            units=pyunits.amp * pyunits.meter**-2,
            doc="Limiting current density across the bipolar membrane",
        )

        self.diffus_mass = Var(
            initialize=2e-9,
            bounds=(1e-16, 1e-6),
            units=pyunits.meter**2 * pyunits.second**-1,
            doc="The mass diffusivity of the solute as molecules (not individual ions)",
        )
        self.conc_water = Var(
            initialize=55 * 1e3,
            bounds=(1e-2, 1e6),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Concentration of water within the channel",
        )

        # Fluxes Vars for constructing mass transfer terms
        self.generation_cel_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component generated by water splitting on the acid channel of the bipolar membrane",
        )
        self.generation_ael_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component generated by water splitting on the base channel of the bipolar membrane",
        )
        self.elec_migration_bpm_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by electrical migration across the bipolar membrane",
        )
        self.nonelec_bpm_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by non-electrical forces across the bipolar membrane",
        )

        self.elec_migration_mono_cem_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by electrical migration across the monopolar CEM membrane",
        )
        self.nonelec_mono_cem_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by non-electrical forces across the monopolar CEM membrane",
        )

        self.elec_migration_mono_aem_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by electrical migration across the monopolar AEM membrane",
        )
        self.nonelec_mono_aem_flux = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Molar flux_in of a component across the membrane driven by non-electrical forces across the monopolar AEM membrane",
        )

        # extension options
        self._make_catalyst()

        if (
            not self.config.pressure_drop_method == PressureDropMethod.none
        ) and self.config.has_pressure_change:
            self._pressure_drop_calculation()

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Pressure drop expression as calculated by the pressure drop data, diluate.",
            )
            def eq_deltaP_diluate(self, t, x):
                return self.diluate.deltaP[t, x] == -self.pressure_drop[t]

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Pressure drop expression as calculated by the pressure drop data, "
                "base channel of the bipolar membrane.",
            )
            def eq_deltaP_basic(self, t, x):
                return self.basic.deltaP[t, x] == -self.pressure_drop[t]

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Pressure drop expression as calculated by the pressure drop data,"
                "  acid channel of the bipolar membrane.",
            )
            def eq_deltaP_acidic(self, t, x):
                return self.acidic.deltaP[t, x] == -self.pressure_drop[t]

        elif self.config.pressure_drop_method == PressureDropMethod.none and (
            not self.config.has_pressure_change
        ):
            pass
        else:
            raise ConfigurationError(
                "A valid (not none) pressure_drop_method and has_pressure_change being True "
                "must be both used or unused at the same time. "
            )

        # To require H2O must be in the component
        if "H2O" not in self.component_set:
            raise ConfigurationError(
                "Property Package MUST constain 'H2O' as a component"
            )

        # Build Constraints
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate flow velocity in a single diluate channel, based on the average of inlet and outlet",
        )
        def eq_get_velocity_diluate(self, t, x):
            return (
                self.velocity_diluate[t, x]
                * self.cell_width
                * self.shadow_factor
                * self.channel_height
                * self.spacer_porosity
                * self.cell_triplet_num
                == self.diluate.properties[t, x].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate flow velocity in a single base channel of the bipolar membrane channel,"
            " based on the average of inlet and outlet",
        )
        def eq_get_velocity_basic(self, t, x):
            return (
                self.velocity_basic[t, x]
                * self.cell_width
                * self.shadow_factor
                * self.channel_height
                * self.spacer_porosity
                * self.cell_triplet_num
                == self.basic.properties[t, x].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate flow velocity in a single acid channel of the bipolar membrane channel,"
            " based on the average of inlet and outlet",
        )
        def eq_get_velocity_acidic(self, t, x):
            return (
                self.velocity_acidic[t, x]
                * self.cell_width
                * self.shadow_factor
                * self.channel_height
                * self.spacer_porosity
                * self.cell_triplet_num
                == self.acidic.properties[t, x].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Evaluate salt concentration on AEM side of the bipolar membrane",
        )
        def eq_salt_aem(self, t, x):
            if self.config.salt_calculation:
                conc_unit = 1 * pyunits.mole * pyunits.meter**-3

                return (
                    self.salt_conc_ael_x[t, x]
                    == smooth_min(
                        self.basic.properties[t, x].conc_mol_phase_comp["Liq", "Na_+"]
                        / conc_unit,
                        self.basic.properties[t, x].conc_mol_phase_comp["Liq", "Cl_-"]
                        / conc_unit,
                    )
                    * conc_unit
                )
            else:
                return self.salt_conc_ael_x[t, x] == self.salt_conc_ael_ref

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Evaluate salt concentration on CEM side of the bipolar membrane",
        )
        def eq_salt_cem(self, t, x):

            if self.config.salt_calculation:
                conc_unit = 1 * pyunits.mole * pyunits.meter**-3

                return (
                    self.salt_conc_cel_x[t, x]
                    == smooth_min(
                        self.acidic.properties[t, x].conc_mol_phase_comp["Liq", "Na_+"]
                        / conc_unit,
                        self.acidic.properties[t, x].conc_mol_phase_comp["Liq", "Cl_-"]
                        / conc_unit,
                    )
                    * conc_unit
                )
            else:
                return self.salt_conc_cel_x[t, x] == self.salt_conc_cel_ref

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate limiting current density across the bipolar membrane",
        )
        def eq_current_dens_lim_bpm(self, t, x):
            if (
                self.config.limiting_current_density_method_bpm
                == LimitingCurrentDensitybpmMethod.InitialValue
            ):
                return self.current_dens_lim_bpm[t, x] == (
                    self.config.limiting_current_density_bpm_data
                    * pyunits.amp
                    * pyunits.meter**-2
                )
            elif (
                self.config.limiting_current_density_method_bpm
                == LimitingCurrentDensitybpmMethod.Empirical
            ):
                return self.current_dens_lim_bpm[
                    t, x
                ] == self.diffus_mass * Constants.faraday_constant * (
                    (self.salt_conc_ael_x[t, x] + self.salt_conc_cel_x[t, x]) * 0.5
                ) ** 2 / (
                    self.membrane_thickness["bpm"] * self.membrane_fixed_charge
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the total areal resistance of a stack",
        )
        def eq_get_total_areal_resistance_x(self, t, x):
            return self.total_areal_resistance_x[t, x] == (
                (
                    self.membrane_areal_resistance_coef_1
                    / (
                        self.acidic.properties[t, x].conc_mass_phase_comp["Liq", "H_+"]
                        + self.acidic.properties[t, x].conc_mass_phase_comp[
                            "Liq", "Cl_-"
                        ]
                        + self.basic.properties[t, x].conc_mass_phase_comp[
                            "Liq", "Na_+"
                        ]
                        + self.basic.properties[t, x].conc_mass_phase_comp[
                            "Liq", "OH_-"
                        ]
                    )
                    + self.membrane_areal_resistance_coef_0
                    + self.channel_height
                    * (
                        self.acidic.properties[t, x].elec_cond_phase["Liq"] ** -1
                        + self.basic.properties[t, x].elec_cond_phase["Liq"] ** -1
                        + self.diluate.properties[t, x].elec_cond_phase["Liq"] ** -1
                    )
                )
                * self.cell_triplet_num
                + self.electrodes_resistance
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate current density from the electrical input",
        )
        def eq_get_current_density(self, t, x):

            @self.Constraint(
                self.flowsheet().time,
                self.diluate.length_domain,
                doc="Calculate total current generated via catalyst action",
            )
            def eq_current_relationship(self, t, x):
                return self.current_density_x[t, x] == (
                    self.current_dens_lim_bpm[t, x]
                    + self.flux_splitting[t, x] * Constants.faraday_constant
                )

            if self.config.operation_mode == ElectricalOperationMode.Constant_Current:
                return (
                    self.current_density_x[t, x]
                    * self.cell_width
                    * self.shadow_factor
                    * self.diluate.length
                    == self.current_applied[t]
                )
            else:
                return (
                    self.current_density_x[t, x] * self.total_areal_resistance_x[t, x]
                    + self.voltage_membrane_drop[t, x] * self.cell_triplet_num
                    == self.voltage_applied[t]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="calculate length_indexed voltage",
        )
        def eq_get_voltage_x(self, t, x):
            return (
                self.current_density_x[t, x] * self.total_areal_resistance_x[t, x]
                + self.voltage_membrane_drop[t, x] * self.cell_triplet_num
                == self.voltage_x[t, x]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for water splitting acid channel of bipolar membrane flux_in",
        )
        def eq_generation_cel_flux(self, t, x, p, j):
            if j == "H_+":
                return self.generation_cel_flux[t, x, p, j] == self.flux_splitting[t, x]

            else:
                if j == "H2O":
                    return (
                        self.generation_cel_flux[t, x, p, j]
                        == -0.5 * self.flux_splitting[t, x]
                    )

                else:
                    self.generation_cel_flux[t, x, p, j].fix(0)
                    return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for water splitting base channel of bipolar membrane flux_in",
        )
        def eq_generation_ael_flux(self, t, x, p, j):
            if j == "OH_-":
                return self.generation_ael_flux[t, x, p, j] == self.flux_splitting[t, x]

            else:
                if j == "H2O":
                    return (
                        self.generation_ael_flux[t, x, p, j]
                        == -0.5 * self.flux_splitting[t, x]
                    )
                else:
                    self.generation_ael_flux[t, x, p, j].fix(0)
                    return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for electrical migration across the monopolar CEM flux_in",
        )
        def eq_elec_migration_mono_cem(self, t, x, p, j):
            if j == "H2O":
                return self.elec_migration_mono_cem_flux[t, x, p, j] == (
                    self.water_trans_number_membrane["cem"]
                ) * (self.current_density_x[t, x] / Constants.faraday_constant)
            elif j in self.ion_set:
                return self.elec_migration_mono_cem_flux[t, x, p, j] == (
                    self.ion_trans_number_membrane["cem", j]
                ) * (self.current_utilization * self.current_density_x[t, x]) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                self.elec_migration_mono_cem_flux[t, x, p, j].fix(0)
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for electrical migration across the monopolar AEM flux_in",
        )
        def eq_elec_migration_mono_aem_flux(self, t, x, p, j):
            if j == "H2O":
                return self.elec_migration_mono_aem_flux[t, x, p, j] == (
                    self.water_trans_number_membrane["aem"]
                ) * (self.current_density_x[t, x] / Constants.faraday_constant)
            elif j in self.ion_set:
                return self.elec_migration_mono_aem_flux[t, x, p, j] == (
                    -self.ion_trans_number_membrane["aem", j]
                ) * (self.current_utilization * self.current_density_x[t, x]) / (
                    self.config.property_package.charge_comp[j]
                    * Constants.faraday_constant
                )
            else:
                self.elec_migration_mono_aem_flux[t, x, p, j].fix(0)
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for electrical migration across the bipolar membrane flux_in",
        )
        def eq_elec_migration_bpm_flux(self, t, x, p, j):
            if j == "H2O":
                return self.elec_migration_bpm_flux[t, x, p, j] == (
                    self.water_trans_number_membrane["bpm"]
                ) * (self.current_density_x[t, x] / Constants.faraday_constant)

            elif j in self.ion_set:
                if not (j == "H_+" or j == "OH_-"):
                    return self.elec_migration_bpm_flux[t, x, p, j] == (
                        self.ion_trans_number_membrane["bpm", j]
                    ) * (self.current_utilization * self.current_dens_lim_bpm[t, x]) / (
                        self.config.property_package.charge_comp[j]
                        * Constants.faraday_constant
                    )
                else:

                    self.elec_migration_bpm_flux[t, x, p, j].fix(
                        0 * pyunits.mol * pyunits.m**-2 * pyunits.s**-1
                    )
                    return Constraint.Skip
            else:
                self.elec_migration_bpm_flux[t, x, p, j].fix(
                    0 * pyunits.mol * pyunits.m**-2 * pyunits.s**-1
                )
                return Constraint.Skip

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux across the monopolar CEM flux_in",
        )
        def eq_nonelec_mono_cem_flux(self, t, x, p, j):
            if j == "H2O":
                return self.nonelec_mono_cem_flux[
                    t, x, p, j
                ] == self.water_density / self.config.property_package.mw_comp[j] * (
                    self.water_permeability_membrane["cem"]
                ) * (
                    self.basic.properties[t, x].pressure_osm_phase[p]
                    - self.diluate.properties[t, x].pressure_osm_phase[p]
                )

            else:
                return self.nonelec_mono_cem_flux[t, x, p, j] == -(
                    self.solute_diffusivity_membrane["cem", j]
                    / self.membrane_thickness["cem"]
                ) * (
                    self.basic.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x].conc_mol_phase_comp[p, j]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux across the monopolar AEM flux_in",
        )
        def eq_nonelec_mono_aem_flux(self, t, x, p, j):
            if j == "H2O":
                return self.nonelec_mono_aem_flux[
                    t, x, p, j
                ] == self.water_density / self.config.property_package.mw_comp[j] * (
                    self.water_permeability_membrane["aem"]
                ) * (
                    self.acidic.properties[t, x].pressure_osm_phase[p]
                    - self.diluate.properties[t, x].pressure_osm_phase[p]
                )

            else:
                return self.nonelec_mono_aem_flux[t, x, p, j] == -(
                    self.solute_diffusivity_membrane["aem", j]
                    / self.membrane_thickness["aem"]
                ) * (
                    self.acidic.properties[t, x].conc_mol_phase_comp[p, j]
                    - self.diluate.properties[t, x].conc_mol_phase_comp[p, j]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Equation for non-electrical flux across the bipolar membrane flux_in",
        )
        def eq_nonelec_bpm_flux(self, t, x, p, j):
            if j == "H2O":
                return self.nonelec_bpm_flux[
                    t, x, p, j
                ] == self.water_density / self.config.property_package.mw_comp[j] * (
                    self.water_permeability_membrane["bpm"]
                ) * (
                    self.basic.properties[t, x].pressure_osm_phase[p]
                    - self.acidic.properties[t, x].pressure_osm_phase[p]
                )

            else:
                self.nonelec_bpm_flux[t, x, p, j].fix(0)
                return Constraint.Skip

        # Add constraints for mass transfer terms (diluate)
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the diluate channel",
        )
        def eq_mass_transfer_term_diluate(self, t, x, p, j):
            return (
                self.diluate.mass_transfer_term[t, x, p, j]
                == -(
                    self.elec_migration_mono_aem_flux[t, x, p, j]
                    + self.elec_migration_mono_cem_flux[t, x, p, j]
                    + self.nonelec_mono_aem_flux[t, x, p, j]
                    + self.nonelec_mono_cem_flux[t, x, p, j]
                )
                * (self.cell_width * self.shadow_factor)
                * self.cell_triplet_num
            )

        # Add constraints for mass transfer terms (base channel of the bipolar membrane)
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the base channel of the bipolar membrane",
        )
        def eq_mass_transfer_term_basic(self, t, x, p, j):
            return (
                self.basic.mass_transfer_term[t, x, p, j]
                == (
                    self.generation_ael_flux[t, x, p, j]
                    - self.elec_migration_bpm_flux[t, x, p, j]
                    - self.nonelec_bpm_flux[t, x, p, j]
                    + self.elec_migration_mono_cem_flux[t, x, p, j]
                    + self.nonelec_mono_cem_flux[t, x, p, j]
                )
                * (self.cell_width * self.shadow_factor)
                * self.cell_triplet_num
            )

        # Add constraints for mass transfer terms (acid channel of the bipolar membrane)
        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term for the acid channel of the bipolar membrane channel",
        )
        def eq_mass_transfer_term_acidic(self, t, x, p, j):
            return (
                self.acidic.mass_transfer_term[t, x, p, j]
                == (
                    self.generation_cel_flux[t, x, p, j]
                    + self.elec_migration_bpm_flux[t, x, p, j]
                    + self.nonelec_bpm_flux[t, x, p, j]
                    + self.elec_migration_mono_aem_flux[t, x, p, j]
                    + self.nonelec_mono_aem_flux[t, x, p, j]
                )
                * (self.cell_width * self.shadow_factor)
                * self.cell_triplet_num
            )

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
                    * self.shadow_factor
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

    def _make_catalyst(self):

        self.flux_splitting = Var(
            self.flowsheet().time,
            self.diluate.length_domain,
            initialize=1,
            domain=NonNegativeReals,
            units=pyunits.mole * pyunits.meter**-2 * pyunits.second**-1,
            doc="Flux generated",
        )
        self.membrane_fixed_catalyst_ael = Var(
            initialize=5e3,
            bounds=(1e-1, 1e5),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Catalyst - AEL of the BPM",
        )
        self.membrane_fixed_catalyst_cel = Var(
            initialize=5e3,
            bounds=(1e-1, 1e5),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Catalyst - CEL of the BPM",
        )

        self.k_a = Var(
            initialize=1e-3,
            bounds=(1e-6, 1e5),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Equilibrium constant of proton disassociation",
        )
        self.k_b = Var(
            initialize=3e-2,
            bounds=(1e-2, 1e5),
            units=pyunits.mole * pyunits.meter**-3,
            doc="Equilibrium constant of hydroxide disassociation",
        )

        const = 0.0936 * pyunits.K**2 * pyunits.volt**-1 * pyunits.meter

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the non-dimensional potential drop across the depletion region",
        )
        def eq_voltage_membrane_drop_non_dim(self, t, x):
            return self.elec_field_non_dim[t, x] == const * self.basic.properties[
                t, x
            ].temperature ** -2 * self.relative_permittivity**-1 * sqrt(
                (
                    Constants.faraday_constant
                    * self.membrane_fixed_charge
                    * self.voltage_membrane_drop[t, x]
                )
                / (Constants.vacuum_electric_permittivity * self.relative_permittivity)
            )

        @self.Constraint(
            self.flowsheet().time,
            self.diluate.length_domain,
            doc="Calculate the potential barrier at limiting current across the bipolar membrane",
        )
        def eq_flux_splitting(self, t, x):
            matrx = 0
            for indx in range(self.config.terms_fE):
                matrx += (
                    2**indx
                    * self.elec_field_non_dim[t, x] ** indx
                    / (math.factorial(indx) * math.factorial(indx + 1))
                )

            matrx *= self.k2_zero * self.conc_water
            matrx_a = matrx * self.membrane_fixed_catalyst_cel / self.k_a
            matrx_b = matrx * self.membrane_fixed_catalyst_ael / self.k_b
            return self.flux_splitting[t, x] == (matrx_a + matrx_b) * sqrt(
                self.voltage_membrane_drop[t, x]
                * Constants.vacuum_electric_permittivity
                * self.relative_permittivity
                / (Constants.faraday_constant * self.membrane_fixed_charge)
            )

    def _get_fluid_dimensionless_quantities(self):
        self.hydraulic_diameter = Var(
            initialize=1e-3,
            bounds=(0, None),
            units=pyunits.meter,
            doc="The hydraulic diameter of the channel",
        )
        self.N_Re = Var(
            initialize=50,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Reynolds Number",
        )
        self.N_Sc = Var(
            initialize=2000,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Schmidt Number",
        )
        self.N_Sh = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Sherwood Number",
        )

        if self.config.hydraulic_diameter_method == HydraulicDiameterMethod.fixed:
            _log.warning("Do not forget to FIX the channel hydraulic diameter in [m]!")
        else:

            @self.Constraint(
                doc="To calculate hydraulic diameter",
            )
            def eq_hydraulic_diameter(self):
                if (
                    self.config.hydraulic_diameter_method
                    == HydraulicDiameterMethod.conventional
                ):
                    return (
                        self.hydraulic_diameter
                        == 2
                        * self.channel_height
                        * self.cell_width
                        * self.shadow_factor
                        * self.spacer_porosity
                        * (self.channel_height + self.cell_width * self.shadow_factor)
                        ** -1
                    )
                else:
                    self.spacer_specific_area = Var(
                        initialize=1e4,
                        bounds=(0, None),
                        units=pyunits.meter**-1,
                        doc="The specific area of the channel",
                    )
                    return (
                        self.hydraulic_diameter
                        == 4
                        * self.spacer_porosity
                        * (
                            2 * self.channel_height**-1
                            + (1 - self.spacer_porosity) * self.spacer_specific_area
                        )
                        ** -1
                    )

        @self.Constraint(
            doc="To calculate Re",
        )
        def eq_Re(self):

            return (
                self.N_Re
                == self.dens_mass
                * self.velocity_diluate[0, 0]
                * self.hydraulic_diameter
                * self.visc_d**-1
            )

        @self.Constraint(
            doc="To calculate Sc",
        )
        def eq_Sc(self):

            return self.N_Sc == self.visc_d * self.dens_mass**-1 * self.diffus_mass**-1

        @self.Constraint(
            doc="To calculate Sh",
        )
        def eq_Sh(self):

            return self.N_Sh == 0.29 * self.N_Re**0.5 * self.N_Sc**0.33

    def _pressure_drop_calculation(self):
        self.pressure_drop = Var(
            self.flowsheet().time,
            initialize=1e4,
            units=pyunits.pascal * pyunits.meter**-1,
            doc="pressure drop per unit of length",
        )
        self.pressure_drop_total = Var(
            self.flowsheet().time,
            initialize=1e6,
            units=pyunits.pascal,
            doc="pressure drop over an entire ED stack",
        )

        if self.config.pressure_drop_method == PressureDropMethod.experimental:
            _log.warning(
                "Do not forget to FIX the experimental pressure drop value in [Pa/m]!"
            )
        else:  # PressureDropMethod.Darcy_Weisbach is used

            self._get_fluid_dimensionless_quantities()

            self.friction_factor = Var(
                initialize=10,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="friction factor of the channel fluid",
            )

            @self.Constraint(
                self.flowsheet().time,
                doc="To calculate pressure drop per unit length",
            )
            def eq_pressure_drop(self, t):
                return (
                    self.pressure_drop[t]
                    == self.dens_mass
                    * self.friction_factor
                    * self.velocity_diluate[0, 0] ** 2
                    * 0.5
                    * self.hydraulic_diameter**-1
                )

            if self.config.friction_factor_method == FrictionFactorMethod.fixed:
                _log.warning("Do not forget to FIX the Darcy's friction factor value!")
            else:

                @self.Constraint(
                    doc="To calculate friction factor",
                )
                def eq_friction_factor(self):
                    if (
                        self.config.friction_factor_method
                        == FrictionFactorMethod.Gurreri
                    ):
                        return (
                            self.friction_factor
                            == 4 * 50.6 * self.spacer_porosity**-7.06 * self.N_Re**-1
                        )
                    elif (
                        self.config.friction_factor_method
                        == FrictionFactorMethod.Kuroda
                    ):
                        return (
                            self.friction_factor
                            == 4 * 9.6 * self.spacer_porosity**-1 * self.N_Re**-0.5
                        )

        @self.Constraint(
            self.flowsheet().time,
            doc="To calculate total pressure drop over a stack",
        )
        def eq_pressure_drop_total(self, t):
            return (
                self.pressure_drop_total[t] == self.pressure_drop[t] * self.cell_length
            )

    # initialize method
    def initialize_build(
        self,
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
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # Set the intial conditions over the 1D length from the state vars -dilate
        for k in self.keys():
            for set in self[k].diluate.properties:
                if ("flow_mol_phase_comp" or "flow_mass_phase_comp") not in self[
                    k
                ].diluate.properties[set].define_state_vars():
                    raise ConfigurationError(
                        "Electrodialysis1D unit model requires "
                        "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                        "state variable basis to apply the 'propogate_initial_state' method"
                    )
                if "temperature" in self[k].diluate.properties[set].define_state_vars():
                    self[k].diluate.properties[set].temperature = value(
                        self[k].diluate.properties[(0.0, 0.0)].temperature
                    )
                if "pressure" in self[k].diluate.properties[set].define_state_vars():
                    self[k].diluate.properties[set].pressure = value(
                        self[k].diluate.properties[(0.0, 0.0)].pressure
                    )
                if (
                    "flow_mol_phase_comp"
                    in self[k].diluate.properties[set].define_state_vars()
                ):
                    for ind in self[k].diluate.properties[set].flow_mol_phase_comp:
                        self[k].diluate.properties[set].flow_mol_phase_comp[ind] = (
                            value(
                                self[k]
                                .diluate.properties[(0.0, 0.0)]
                                .flow_mol_phase_comp[ind]
                            )
                        )
                if (
                    "flow_mass_phase_comp"
                    in self[k].diluate.properties[set].define_state_vars()
                ):
                    for ind in self[k].diluate.properties[set].flow_mass_phase_comp:
                        self[k].diluate.properties[set].flow_mass_phase_comp[ind] = (
                            value(
                                self[k]
                                .diluate.properties[(0.0, 0.0)]
                                .flow_mass_phase_comp[ind]
                            )
                        )
                self[k].total_areal_resistance_x[set].set_value(
                    (
                        pyunits.ohm
                        * pyunits.meter**2
                        * (
                            (
                                0.108
                                * pyunits.kg
                                * pyunits.meter**-3
                                / (
                                    self.acidic.properties[set].conc_mass_phase_comp[
                                        "Liq", "H_+"
                                    ]
                                    + self.acidic.properties[set].conc_mass_phase_comp[
                                        "Liq", "Cl_-"
                                    ]
                                    + self.basic.properties[set].conc_mass_phase_comp[
                                        "Liq", "Na_+"
                                    ]
                                    + self.basic.properties[set].conc_mass_phase_comp[
                                        "Liq", "OH_-"
                                    ]
                                )
                                + 0.0492
                            )
                            / 5
                        )
                        + self[k].channel_height
                        * (
                            self[k].acidic.properties[set].elec_cond_phase["Liq"] ** -1
                            + self[k].basic.properties[set].elec_cond_phase["Liq"] ** -1
                            + self[k].diluate.properties[set].elec_cond_phase["Liq"]
                            ** -1
                        )
                    )
                    * self[k].cell_triplet_num
                    + self[k].electrodes_resistance
                )

        # Set the intial conditions over the 1D length from the state vars - basic
        for k in self.keys():
            for set in self[k].basic.properties:
                if ("flow_mol_phase_comp" or "flow_mass_phase_comp") not in self[
                    k
                ].basic.properties[set].define_state_vars():
                    raise ConfigurationError(
                        "Electrodialysis1D unit model requires "
                        "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                        "state variable basis to apply the 'propogate_initial_state' method"
                    )
                if "temperature" in self[k].basic.properties[set].define_state_vars():
                    self[k].basic.properties[set].temperature = value(
                        self[k].basic.properties[(0.0, 0.0)].temperature
                    )
                if "pressure" in self[k].basic.properties[set].define_state_vars():
                    self[k].basic.properties[set].pressure = value(
                        self[k].basic.properties[(0.0, 0.0)].pressure
                    )
                if (
                    "flow_mol_phase_comp"
                    in self[k].basic.properties[set].define_state_vars()
                ):
                    for ind in self[k].basic.properties[set].flow_mol_phase_comp:
                        self[k].basic.properties[set].flow_mol_phase_comp[ind] = value(
                            self[k]
                            .basic.properties[(0.0, 0.0)]
                            .flow_mol_phase_comp[ind]
                        )
                if (
                    "flow_mass_phase_comp"
                    in self[k].basic.properties[set].define_state_vars()
                ):
                    for ind in self[k].basic.properties[set].flow_mass_phase_comp:
                        self[k].basic.properties[set].flow_mass_phase_comp[ind] = value(
                            self[k]
                            .basic.properties[(0.0, 0.0)]
                            .flow_mass_phase_comp[ind]
                        )

                # Set the intial conditions over the 1D length from the state vars - acidic
                for k in self.keys():
                    for set in self[k].acidic.properties:
                        if (
                            "flow_mol_phase_comp" or "flow_mass_phase_comp"
                        ) not in self[k].acidic.properties[set].define_state_vars():
                            raise ConfigurationError(
                                "Electrodialysis1D unit model requires "
                                "either a 'flow_mol_phase_comp' or 'flow_mass_phase_comp' "
                                "state variable basis to apply the 'propogate_initial_state' method"
                            )
                        if (
                            "temperature"
                            in self[k].acidic.properties[set].define_state_vars()
                        ):
                            self[k].acidic.properties[set].temperature = value(
                                self[k].acidic.properties[(0.0, 0.0)].temperature
                            )
                        if (
                            "pressure"
                            in self[k].acidic.properties[set].define_state_vars()
                        ):
                            self[k].acidic.properties[set].pressure = value(
                                self[k].acidic.properties[(0.0, 0.0)].pressure
                            )
                        if (
                            "flow_mol_phase_comp"
                            in self[k].acidic.properties[set].define_state_vars()
                        ):
                            for ind in (
                                self[k].acidic.properties[set].flow_mol_phase_comp
                            ):
                                self[k].acidic.properties[set].flow_mol_phase_comp[
                                    ind
                                ] = value(
                                    self[k]
                                    .acidic.properties[(0.0, 0.0)]
                                    .flow_mol_phase_comp[ind]
                                )
                        if (
                            "flow_mass_phase_comp"
                            in self[k].acidic.properties[set].define_state_vars()
                        ):
                            for ind in (
                                self[k].acidic.properties[set].flow_mass_phase_comp
                            ):
                                self[k].acidic.properties[set].flow_mass_phase_comp[
                                    ind
                                ] = value(
                                    self[k]
                                    .acidic.properties[(0.0, 0.0)]
                                    .flow_mass_phase_comp[ind]
                                )

        # ---------------------------------------------------------------------

        # Initialize diluate block
        flags_diluate = self.diluate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        if not ignore_dof:
            check_dof(self, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Initialize concentrate_basic_side block
        flags_basic = self.basic.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------
        # Initialize concentrate_acidic_side block
        flags_acidic = self.acidic.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,  # inlet var
            hold_state=True,
        )
        init_log.info_high("Initialization Step 3 Complete.")
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))
        check_solve(
            res,
            logger=init_log,
            fail_flag=fail_on_warning,
            checkpoint="Initialization Step 4",
        )
        # ---------------------------------------------------------------------
        # Release state
        self.diluate.release_state(flags_diluate, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        self.basic.release_state(flags_basic, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        self.acidic.release_state(flags_acidic, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

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
        if iscale.get_scaling_factor(self.cell_triplet_num, warning=True) is None:
            iscale.set_scaling_factor(self.cell_triplet_num, 0.1)
        if iscale.get_scaling_factor(self.cell_length, warning=True) is None:
            iscale.set_scaling_factor(self.cell_length, 1e1)
        if iscale.get_scaling_factor(self.cell_width, warning=True) is None:
            iscale.set_scaling_factor(self.cell_width, 1e2)
        if iscale.get_scaling_factor(self.shadow_factor, warning=True) is None:
            iscale.set_scaling_factor(self.shadow_factor, 1e0)
        if iscale.get_scaling_factor(self.channel_height, warning=True) is None:
            iscale.set_scaling_factor(self.channel_height, 1e5)
        if iscale.get_scaling_factor(self.spacer_porosity, warning=True) is None:
            iscale.set_scaling_factor(self.spacer_porosity, 1)
        if iscale.get_scaling_factor(self.electrodes_resistance, warning=True) is None:
            iscale.set_scaling_factor(self.electrodes_resistance, 1e1)
        if hasattr(self, "voltage_applied") and (
            iscale.get_scaling_factor(self.voltage_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.voltage_applied, 1)
        if hasattr(self, "current_applied") and (
            iscale.get_scaling_factor(self.current_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.current_applied, 1)
        if hasattr(self, "conc_water") and (
            iscale.get_scaling_factor(self.conc_water, warning=True) is None
        ):
            iscale.set_scaling_factor(self.conc_water, 1e-4)

        if hasattr(self, "voltage_applied") and (
            iscale.get_scaling_factor(self.voltage_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.voltage_applied, 1)
        if hasattr(self, "current_applied") and (
            iscale.get_scaling_factor(self.current_applied, warning=True) is None
        ):
            iscale.set_scaling_factor(self.current_applied, 1)

        if (
            iscale.get_scaling_factor(self.membrane_fixed_catalyst_cel, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.membrane_fixed_catalyst_cel, 1e-3)
        if (
            iscale.get_scaling_factor(self.membrane_fixed_catalyst_ael, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.membrane_fixed_catalyst_ael, 1e-3)

        if iscale.get_scaling_factor(self.k_a, warning=True) is None:
            iscale.set_scaling_factor(self.k_a, 1e6)

        if iscale.get_scaling_factor(self.k_b, warning=True) is None:
            iscale.set_scaling_factor(self.k_b, 1e2)

        if iscale.get_scaling_factor(self.elec_field_non_dim, warning=True) is None:
            iscale.set_scaling_factor(self.elec_field_non_dim, 1e-1)

        if iscale.get_scaling_factor(self.voltage_membrane_drop, warning=True) is None:
            iscale.set_scaling_factor(self.voltage_membrane_drop, 1)

        if hasattr(self, "membrane_fixed_charge") and (
            iscale.get_scaling_factor(self.membrane_fixed_charge, warning=True) is None
        ):
            iscale.set_scaling_factor(self.membrane_fixed_charge, 1e-3)
        if hasattr(self, "diffus_mass") and (
            iscale.get_scaling_factor(self.diffus_mass, warning=True) is None
        ):
            iscale.set_scaling_factor(self.diffus_mass, 1e9)
        if hasattr(self, "salt_conc_ael_x") and (
            iscale.get_scaling_factor(self.salt_conc_ael_x, warning=True) is None
        ):
            if self.config.salt_calculation:
                sf = (
                    iscale.get_scaling_factor(
                        self.basic.properties[0, 0].conc_mol_phase_comp["Liq", "Na_+"]
                    )
                    ** 2
                    + iscale.get_scaling_factor(
                        self.basic.properties[0, 0].conc_mol_phase_comp["Liq", "Cl_-"]
                    )
                    ** 2
                ) ** 0.5
            else:
                sf = value(self.salt_conc_ael_ref) ** -1
            iscale.set_scaling_factor(self.salt_conc_ael_x, sf)
        if hasattr(self, "salt_conc_cel_x") and (
            iscale.get_scaling_factor(self.salt_conc_cel_x, warning=True) is None
        ):
            if self.config.salt_calculation:
                sf = (
                    iscale.get_scaling_factor(
                        self.acidic.properties[0, 0].conc_mol_phase_comp["Liq", "Na_+"]
                    )
                    ** 2
                    + iscale.get_scaling_factor(
                        self.acidic.properties[0, 0].conc_mol_phase_comp["Liq", "Cl_-"]
                    )
                    ** 2
                ) ** 0.5
            else:
                sf = value(self.salt_conc_cel_ref) ** -1
            iscale.set_scaling_factor(self.salt_conc_cel_x, sf)

        if hasattr(self, "relative_permittivity") and (
            iscale.get_scaling_factor(self.relative_permittivity, warning=True) is None
        ):
            iscale.set_scaling_factor(self.relative_permittivity, 1e-1)

        if (
            hasattr(self, "k2_zero")
            and iscale.get_scaling_factor(self.k2_zero, warning=True) is None
        ):
            iscale.set_scaling_factor(self.k2_zero, 1e5)

        if (
            hasattr(self, "voltage_membrane_drop")
            and iscale.get_scaling_factor(self.voltage_membrane_drop, warning=True)
            is None
        ):
            iscale.set_scaling_factor(self.voltage_membrane_drop, 1e0)

        # The folloing Vars are built for constructing constraints and their sf are computed from other Vars.

        if iscale.get_scaling_factor(self.flux_splitting, warning=True) is None:

            sf = 0
            for indx in range(self.config.terms_fE):
                sf += (
                    2**indx
                    * iscale.get_scaling_factor(self.elec_field_non_dim) ** -indx
                    / (math.factorial(indx) * math.factorial(indx + 1))
                )

            sf **= -1
            sf *= iscale.get_scaling_factor(self.k2_zero) * iscale.get_scaling_factor(
                self.conc_water
            )
            sf_a = (
                sf
                * iscale.get_scaling_factor(self.membrane_fixed_catalyst_cel)
                / iscale.get_scaling_factor(self.k_a)
            )
            sf_b = (
                sf
                * iscale.get_scaling_factor(self.membrane_fixed_catalyst_ael)
                / iscale.get_scaling_factor(self.k_b)
            )

            sf = (sf_a**-1 + sf_b**-1) ** -1 * sqrt(
                iscale.get_scaling_factor(self.voltage_membrane_drop)
                * value(Constants.vacuum_electric_permittivity) ** -1
                * iscale.get_scaling_factor(self.relative_permittivity)
                / (
                    value(Constants.faraday_constant) ** -1
                    * iscale.get_scaling_factor(self.membrane_fixed_charge)
                )
            )

            iscale.set_scaling_factor(self.flux_splitting, sf)

        for ind in self.total_areal_resistance_x:
            if (
                iscale.get_scaling_factor(
                    self.total_areal_resistance_x[ind], warning=False
                )
                is None
            ):
                iscale.set_scaling_factor(self.total_areal_resistance_x[ind], 1e5)

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
                        / iscale.get_scaling_factor(self.shadow_factor)
                        / iscale.get_scaling_factor(self.cell_length)
                    )
                    iscale.set_scaling_factor(self.current_density_x[ind], sf)
                else:
                    sf = iscale.get_scaling_factor(
                        self.voltage_applied
                    ) / iscale.get_scaling_factor(self.total_areal_resistance_x[ind])
                    iscale.set_scaling_factor(self.current_density_x[ind], sf)

        for ind in self.elec_migration_mono_cem_flux:
            iscale.set_scaling_factor(
                self.elec_migration_mono_cem_flux[ind],
                iscale.get_scaling_factor(self.current_density_x[ind[0], ind[1]]) * 1e5,
            )
        for ind in self.elec_migration_mono_aem_flux:
            iscale.set_scaling_factor(
                self.elec_migration_mono_aem_flux[ind],
                iscale.get_scaling_factor(self.current_density_x[ind[0], ind[1]]) * 1e5,
            )
        for ind in self.elec_migration_bpm_flux:
            iscale.set_scaling_factor(
                self.elec_migration_bpm_flux[ind],
                iscale.get_scaling_factor(self.current_density_x[ind[0], ind[1]]) * 1e5,
            )

        for ind in self.generation_cel_flux:
            if ind[3] == "H_+" or "H2O":
                sf = 0.5 * iscale.get_scaling_factor(self.flux_splitting)
            else:
                sf = 1

            iscale.set_scaling_factor(self.generation_cel_flux[ind], sf)

        for ind in self.generation_ael_flux:
            if ind[3] == "OH_-" or "H2O":
                sf = iscale.get_scaling_factor(self.flux_splitting)
            else:
                sf = 1

            iscale.set_scaling_factor(self.generation_ael_flux[ind], sf)

        for ind in self.nonelec_mono_cem_flux:
            if ind[3] == "H2O":
                sf = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * iscale.get_scaling_factor(
                        self.basic.properties[ind[0], ind[1]].pressure_osm_phase[ind[2]]
                    )
                )
            else:
                sf = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * iscale.get_scaling_factor(
                        self.basic.properties[ind[0], ind[1]].conc_mol_phase_comp[
                            ind[2], ind[3]
                        ]
                    )
                )
            iscale.set_scaling_factor(self.nonelec_mono_cem_flux[ind], sf)

        for ind in self.nonelec_mono_aem_flux:
            if ind[3] == "H2O":
                sf = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * iscale.get_scaling_factor(
                        self.acidic.properties[ind[0], ind[1]].pressure_osm_phase[
                            ind[2]
                        ]
                    )
                )
            else:
                sf = (
                    iscale.get_scaling_factor(self.solute_diffusivity_membrane)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * iscale.get_scaling_factor(
                        self.acidic.properties[ind[0], ind[1]].conc_mol_phase_comp[
                            ind[2], ind[3]
                        ]
                    )
                )
            iscale.set_scaling_factor(self.nonelec_mono_aem_flux[ind], sf)

        for ind in self.nonelec_bpm_flux:
            if ind[3] == "H2O":
                sf = (
                    1e-3
                    * 0.018
                    * iscale.get_scaling_factor(self.water_permeability_membrane)
                    * iscale.get_scaling_factor(
                        self.acidic.properties[ind[0], ind[1]].pressure_osm_phase[
                            ind[2]
                        ]
                    )
                )
            else:
                sf = 1
            iscale.set_scaling_factor(self.nonelec_bpm_flux[ind], sf)

        for ind in self.acidic.mass_transfer_term:
            if ind[3] == "H_+":
                sf = iscale.get_scaling_factor(self.generation_cel_flux[ind])
            else:
                if ind[3] == "H2O":
                    sf = iscale.get_scaling_factor(self.nonelec_bpm_flux[ind])
                else:
                    sf = iscale.get_scaling_factor(self.elec_migration_bpm_flux[ind])

            sf *= (
                iscale.get_scaling_factor(self.cell_width)
                * iscale.get_scaling_factor(self.shadow_factor)
                * iscale.get_scaling_factor(self.cell_length)
                * iscale.get_scaling_factor(self.cell_triplet_num)
            )
            iscale.set_scaling_factor(self.acidic.mass_transfer_term[ind], sf)
        #
        for ind in self.basic.mass_transfer_term:
            if ind[3] == "OH_-":
                sf = iscale.get_scaling_factor(self.generation_ael_flux[ind])
            else:
                if ind[3] == "H2O":
                    sf = iscale.get_scaling_factor(self.nonelec_bpm_flux[ind])
                else:
                    sf = iscale.get_scaling_factor(self.elec_migration_bpm_flux[ind])

            sf *= (
                iscale.get_scaling_factor(self.cell_width)
                * iscale.get_scaling_factor(self.shadow_factor)
                * iscale.get_scaling_factor(self.cell_length)
                * iscale.get_scaling_factor(self.cell_triplet_num)
            )
            iscale.set_scaling_factor(self.basic.mass_transfer_term[ind], sf)

        for ind in self.velocity_diluate:
            if (
                iscale.get_scaling_factor(self.velocity_diluate[ind], warning=False)
                is None
            ):
                sf = (
                    iscale.get_scaling_factor(
                        self.diluate.properties[ind].flow_vol_phase["Liq"]
                    )
                    * iscale.get_scaling_factor(self.cell_width) ** -1
                    * iscale.get_scaling_factor(self.shadow_factor) ** -1
                    * iscale.get_scaling_factor(self.channel_height) ** -1
                    * iscale.get_scaling_factor(self.spacer_porosity) ** -1
                    * iscale.get_scaling_factor(self.cell_triplet_num) ** -1
                )

                iscale.set_scaling_factor(self.velocity_diluate[ind], sf)

            for ind in self.velocity_basic:
                if (
                    iscale.get_scaling_factor(self.velocity_basic[ind], warning=False)
                    is None
                ):
                    sf = (
                        iscale.get_scaling_factor(
                            self.basic.properties[ind].flow_vol_phase["Liq"]
                        )
                        * iscale.get_scaling_factor(self.cell_width) ** -1
                        * iscale.get_scaling_factor(self.shadow_factor) ** -1
                        * iscale.get_scaling_factor(self.channel_height) ** -1
                        * iscale.get_scaling_factor(self.spacer_porosity) ** -1
                        * iscale.get_scaling_factor(self.cell_triplet_num) ** -1
                    )

                    iscale.set_scaling_factor(self.velocity_basic[ind], sf)

            for ind in self.velocity_acidic:
                if (
                    iscale.get_scaling_factor(self.velocity_diluate[ind], warning=False)
                    is None
                ):
                    sf = (
                        iscale.get_scaling_factor(
                            self.acidic.properties[ind].flow_vol_phase["Liq"]
                        )
                        * iscale.get_scaling_factor(self.cell_width) ** -1
                        * iscale.get_scaling_factor(self.shadow_factor) ** -1
                        * iscale.get_scaling_factor(self.channel_height) ** -1
                        * iscale.get_scaling_factor(self.spacer_porosity) ** -1
                        * iscale.get_scaling_factor(self.cell_triplet_num) ** -1
                    )

                    iscale.set_scaling_factor(self.velocity_acidic[ind], sf)

        for ind in self.voltage_x:
            if iscale.get_scaling_factor(self.voltage_x[ind], warning=False) is None:
                sf = iscale.get_scaling_factor(
                    self.current_density_x[ind]
                ) * iscale.get_scaling_factor(self.total_areal_resistance_x[ind])
                iscale.set_scaling_factor(self.voltage_x[ind], sf)

        if iscale.get_scaling_factor(self.spacer_porosity, warning=False) is None:
            iscale.set_scaling_factor(self.spacer_porosity, 1)

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
                    * iscale.get_scaling_factor(self.shadow_factor)
                    * iscale.get_scaling_factor(self.cell_length),
                )
        for ind in self.diluate.Dpower_electrical_Dx:
            if (
                iscale.get_scaling_factor(
                    self.diluate.Dpower_electrical_Dx[ind], warning=False
                )
                is None
            ):
                iscale.set_scaling_factor(
                    self.diluate.Dpower_electrical_Dx[ind],
                    iscale.get_scaling_factor(self.diluate.power_electrical_x[ind]),
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
                    * iscale.get_scaling_factor(self.cell_triplet_num)
                )
                ** -1,
            )

        if hasattr(self, "spacer_specific_area") and (
            iscale.get_scaling_factor(self.spacer_specific_area, warning=True) is None
        ):
            iscale.set_scaling_factor(self.spacer_specific_area, 1e-4)
        if hasattr(self, "hydraulic_diameter") and (
            iscale.get_scaling_factor(self.hydraulic_diameter, warning=True) is None
        ):
            iscale.set_scaling_factor(self.hydraulic_diameter, 1e4)
        if hasattr(self, "dens_mass") and (
            iscale.get_scaling_factor(self.dens_mass, warning=True) is None
        ):
            iscale.set_scaling_factor(self.dens_mass, 1e-3)
        if hasattr(self, "N_Re") and (
            iscale.get_scaling_factor(self.N_Re, warning=True) is None
        ):
            sf = (
                iscale.get_scaling_factor(self.dens_mass)
                * iscale.get_scaling_factor(self.velocity_diluate[0, 0])
                * iscale.get_scaling_factor(self.hydraulic_diameter)
                * iscale.get_scaling_factor(self.visc_d) ** -1
            )
            iscale.set_scaling_factor(self.N_Re, sf)
        if hasattr(self, "N_Sc") and (
            iscale.get_scaling_factor(self.N_Sc, warning=True) is None
        ):
            sf = (
                iscale.get_scaling_factor(self.visc_d)
                * iscale.get_scaling_factor(self.dens_mass) ** -1
                * iscale.get_scaling_factor(self.diffus_mass) ** -1
            )
            iscale.set_scaling_factor(self.N_Sc, sf)
        if hasattr(self, "N_Sh") and (
            iscale.get_scaling_factor(self.N_Sh, warning=True) is None
        ):
            sf = (
                10
                * iscale.get_scaling_factor(self.N_Re) ** 0.5
                * iscale.get_scaling_factor(self.N_Sc) ** 0.33
            )
            iscale.set_scaling_factor(self.N_Sh, sf)
        if hasattr(self, "friction_factor") and (
            iscale.get_scaling_factor(self.friction_factor, warning=True) is None
        ):
            if self.config.friction_factor_method == FrictionFactorMethod.fixed:
                sf = 0.1
            elif self.config.friction_factor_method == FrictionFactorMethod.Gurreri:
                sf = (
                    (4 * 50.6) ** -1
                    * (iscale.get_scaling_factor(self.spacer_porosity)) ** -7.06
                    * iscale.get_scaling_factor(self.N_Re) ** -1
                )
            elif self.config.friction_factor_method == FrictionFactorMethod.Kuroda:
                sf = (4 * 9.6) ** -1 * iscale.get_scaling_factor(self.N_Re) ** -0.5
            iscale.set_scaling_factor(self.friction_factor, sf)

        if hasattr(self, "pressure_drop") and (
            iscale.get_scaling_factor(self.pressure_drop, warning=True) is None
        ):
            if self.config.pressure_drop_method == PressureDropMethod.experimental:
                sf = 1e-5
            else:
                sf = (
                    iscale.get_scaling_factor(self.dens_mass)
                    * iscale.get_scaling_factor(self.friction_factor)
                    * iscale.get_scaling_factor(self.velocity_diluate[0, 0]) ** 2
                    * 2
                    * iscale.get_scaling_factor(self.hydraulic_diameter) ** -1
                )
            iscale.set_scaling_factor(self.pressure_drop, sf)

        if hasattr(self, "pressure_drop_total") and (
            iscale.get_scaling_factor(self.pressure_drop_total, warning=True) is None
        ):
            if self.config.pressure_drop_method == PressureDropMethod.experimental:
                sf = 1e-5 * iscale.get_scaling_factor(self.cell_length)
            else:
                sf = (
                    iscale.get_scaling_factor(self.dens_mass)
                    * iscale.get_scaling_factor(self.friction_factor)
                    * iscale.get_scaling_factor(self.velocity_diluate[0, 0]) ** 2
                    * 2
                    * iscale.get_scaling_factor(self.hydraulic_diameter) ** -1
                    * iscale.get_scaling_factor(self.cell_length)
                )
            iscale.set_scaling_factor(self.pressure_drop_total, sf)

        if hasattr(self, "current_dens_lim_bpm"):
            if iscale.get_scaling_factor(self.current_dens_lim_bpm) is None:
                if (
                    self.config.limiting_current_density_method_bpm
                    == LimitingCurrentDensitybpmMethod.InitialValue
                ):
                    sf = self.config.limiting_current_density_data**-1
                    iscale.set_scaling_factor(self.current_dens_lim_bpm, sf)
                elif (
                    self.config.limiting_current_density_method_bpm
                    == LimitingCurrentDensitybpmMethod.Empirical
                ):
                    sf = (
                        iscale.get_scaling_factor(self.diffus_mass)
                        * value(Constants.faraday_constant) ** -1
                        * (2 * iscale.get_scaling_factor(self.salt_conc_ael_x)) ** 2
                        / (
                            iscale.get_scaling_factor(self.membrane_thickness)
                            * iscale.get_scaling_factor(self.membrane_fixed_charge)
                        )
                    )

                    iscale.set_scaling_factor(self.current_dens_lim_bpm, sf)

        for ind, c in self.eq_get_current_density.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.current_density_x[ind])
            )

        for ind, c in self.eq_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.diluate.power_electrical_x[ind])
            )
        for ind, c in self.eq_specific_power_electrical.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.specific_power_electrical[ind])
            )

        for ind, c in self.eq_elec_migration_mono_cem.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.elec_migration_mono_cem_flux[ind])
            )
        for ind, c in self.eq_elec_migration_mono_aem_flux.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.elec_migration_mono_aem_flux[ind])
            )
        for ind, c in self.eq_elec_migration_bpm_flux.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.elec_migration_bpm_flux[ind])
            )

        for ind, c in self.eq_nonelec_mono_cem_flux.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.nonelec_mono_cem_flux[ind])
            )

        for ind, c in self.eq_nonelec_mono_aem_flux.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.nonelec_mono_aem_flux[ind])
            )

        for ind, c in self.eq_nonelec_bpm_flux.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(self.nonelec_bpm_flux[ind])
            )

        for ind, c in self.eq_mass_transfer_term_diluate.items():
            iscale.constraint_scaling_transform(
                c,
                min(
                    iscale.get_scaling_factor(self.elec_migration_mono_cem_flux[ind]),
                    iscale.get_scaling_factor(self.nonelec_mono_cem_flux[ind]),
                ),
            )

        for ind, c in self.eq_mass_transfer_term_basic.items():
            iscale.constraint_scaling_transform(
                c,
                min(
                    iscale.get_scaling_factor(self.generation_ael_flux[ind]),
                    iscale.get_scaling_factor(self.elec_migration_bpm_flux[ind]),
                    iscale.get_scaling_factor(self.nonelec_bpm_flux[ind]),
                )
                * iscale.get_scaling_factor(self.cell_width)
                * iscale.get_scaling_factor(self.shadow_factor)
                * iscale.get_scaling_factor(self.cell_length)
                * iscale.get_scaling_factor(self.cell_triplet_num),
            )
        for ind, c in self.eq_mass_transfer_term_acidic.items():
            iscale.constraint_scaling_transform(
                c,
                min(
                    iscale.get_scaling_factor(self.generation_cel_flux[ind]),
                    iscale.get_scaling_factor(self.nonelec_bpm_flux[ind]),
                )
                * iscale.get_scaling_factor(self.cell_width)
                * iscale.get_scaling_factor(self.shadow_factor)
                * iscale.get_scaling_factor(self.cell_length)
                * iscale.get_scaling_factor(self.cell_triplet_num),
            )

        if hasattr(self, "eq_flux_splitting"):
            for ind, c in self.eq_flux_splitting.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.flux_splitting[ind]),
                )

        if hasattr(self, "eq_salt_cem"):
            for ind, c in self.eq_salt_cem.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.salt_conc_cel_x[ind]),
                )

        if hasattr(self, "eq_salt_aem"):
            for ind, c in self.eq_salt_aem.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.salt_conc_ael_x[ind]),
                )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Diluate Channel Inlet": self.inlet_diluate,
                "base channel of the bipolar membrane Channel Inlet": self.inlet_basic,
                "acid channel of the bipolar membrane Channel Inlet": self.inlet_acidic,
                "Diluate Channel Outlet": self.outlet_diluate,
                "base channel of the bipolar membrane Channel Outlet": self.outlet_basic,
                "acid channel of the bipolar membrane Channel Outlet": self.outlet_acidic,
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
            },
            "exprs": {},
            "params": {},
        }

    def get_power_electrical(self, time_point=0):
        return self.diluate.power_electrical_x[
            time_point, self.diluate.length_domain.last()
        ]

    @property
    def default_costing_method(self):
        return cost_electrodialysis
