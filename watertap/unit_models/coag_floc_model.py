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
    PositiveIntegers,
    Reference,
    value,
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
)
from idaes.core.util.constants import Constants
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Austin Ladshaw"

_log = idaeslog.getLogger(__name__)

# Name of the unit model
@declare_process_block_class("CoagulationFlocculation")
class CoagulationFlocculationData(UnitModelBlockData):
    """
    Zero order Coagulation-Flocculation model based on Jar Tests
    """

    # CONFIG are options for the unit model
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag; must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default**: False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag; must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default**: False. The filtration unit does not have defined volume, thus
    this must be False.""",
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

    # NOTE: This option is temporarily disabled
    '''
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

    CONFIG.declare(
        "chemical_additives",
        ConfigValue(
            default={},
            domain=dict,
            description="Specification of chemical additives used in coagulation process",
            doc="""
        A dict of chemical additives used in coagulation process
        along with their molecular weights, the moles of salt produced per mole of
        chemical added, and the molecular weights of the salt produced by the chemical
        additive with the format of::

            {'chem_name_1':
                {'parameter_data':
                    {
                    'mw_additive': (value, units),
                    'moles_salt_per_mole_additive': value,
                    'mw_salt': (value, units)
                    }
                },
            'chem_name_2':
                {'parameter_data':
                    {
                    'mw_additive': (value, units),
                    'moles_salt_per_mole_additive': value,
                    'mw_salt': (value, units)
                    }
                },
            }

        """,
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

        # check the optional config arg 'chemical_additives'
        common_msg = (
            "The 'chemical_additives' dict MUST contain a dict of 'parameter_data' for "
            + "each chemical name. That 'parameter_data' dict MUST contain 'mw_chem', "
            + "'moles_salt_per_mole_additive', and 'mw_salt' as keys. Users are also "
            + "required to provide the values for the molecular weights and the units "
            + "within a tuple arg. Example format provided below.\n\n"
            + "{'chem_name_1': \n"
            + "     {'parameter_data': \n"
            + "        {'mw_additive': (value, units), \n"
            + "         'moles_salt_per_mole_additive': value, \n"
            + "         'mw_salt': (value, units)} \n"
            + "     }, \n"
            + "}\n\n"
        )
        mw_adds = {}
        mw_salts = {}
        molar_rat = {}
        for j in self.config.chemical_additives:
            if type(self.config.chemical_additives[j]) != dict:
                raise ConfigurationError(
                    "\n Did not provide a 'dict' for chemical \n" + common_msg
                )
            if "parameter_data" not in self.config.chemical_additives[j]:
                raise ConfigurationError(
                    "\n Did not provide a 'parameter_data' for chemical \n" + common_msg
                )
            if "mw_additive" not in self.config.chemical_additives[j]["parameter_data"]:
                raise ConfigurationError(
                    "\n Did not provide a 'mw_additive' for chemical \n" + common_msg
                )
            if (
                "moles_salt_per_mole_additive"
                not in self.config.chemical_additives[j]["parameter_data"]
            ):
                raise ConfigurationError(
                    "\n Did not provide a 'moles_salt_per_mole_additive' for chemical \n"
                    + common_msg
                )
            if "mw_salt" not in self.config.chemical_additives[j]["parameter_data"]:
                raise ConfigurationError(
                    "\n Did not provide a 'mw_salt' for chemical \n" + common_msg
                )
            if (
                type(self.config.chemical_additives[j]["parameter_data"]["mw_additive"])
                != tuple
            ):
                raise ConfigurationError(
                    "\n Did not provide a tuple for 'mw_additive' \n" + common_msg
                )
            if (
                type(self.config.chemical_additives[j]["parameter_data"]["mw_salt"])
                != tuple
            ):
                raise ConfigurationError(
                    "\n Did not provide a tuple for 'mw_salt' \n" + common_msg
                )
            if not isinstance(
                self.config.chemical_additives[j]["parameter_data"][
                    "moles_salt_per_mole_additive"
                ],
                (int, float),
            ):
                raise ConfigurationError(
                    "\n Did not provide a number for 'moles_salt_per_mole_additive' \n"
                    + common_msg
                )

            # Populate temp dicts for parameter and variable setting
            mw_adds[j] = pyunits.convert_value(
                self.config.chemical_additives[j]["parameter_data"]["mw_additive"][0],
                from_units=self.config.chemical_additives[j]["parameter_data"][
                    "mw_additive"
                ][1],
                to_units=pyunits.kg / pyunits.mol,
            )
            mw_salts[j] = pyunits.convert_value(
                self.config.chemical_additives[j]["parameter_data"]["mw_salt"][0],
                from_units=self.config.chemical_additives[j]["parameter_data"][
                    "mw_salt"
                ][1],
                to_units=pyunits.kg / pyunits.mol,
            )
            molar_rat[j] = self.config.chemical_additives[j]["parameter_data"][
                "moles_salt_per_mole_additive"
            ]

        # Add unit variables
        # Linear relationship between TSS (mg/L) and Turbidity (NTU)
        #           TSS (mg/L) = Turbidity (NTU) * slope + intercept
        #   Default values come from the following paper:
        #       H. Rugner, M. Schwientek,B. Beckingham, B. Kuch, P. Grathwohl,
        #       Environ. Earth Sci. 69 (2013) 373-380. DOI: 10.1007/s12665-013-2307-1
        self.slope = Var(
            self.flowsheet().config.time,
            initialize=1.86,
            bounds=(0.0, 10),
            domain=NonNegativeReals,
            units=pyunits.mg / pyunits.L,
            doc="Slope relation between TSS (mg/L) and Turbidity (NTU)",
        )

        self.intercept = Var(
            self.flowsheet().config.time,
            initialize=0,
            bounds=(0, 10),
            domain=NonNegativeReals,
            units=pyunits.mg / pyunits.L,
            doc="Intercept relation between TSS (mg/L) and Turbidity (NTU)",
        )

        self.initial_turbidity_ntu = Var(
            self.flowsheet().config.time,
            initialize=50,
            bounds=(0, 10000),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Initial measured Turbidity (NTU) from Jar Test",
        )

        self.final_turbidity_ntu = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 10000),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Final measured Turbidity (NTU) from Jar Test",
        )

        self.chemical_doses = Var(
            self.flowsheet().config.time,
            self.config.chemical_additives.keys(),
            initialize=0,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.mg / pyunits.L,
            doc="Dosages of the set of chemical additives",
        )

        self.chemical_mw = Param(
            self.config.chemical_additives.keys(),
            mutable=True,
            initialize=mw_adds,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weights of the set of chemical additives",
        )

        self.salt_mw = Param(
            self.config.chemical_additives.keys(),
            mutable=True,
            initialize=mw_salts,
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weights of the produced salts from chemical additives",
        )

        self.salt_from_additive_mole_ratio = Param(
            self.config.chemical_additives.keys(),
            mutable=True,
            initialize=molar_rat,
            domain=NonNegativeReals,
            units=pyunits.mol / pyunits.mol,
            doc="Moles of the produced salts from 1 mole of chemical additives",
        )

        # Below set of Vars are for the power usage of the unit
        #       User's will need to provide scaling factors for these
        # -----------------------------------------------------------
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.rapid_mixing_retention_time = Var(
            self.flowsheet().config.time,
            initialize=30,
            bounds=(0.1, 10000),
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Hydraulic retention time of each rapid mixing basin in seconds",
        )

        self.num_rapid_mixing_basins = Var(
            initialize=1,
            bounds=(1, 10),
            domain=PositiveIntegers,
            units=pyunits.dimensionless,
            doc="Number of rapid mixing basins in series",
        )

        self.rapid_mixing_vel_grad = Var(
            self.flowsheet().config.time,
            initialize=250,
            bounds=(0.1, 10000),
            domain=NonNegativeReals,
            units=pyunits.s**-1,
            doc="Velocity gradient in each rapid mixing basin in (m/s)/m",
        )

        # NOTE: There are 2 modes for flocculation mixing discussed in literature
        #       Here we are only intially defining the 'Paddle-Wheel' mode. Other
        #       modes can be added later (if needed). The 'Paddle-Wheel' configuration
        #       is the most common used for conventional water treatment.
        self.floc_retention_time = Var(
            self.flowsheet().config.time,
            initialize=1800,
            bounds=(10, 10000),
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Hydraulic retention time of the flocculation mixing basin in seconds",
        )

        self.single_paddle_length = Var(
            initialize=2,
            bounds=(0.1, 100),
            domain=NonNegativeReals,
            units=pyunits.m,
            doc="Length of a single paddle blade (from center of rotation to the edge) in meters",
        )

        self.single_paddle_width = Var(
            initialize=0.5,
            bounds=(0.01, 100),
            domain=NonNegativeReals,
            units=pyunits.m,
            doc="Width of a single paddle blade in meters",
        )

        self.paddle_rotational_speed = Var(
            self.flowsheet().config.time,
            initialize=100,
            bounds=(0.01, 10000),
            domain=NonNegativeReals,
            units=pyunits.s**-1,
            doc="Rotational speed of the paddles in revolutions per second",
        )

        self.paddle_drag_coef = Var(
            self.flowsheet().config.time,
            initialize=1.5,
            bounds=(0.1, 10),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Drag coefficient for the paddles in the flocculation basin",
        )

        self.vel_fraction = Var(
            initialize=0.7,
            bounds=(0.6, 0.9),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Fraction of actual paddle velocity relative to local water velocity",
        )

        self.num_paddle_wheels = Var(
            initialize=4,
            bounds=(1, 10),
            domain=PositiveIntegers,
            units=pyunits.dimensionless,
            doc="Number of rotating paddle wheels in the flocculation basin",
        )

        self.num_paddles_per_wheel = Var(
            initialize=4,
            bounds=(1, 10),
            domain=PositiveIntegers,
            units=pyunits.dimensionless,
            doc="Number of paddles attached to each rotating wheel in the flocculation basin",
        )

        # Build control volume for feed side
        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        # NOTE: This checks for if an energy_balance_type is defined
        if hasattr(self.config, "energy_balance_type"):
            self.control_volume.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=False
        )

        # Add ports
        self.add_inlet_port(name="inlet", block=self.control_volume)
        self.add_outlet_port(name="outlet", block=self.control_volume)

        # Check _phase_component_set for required items
        if ("Liq", "TDS") not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','TDS') as a component, but "
                "the property package has only specified the following components {}".format(
                    [p for p in self.config.property_package._phase_component_set]
                )
            )
        if ("Liq", "Sludge") not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','Sludge') as a component, but "
                "the property package has only specified the following components {}".format(
                    [p for p in self.config.property_package._phase_component_set]
                )
            )
        if ("Liq", "TSS") not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','TSS') as a component, but "
                "the property package has only specified the following components {}".format(
                    [p for p in self.config.property_package._phase_component_set]
                )
            )

        # -------- Add constraints ---------
        # Adds isothermal constraint if no energy balance present
        if not hasattr(self.config, "energy_balance_type"):

            @self.Constraint(self.flowsheet().config.time, doc="Isothermal condition")
            def eq_isothermal(self, t):
                return (
                    self.control_volume.properties_out[t].temperature
                    == self.control_volume.properties_in[t].temperature
                )

        # Constraint for tss loss rate based on measured final turbidity
        self.tss_loss_rate = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="Mass per time loss rate of TSS based on the measured final turbidity",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the loss rate of TSS to be used in mass_transfer_term",
        )
        def eq_tss_loss_rate(self, t):
            tss_out = pyunits.convert(
                self.slope[t] * self.final_turbidity_ntu[t] + self.intercept[t],
                to_units=units_meta("mass") * units_meta("length") ** -3,
            )
            input_rate = self.control_volume.properties_in[t].flow_mass_phase_comp[
                "Liq", "TSS"
            ]
            exit_rate = (
                self.control_volume.properties_out[t].flow_vol_phase["Liq"] * tss_out
            )

            return self.tss_loss_rate[t] == input_rate - exit_rate

        # Constraint for tds gain rate based on 'chemical_doses' and 'chemical_additives'
        if self.config.chemical_additives:
            self.tds_gain_rate = Var(
                self.flowsheet().config.time,
                initialize=0,
                bounds=(0, 100),
                domain=NonNegativeReals,
                units=units_meta("mass") * units_meta("time") ** -1,
                doc="Mass per time gain rate of TDS based on the chemicals added for coagulation",
            )

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Constraint for the loss rate of TSS to be used in mass_transfer_term",
            )
            def eq_tds_gain_rate(self, t):
                sum = 0
                for j in self.config.chemical_additives.keys():
                    chem_dose = pyunits.convert(
                        self.chemical_doses[t, j],
                        to_units=units_meta("mass") * units_meta("length") ** -3,
                    )
                    chem_dose = (
                        chem_dose
                        / self.chemical_mw[j]
                        * self.salt_from_additive_mole_ratio[j]
                        * self.salt_mw[j]
                        * self.control_volume.properties_out[t].flow_vol_phase["Liq"]
                    )
                    sum = sum + chem_dose

                return self.tds_gain_rate[t] == sum

        # Add constraints for mass transfer terms
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(self, t, p, j):
            if (p, j) == ("Liq", "TSS"):
                return (
                    self.control_volume.mass_transfer_term[t, p, j]
                    == -self.tss_loss_rate[t]
                )
            elif (p, j) == ("Liq", "Sludge"):
                return (
                    self.control_volume.mass_transfer_term[t, p, j]
                    == self.tss_loss_rate[t]
                )
            elif (p, j) == ("Liq", "TDS"):
                if self.config.chemical_additives:
                    return (
                        self.control_volume.mass_transfer_term[t, p, j]
                        == self.tds_gain_rate[t]
                    )
                else:
                    return self.control_volume.mass_transfer_term[t, p, j] == 0.0
            else:
                return self.control_volume.mass_transfer_term[t, p, j] == 0.0

        # Constraint for the volume of each rapid mixing basin in series
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.rapid_mixing_basin_vol = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=units_meta("length") ** 3,
            doc="Volume of each rapid mixing basin in the series",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the volume of each rapid mixing basin",
        )
        def eq_rapid_mixing_basin_vol(self, t):
            flow_rate = pyunits.convert(
                self.control_volume.properties_in[t].flow_vol_phase["Liq"],
                to_units=units_meta("length") ** 3 / pyunits.s,
            )

            return (
                self.rapid_mixing_basin_vol[t]
                == flow_rate * self.rapid_mixing_retention_time[t]
            )

        # Constraint for the power usage of the rapid mixers
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.rapid_mixing_power = Var(
            self.flowsheet().config.time,
            initialize=0.01,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.kW,
            doc="Power usage of the rapid mixing basins in kW",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the power usage of the rapid mixing basins",
        )
        def eq_rapid_mixing_power(self, t):
            vel_grad = pyunits.convert(
                self.rapid_mixing_vel_grad[t], to_units=units_meta("time") ** -1
            )
            power_usage = pyunits.convert(
                vel_grad**2
                * self.control_volume.properties_out[t].visc_d_phase["Liq"]
                * self.rapid_mixing_basin_vol[t]
                * self.num_rapid_mixing_basins,
                to_units=pyunits.kW,
            )

            return self.rapid_mixing_power[t] == power_usage

        # Constraint for the volume of the flocculation basin
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.floc_basin_vol = Var(
            self.flowsheet().config.time,
            initialize=10,
            bounds=(0, 10000),
            domain=NonNegativeReals,
            units=units_meta("length") ** 3,
            doc="Volume of the flocculation basin",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the volume of the flocculation basin",
        )
        def eq_floc_basin_vol(self, t):
            flow_rate = pyunits.convert(
                self.control_volume.properties_in[t].flow_vol_phase["Liq"],
                to_units=units_meta("length") ** 3 / pyunits.s,
            )

            return self.floc_basin_vol[t] == flow_rate * self.floc_retention_time[t]

        # Constraint for the velocity of the paddle wheels
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.floc_wheel_speed = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Velocity of the wheels in the flocculation basin",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the velocity of the wheels in the flocculation basin",
        )
        def eq_floc_wheel_speed(self, t):
            wheel_rate = pyunits.convert(
                Constants.pi
                * self.single_paddle_length
                * self.paddle_rotational_speed[t],
                to_units=units_meta("length") * units_meta("time") ** -1,
            )

            return self.floc_wheel_speed[t] == wheel_rate

        # Constraint for the power usage of the flocculation mixer
        #   Mines, R.O., Environmental engineering: Principles
        #       and Practice, 1st Ed, John Wiley & Sons, 2014.
        #       Ch. 6.
        self.flocculation_power = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.kW,
            doc="Power usage of the flocculation basin in kW",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the power usage of the flocculation basin",
        )
        def eq_flocculation_power(self, t):
            total_area = pyunits.convert(
                self.single_paddle_width
                * self.single_paddle_length
                * self.num_paddle_wheels
                * self.num_paddles_per_wheel,
                to_units=units_meta("length") ** 2,
            )
            power_usage = pyunits.convert(
                0.5
                * self.paddle_drag_coef[t]
                * total_area
                * self.control_volume.properties_out[t].dens_mass_phase["Liq"]
                * self.vel_fraction**3
                * self.floc_wheel_speed[t] ** 3,
                to_units=pyunits.kW,
            )

            return self.flocculation_power[t] == power_usage

        self.total_power = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.kW,
            doc="Power usage of the full unit model in kW",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint for the power usage of the full unit model",
        )
        def eq_total_power(self, t):
            return (
                self.total_power[t]
                == self.flocculation_power[t] + self.rapid_mixing_power[t]
            )

    # Return a scalar expression for the inlet concentration of TSS
    def compute_inlet_tss_mass_concentration(self, t):
        """
        Function to generate an expression that would represent the mass
        concentration of TSS at the inlet port of the unit. Inlet ports
        are generally established upstream, but this will be useful for
        establishing the inlet TSS when an upstream TSS is unknown. This
        level of inlet TSS is based off of measurements made of Turbidity
        during the Jar Test.

        Keyword Arguments:
            self : this unit model object
            t : time index on the flowsheet

        Returns: Expression

        Recover the numeric value by using 'value(Expression)'
        """
        units_meta = self.config.property_package.get_metadata().get_derived_units
        return pyunits.convert(
            self.slope[t] * self.initial_turbidity_ntu[t] + self.intercept[t],
            to_units=units_meta("mass") * units_meta("length") ** -3,
        )

    # Return a scale expression for the inlet mass flow rate of TSS
    def compute_inlet_tss_mass_flow(self, t):
        """
        Function to generate an expression that would represent the mass
        flow rate of TSS at the inlet port of the unit. Inlet ports
        are generally established upstream, but this will be useful for
        establishing the inlet TSS when an upstream TSS is unknown. This
        level of inlet TSS is based off of measurements made of Turbidity
        during the Jar Test.

        Keyword Arguments:
            self : this unit model object
            t : time index on the flowsheet

        Returns: Expression

        Recover the numeric value by using 'value(Expression)'
        """
        return self.control_volume.properties_in[t].flow_vol_phase[
            "Liq"
        ] * self.compute_inlet_tss_mass_concentration(t)

    # Function to automate fixing of the Turbidity v TSS relation params to defaults
    def fix_tss_turbidity_relation_defaults(self):
        self.slope.fix()
        self.intercept.fix()

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
        # Initialize holdup block
        flags = blk.control_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # scaling factors for turbidity relationship
        #       Supressing warning (these factors are not very important)
        if iscale.get_scaling_factor(self.slope) is None:
            sf = iscale.get_scaling_factor(self.slope, default=1, warning=False)
            iscale.set_scaling_factor(self.slope, sf)
        if iscale.get_scaling_factor(self.intercept) is None:
            sf = iscale.get_scaling_factor(self.intercept, default=1, warning=False)
            iscale.set_scaling_factor(self.intercept, sf)

        # scaling factors for turbidity measurements and chemical doses
        #       Supressing warning
        if iscale.get_scaling_factor(self.initial_turbidity_ntu) is None:
            sf = iscale.get_scaling_factor(
                self.initial_turbidity_ntu, default=1, warning=False
            )
            iscale.set_scaling_factor(self.initial_turbidity_ntu, sf)
        if iscale.get_scaling_factor(self.final_turbidity_ntu) is None:
            sf = iscale.get_scaling_factor(
                self.final_turbidity_ntu, default=1, warning=False
            )
            iscale.set_scaling_factor(self.final_turbidity_ntu, sf)
        if iscale.get_scaling_factor(self.chemical_doses) is None:
            sf = iscale.get_scaling_factor(
                self.chemical_doses, default=1, warning=False
            )
            iscale.set_scaling_factor(self.chemical_doses, sf)

        # scaling factors for power usage in rapid mixing
        #       Display warning
        if iscale.get_scaling_factor(self.rapid_mixing_retention_time) is None:
            sf = iscale.get_scaling_factor(
                self.rapid_mixing_retention_time, default=1e-1, warning=True
            )
            iscale.set_scaling_factor(self.rapid_mixing_retention_time, sf)
        if iscale.get_scaling_factor(self.num_rapid_mixing_basins) is None:
            sf = iscale.get_scaling_factor(
                self.num_rapid_mixing_basins, default=1, warning=True
            )
            iscale.set_scaling_factor(self.num_rapid_mixing_basins, sf)
        if iscale.get_scaling_factor(self.rapid_mixing_vel_grad) is None:
            sf = iscale.get_scaling_factor(
                self.rapid_mixing_vel_grad, default=1e-2, warning=True
            )
            iscale.set_scaling_factor(self.rapid_mixing_vel_grad, sf)
        if iscale.get_scaling_factor(self.floc_retention_time) is None:
            sf = iscale.get_scaling_factor(
                self.floc_retention_time, default=1e-3, warning=True
            )
            iscale.set_scaling_factor(self.floc_retention_time, sf)
        if iscale.get_scaling_factor(self.single_paddle_length) is None:
            sf = iscale.get_scaling_factor(
                self.single_paddle_length, default=1, warning=True
            )
            iscale.set_scaling_factor(self.single_paddle_length, sf)
        if iscale.get_scaling_factor(self.single_paddle_width) is None:
            sf = iscale.get_scaling_factor(
                self.single_paddle_width, default=1, warning=True
            )
            iscale.set_scaling_factor(self.single_paddle_width, sf)
        if iscale.get_scaling_factor(self.paddle_rotational_speed) is None:
            sf = iscale.get_scaling_factor(
                self.paddle_rotational_speed, default=10, warning=True
            )
            iscale.set_scaling_factor(self.paddle_rotational_speed, sf)
        if iscale.get_scaling_factor(self.paddle_drag_coef) is None:
            sf = iscale.get_scaling_factor(
                self.paddle_drag_coef, default=1, warning=True
            )
            iscale.set_scaling_factor(self.paddle_drag_coef, sf)
        if iscale.get_scaling_factor(self.vel_fraction) is None:
            sf = iscale.get_scaling_factor(self.vel_fraction, default=1, warning=True)
            iscale.set_scaling_factor(self.vel_fraction, sf)
        if iscale.get_scaling_factor(self.num_paddle_wheels) is None:
            sf = iscale.get_scaling_factor(
                self.num_paddle_wheels, default=1, warning=True
            )
            iscale.set_scaling_factor(self.num_paddle_wheels, sf)
        if iscale.get_scaling_factor(self.num_paddles_per_wheel) is None:
            sf = iscale.get_scaling_factor(
                self.num_paddles_per_wheel, default=1, warning=True
            )
            iscale.set_scaling_factor(self.num_paddles_per_wheel, sf)

        # set scaling for tss_loss_rate
        if iscale.get_scaling_factor(self.tss_loss_rate) is None:
            sf = 0
            for t in self.control_volume.properties_in:
                sf += value(
                    self.control_volume.properties_in[t].flow_mass_phase_comp[
                        "Liq", "TSS"
                    ]
                )
            sf = sf / len(self.control_volume.properties_in)
            if sf < 0.01:
                sf = 0.01
            iscale.set_scaling_factor(self.tss_loss_rate, 1 / sf)

            for ind, c in self.eq_tss_loss_rate.items():
                iscale.constraint_scaling_transform(c, 1 / sf)

        # set scaling for tds_gain_rate
        if self.config.chemical_additives:
            if iscale.get_scaling_factor(self.tds_gain_rate) is None:
                sf = 0
                for t in self.control_volume.properties_in:
                    sum = 0
                    for j in self.config.chemical_additives.keys():
                        chem_dose = pyunits.convert(
                            self.chemical_doses[t, j],
                            to_units=units_meta("mass") * units_meta("length") ** -3,
                        )
                        chem_dose = (
                            chem_dose
                            / self.chemical_mw[j]
                            * self.salt_from_additive_mole_ratio[j]
                            * self.salt_mw[j]
                            * self.control_volume.properties_in[t].flow_vol_phase["Liq"]
                        )
                        sum = sum + chem_dose
                    sf += value(sum)
                sf = sf / len(self.control_volume.properties_in)
                if sf < 0.001:
                    sf = 0.001
                iscale.set_scaling_factor(self.tds_gain_rate, 1 / sf)

                for ind, c in self.eq_tds_gain_rate.items():
                    iscale.constraint_scaling_transform(c, 1 / sf)

        # set scaling for mass transfer terms
        for ind, c in self.eq_mass_transfer_term.items():
            if ind[2] == "TDS":
                if self.config.chemical_additives:
                    sf = iscale.get_scaling_factor(self.tds_gain_rate)
                else:
                    sf = 1
            elif ind[2] == "TSS":
                sf = iscale.get_scaling_factor(self.tss_loss_rate)
            elif ind[2] == "Sludge":
                sf = iscale.get_scaling_factor(self.tss_loss_rate)
            else:
                sf = 1
            iscale.constraint_scaling_transform(c, sf)
            iscale.set_scaling_factor(self.control_volume.mass_transfer_term[ind], sf)

        # set scaling for rapid_mixing_basin_vol
        if iscale.get_scaling_factor(self.rapid_mixing_basin_vol[t]) is None:
            sf1 = 0
            sf2 = 0
            for t in self.control_volume.properties_out:
                sf1 += iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].flow_vol_phase["Liq"]
                )
            sf2 = iscale.get_scaling_factor(self.rapid_mixing_retention_time)
            sf1 = sf1 / len(self.control_volume.properties_in)
            sf = sf1 * sf2
            iscale.set_scaling_factor(self.rapid_mixing_basin_vol, sf)

            for ind, c in self.eq_rapid_mixing_basin_vol.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling for rapid_mixing_power
        if iscale.get_scaling_factor(self.rapid_mixing_power[t]) is None:
            sf1 = 0
            sf2 = 0
            for t in self.control_volume.properties_out:
                sf1 += iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].visc_d_phase["Liq"]
                )
            sf2 = iscale.get_scaling_factor(self.rapid_mixing_vel_grad)
            sf3 = iscale.get_scaling_factor(self.rapid_mixing_basin_vol)
            sf4 = iscale.get_scaling_factor(self.num_rapid_mixing_basins)
            sf1 = sf1 / len(self.control_volume.properties_in)
            sf = sf1 * sf2**2 * sf3 * sf4
            iscale.set_scaling_factor(self.rapid_mixing_power, sf)

            for ind, c in self.eq_rapid_mixing_power.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling for floc_basin_vol
        if iscale.get_scaling_factor(self.floc_basin_vol[t]) is None:
            sf1 = 0
            sf2 = 0
            for t in self.control_volume.properties_out:
                sf1 += iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].flow_vol_phase["Liq"]
                )
            sf2 = iscale.get_scaling_factor(self.floc_retention_time)
            sf1 = sf1 / len(self.control_volume.properties_in)
            sf = sf1 * sf2
            iscale.set_scaling_factor(self.floc_basin_vol, sf)

            for ind, c in self.eq_floc_basin_vol.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling for floc_wheel_speed
        if iscale.get_scaling_factor(self.floc_wheel_speed[t]) is None:
            sf1 = iscale.get_scaling_factor(self.paddle_rotational_speed)
            sf2 = iscale.get_scaling_factor(self.single_paddle_length)
            sf = sf1 * sf2 * Constants.pi / 10
            iscale.set_scaling_factor(self.floc_wheel_speed, sf)

            for ind, c in self.eq_floc_wheel_speed.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling for flocculation_power
        if iscale.get_scaling_factor(self.flocculation_power[t]) is None:
            sf1 = iscale.get_scaling_factor(self.floc_wheel_speed)
            sf2 = iscale.get_scaling_factor(self.vel_fraction)
            sf3 = 0
            for t in self.control_volume.properties_out:
                sf3 += iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].dens_mass_phase["Liq"]
                )
            sf3 = sf3 / len(self.control_volume.properties_in)
            sf4 = iscale.get_scaling_factor(self.single_paddle_length)
            sf5 = iscale.get_scaling_factor(self.single_paddle_width)
            sf6 = iscale.get_scaling_factor(self.num_paddle_wheels)
            sf7 = iscale.get_scaling_factor(self.num_paddles_per_wheel)
            sf8 = iscale.get_scaling_factor(self.paddle_drag_coef)
            sf = 0.5 * sf8 * (sf4 * sf5 * sf6 * sf7) * sf3 * sf2**3 * sf1**3 * 500
            iscale.set_scaling_factor(self.flocculation_power, sf)

            for ind, c in self.eq_flocculation_power.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling for total_power
        if iscale.get_scaling_factor(self.total_power[t]) is None:
            sf1 = iscale.get_scaling_factor(self.flocculation_power)
            sf2 = iscale.get_scaling_factor(self.rapid_mixing_power)
            sf = (sf1 + sf2) / 2
            iscale.set_scaling_factor(self.total_power, sf)

            for ind, c in self.eq_total_power.items():
                iscale.constraint_scaling_transform(c, sf)

        # set scaling factors for control_volume.properties_in based on control_volume.properties_out
        for t in self.control_volume.properties_in:
            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].dens_mass_phase
                )
                is None
            ):
                sf = iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].dens_mass_phase
                )
                iscale.set_scaling_factor(
                    self.control_volume.properties_in[t].dens_mass_phase, sf
                )

            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].flow_mass_phase_comp
                )
                is None
            ):
                for ind in self.control_volume.properties_in[t].flow_mass_phase_comp:
                    sf = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].flow_mass_phase_comp[ind]
                    )
                    iscale.set_scaling_factor(
                        self.control_volume.properties_in[t].flow_mass_phase_comp[ind],
                        sf,
                    )

            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].mass_frac_phase_comp
                )
                is None
            ):
                for ind in self.control_volume.properties_in[t].mass_frac_phase_comp:
                    sf = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].mass_frac_phase_comp[ind]
                    )
                    iscale.set_scaling_factor(
                        self.control_volume.properties_in[t].mass_frac_phase_comp[ind],
                        sf,
                    )

            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].flow_vol_phase
                )
                is None
            ):
                for ind in self.control_volume.properties_in[t].flow_vol_phase:
                    sf = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].flow_vol_phase[ind]
                    )
                    iscale.set_scaling_factor(
                        self.control_volume.properties_in[t].flow_vol_phase[ind], sf
                    )

        # update scaling for control_volume.properties_out
        for t in self.control_volume.properties_out:
            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].dens_mass_phase
                )
                is None
            ):
                iscale.set_scaling_factor(
                    self.control_volume.properties_out[t].dens_mass_phase, 1e-3
                )
            if (
                iscale.get_scaling_factor(
                    self.control_volume.properties_out[t].visc_d_phase
                )
                is None
            ):
                iscale.set_scaling_factor(
                    self.control_volume.properties_out[t].visc_d_phase, 1e3
                )

            # need to update scaling factors for TSS, Sludge, and TDS to account for the
            #   expected change in their respective values from the loss/gain rates
            for ind in self.control_volume.properties_out[t].flow_mass_phase_comp:
                if ind[1] == "TSS":
                    sf_og = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].flow_mass_phase_comp[ind]
                    )
                    sf_new = iscale.get_scaling_factor(self.tss_loss_rate)
                    iscale.set_scaling_factor(
                        self.control_volume.properties_out[t].flow_mass_phase_comp[ind],
                        100 * sf_new * (sf_new / sf_og),
                    )
                if ind[1] == "Sludge":
                    sf_og = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].flow_mass_phase_comp[ind]
                    )
                    sf_new = iscale.get_scaling_factor(self.tss_loss_rate)
                    iscale.set_scaling_factor(
                        self.control_volume.properties_out[t].flow_mass_phase_comp[ind],
                        100 * sf_new * (sf_new / sf_og),
                    )

            for ind in self.control_volume.properties_out[t].mass_frac_phase_comp:
                if ind[1] == "TSS":
                    sf_og = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].mass_frac_phase_comp[ind]
                    )
                    sf_new = iscale.get_scaling_factor(self.tss_loss_rate)
                    iscale.set_scaling_factor(
                        self.control_volume.properties_out[t].mass_frac_phase_comp[ind],
                        100 * sf_new * (sf_new / sf_og),
                    )
                if ind[1] == "Sludge":
                    sf_og = iscale.get_scaling_factor(
                        self.control_volume.properties_out[t].mass_frac_phase_comp[ind]
                    )
                    sf_new = iscale.get_scaling_factor(self.tss_loss_rate)
                    iscale.set_scaling_factor(
                        self.control_volume.properties_out[t].mass_frac_phase_comp[ind],
                        100 * sf_new * (sf_new / sf_og),
                    )

    def _get_performance_contents(self, time_point=0):
        t = time_point
        return {
            "vars": {
                "Total Power Usage  (kW)": self.total_power[t],
                "Rapid Mixing Power (kW)": self.rapid_mixing_power[t],
                "Flocc Mixing Power (kW)": self.flocculation_power[t],
            },
            "exprs": {},
            "params": {},
        }
