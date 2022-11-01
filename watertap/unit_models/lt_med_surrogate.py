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
    Suffix,
    NonNegativeReals,
    Reference,
    Constraint,
    units as pyunits,
    exp,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


# BPE functions
def BPE(dT, dS):
    T = dT
    S = dS / 1000
    a1, a2, a3, a4, a5, a6 = (
        -0.00045838530457,
        0.28230948284,
        17.945189194,
        0.00015361752708,
        0.052669058133,
        6.5604855793,
    )
    AA = a1 * T**2 + a2 * T + a3
    BB = a4 * T**2 + a5 * T + a6
    SW_BPE = AA * S**2 + BB * S

    return SW_BPE


# TO-DO:
# Add constraints to the input variables


@declare_process_block_class("LT_MED_surrogate")
class LTMEDData(UnitModelBlockData):
    """
    Zero order filtration model
    """

    CONFIG = ConfigBlock()

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
        "property_package2",
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
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        """
        Add system configurations
        """
        self.Nef = Param(
            initialize=12, units=pyunits.dimensionless, doc="Number of effects"
        )

        self.RR = Var(
            initialize=0.50,
            bounds=(0.30, 0.50),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        self.Q_loss = Param(
            initialize=0.054, units=pyunits.dimensionless, doc="System thermal loss"
        )

        """
        Add block for feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict
        )

        # Set alias for feed water salinity and convert to mg/L for surrogate model
        self.Xf = pyunits.convert(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.mg / pyunits.L,
        )

        # Set alias for feed water temperature and convert to degree C for surrogate model
        self.Tin = self.feed_props[0].temperature - 273.15 * pyunits.K

        # Set alias for feed water enthalpy and convert to kJ/kg
        self.h_sw = pyunits.convert(
            self.feed_props[0].enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
        )

        # Set alias for feed water mass flow rate (kg/s)
        self.m_f = sum(
            self.feed_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.feed_props.component_list
        )

        """
        Add block for distillate
        """
        self.distillate_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of distillate",
            **tmp_dict
        )

        # distillate salinity is 0
        @self.Constraint(doc="distillate salinity")
        def distillate_s(b):
            return b.distillate_props[0].conc_mass_phase_comp["Liq", "TDS"] == 0

        # Connect the plant capacity with the distillate mass flow rate
        self.Capacity = Var(
            initialize=2000,
            bounds=(2000, None),
            units=pyunits.m**3 / pyunits.d,
            doc="Capacity of the plant (m3/day)",
        )

        @self.Constraint(doc="distillate mass flow rate")
        def distillate_mfr(b):
            distillate_mfr = b.distillate_props[0].flow_vol_phase["Liq"]
            return b.Capacity == pyunits.convert(
                distillate_mfr, to_units=pyunits.m**3 / pyunits.d
            )

        # distillate temperature is the same as last effect vapor temperature, which is 10 deg higher than condenser inlet seawater temperature
        @self.Constraint(doc="distillate temperature")
        def distillate_temp(b):
            return (
                b.distillate_props[0].temperature == b.Tin + (273.15 + 10) * pyunits.K
            )

        # Set alias for distillate temperature and convert to degree C
        # To-do: may not need T_d anymore
        self.T_d = self.distillate_props[0].temperature - 273.15 * pyunits.K

        # Set alias for distillate enthalpy and convert to kJ/kg
        self.h_d = pyunits.convert(
            self.distillate_props[0].enth_mass_phase["Liq"],
            to_units=pyunits.kJ / pyunits.kg,
        )

        # Set alias for distillate volume (m3/hr) and mass (kg/s) flow rate
        self.m_d = sum(
            self.distillate_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.distillate_props.component_list
        )
        self.q_d = pyunits.convert(
            self.distillate_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

        self.brine_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        # Relationship between brine salinity and feed salinity
        @self.Constraint(doc="Brine salinity")
        def s_b_cal(b):
            return b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"] == b.feed_props[
                0
            ].conc_mass_phase_comp["Liq", "TDS"] / (1 - b.RR)

        # Brine temperature
        @self.Constraint(doc="Brine temperature")
        def T_b_cal(b):
            return (
                b.brine_props[0].temperature
                == b.T_d
                + (
                    273.15
                    + BPE(b.T_d, b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"])
                )
                * pyunits.K
            )

        # Set alias for brine temperature and convert to deg C
        self.T_b = self.brine_props[0].temperature - 273.15 * pyunits.K
        # Set alias for brine enthalpy and convert to kJ/kg
        self.h_b = pyunits.convert(
            self.brine_props[0].enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
        )
        # Set alias for brine volume (m3/hr) and mass (kg/s) flow rate
        self.m_b = sum(
            self.brine_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.brine_props.component_list
        )
        self.q_b = pyunits.convert(
            self.brine_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        """
        Add block for reject cooling water
        """
        self.cooling_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of cooling reject",
            **tmp_dict
        )

        # Use the source water for cooling (Same salinity)
        @self.Constraint(doc="Cooling reject salinity")
        def s_cooling_cal(b):
            return (
                b.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
            )

        # Set alias for enthalpy and convert to kJ/kg
        self.h_cool = pyunits.convert(
            self.cooling_out_props[0].enth_mass_phase["Liq"],
            to_units=pyunits.kJ / pyunits.kg,
        )
        # Assumption: the temperature of cooling reject is 3 degC lower than in the condenser
        @self.Constraint(doc="Cooling reject temperature")
        def T_cool_cal(b):
            return (
                b.cooling_out_props[0].temperature
                == b.distillate_props[0].temperature - 3 * pyunits.K
            )

        # Set alias for cooling reject volume flow rate (m3/hr)
        self.q_cooling = pyunits.convert(
            self.cooling_out_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        """
        Add block for heating steam
        """
        tmp_dict["parameters"] = self.config.property_package2
        tmp_dict["defined_state"] = True

        self.steam_props = self.config.property_package2.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of heating steam",
            **tmp_dict
        )

        # Set alias for temperature and convert to degree C for surrogate model
        self.Ts = self.steam_props[0].temperature - 273.15 * pyunits.K

        # All in steam
        @self.Constraint(doc="Steam inlet")
        def steam_phase(b):
            return b.steam_props[0].flow_mass_phase_comp["Liq", "H2O"] == 0

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="distillate", block=self.distillate_props)
        self.add_port(name="brine", block=self.brine_props)

        """
        Mass balances
        """
        # Feed flow rate calculation
        @self.Constraint(doc="Feed water volume flow rate")
        def qF_cal(b):
            return b.feed_props[0].flow_vol_phase["Liq"] == pyunits.convert(
                b.Capacity / b.RR, to_units=pyunits.m**3 / pyunits.s
            )

        # Set alias for feed volume flow rate and covert to m3/hr
        self.qF = pyunits.convert(
            self.feed_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        # Brine flow rate calculation
        @self.Constraint(doc="Brine volume flow rate")
        def q_b_cal(b):
            return (
                b.brine_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                - b.distillate_props[0].flow_vol_phase["Liq"]
            )

        """
        Add Vars for model outputs
        """
        self.P_req = Var(
            initialize=5000,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW)",
        )

        self.sA = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**2 / (pyunits.m**3 / pyunits.d),
            doc="Specific area (m2/m3/day))",
        )

        self.STEC = Var(
            initialize=65,
            bounds=(0, None),
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific thermal power consumption (kWh/m3)",
        )

        self.GOR = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg / pyunits.kg,
            doc="Gained output ratio (kg of distillate water per kg of heating steam",
        )

        self.qs = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.kg / pyunits.h,
            doc="Steam flow rate (kg/h)",
        )

        """
        Add Vars for intermediate model variables
        """
        self.TN = Var(
            initialize=35,
            bounds=(25, 45),
            units=pyunits.C,
            doc="Last effect vapor temperature",
        )

        self.m_sw = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.kg / pyunits.s,
            doc="Feed and cooling water mass flow rate (kg/s)",
        )

        self.q_sw = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.h,
            doc="Feed and cooling water volume flow rate (m3/h)",
        )

        # Surrogate equations for calculating GOR
        Nef_vals = [3, 6, 9, 12, 14]  # Number of effects
        coeffs_set = range(15)  # Number of coefficients

        GOR_coeffs = {
            3: [
                1.60e-07,
                0.826895712,
                -2.04e-07,
                0.003340838,
                -5.56e-09,
                0.000666667,
                -0.003295958,
                1.17e-10,
                -0.000549708,
                -2.46e-06,
                2.662545127,
                -1.98e-07,
                7.41e-07,
                -0.675925926,
                -4.12e-13,
            ],
            6: [
                5.86e-07,
                2.940942982,
                -1.44e-06,
                0.007985234,
                -2.26e-08,
                0.001472222,
                -0.007157144,
                8.73e-09,
                -0.001540936,
                4.88e-06,
                4.741753363,
                -1.04e-06,
                -1.67e-05,
                -2.333333333,
                -7.41e-13,
            ],
            9: [
                1.67e-06,
                5.9507846,
                -3.94e-06,
                0.012607651,
                -5.50e-08,
                0.002222222,
                -0.010203548,
                2.47e-08,
                -0.002549708,
                2.94e-05,
                6.350104873,
                -1.05e-05,
                -5.09e-05,
                -4.653703704,
                -2.88e-12,
            ],
            12: [
                3.30e-06,
                9.621851852,
                -7.98e-06,
                0.016637037,
                -1.09e-07,
                0.002,
                -0.012637326,
                5.13e-08,
                -0.003277778,
                8.28e-05,
                7.592772368,
                -3.09e-05,
                -9.98e-05,
                -7.425925926,
                -6.09e-12,
            ],
            14: [
                5.27e-06,
                12.44928443,
                -1.15e-05,
                0.019398098,
                -1.59e-07,
                0.001666667,
                -0.013396636,
                7.88e-08,
                -0.003333333,
                0.000121663,
                8.195669495,
                -5.34e-05,
                -0.00013251,
                -9.627577763,
                -1.80e-11,
            ],
        }

        @self.Constraint(doc="GOR surrogate equation")
        def GOR_cal(b):
            return (
                b.GOR
                == b.Xf * GOR_coeffs[b.Nef.value][0]
                + b.RR * GOR_coeffs[b.Nef.value][1]
                + b.Xf * b.RR * GOR_coeffs[b.Nef.value][2]
                + b.TN * GOR_coeffs[b.Nef.value][3]
                + b.TN * b.Xf * GOR_coeffs[b.Nef.value][4]
                + b.TN * b.RR * GOR_coeffs[b.Nef.value][5]
                + b.Ts * GOR_coeffs[b.Nef.value][6]
                + b.Ts * b.Xf * GOR_coeffs[b.Nef.value][7]
                + b.Ts * b.RR * GOR_coeffs[b.Nef.value][8]
                + b.Ts * b.TN * GOR_coeffs[b.Nef.value][9]
                + 1 * GOR_coeffs[b.Nef.value][10]
                + b.Ts**2 * GOR_coeffs[b.Nef.value][11]
                + b.TN**2 * GOR_coeffs[b.Nef.value][12]
                + b.RR**2 * GOR_coeffs[b.Nef.value][13]
                + b.Xf**2 * GOR_coeffs[b.Nef.value][14]
            )

        # Surrogate equations for calculating sA
        Nef_vals1 = [3, 6, 9]
        Nef_vals2 = [12, 14]

        coeffs_set1 = range(35)  # Number of coefficients for N=3, 6, 9
        coeffs_set2 = range(70)  # Number of coefficients for N=12, 14

        sA_coeffs = {
            3: [
                0.000596217,
                -3.66e-09,
                0,
                -2.44e-05,
                1.93e-09,
                0,
                5.60e-05,
                0,
                -2.95e-07,
                5.30e-11,
                0,
                7.14e-06,
                0,
                0.064807392,
                2.06e-07,
                0.00974051,
                0,
                -1.16e-05,
                -3.96e-11,
                0,
                -5.27e-06,
                0,
                -0.05718687,
                -2.61e-07,
                -0.011936049,
                -0.000702529,
                0.013464849,
                1.65e-07,
                0.003686623,
                0.000759933,
                0,
                -0.00019293,
                -0.000182949,
                0,
                3.20e-14,
            ],
            6: [
                0.00040105,
                -6.57e-09,
                0,
                -1.56e-05,
                3.67e-10,
                0,
                2.62e-05,
                0,
                7.08e-07,
                8.73e-12,
                0,
                1.46e-06,
                0,
                0.032775092,
                5.04e-08,
                0.002499309,
                0,
                -3.30e-06,
                -6.53e-12,
                0,
                -1.02e-06,
                0,
                -0.028641745,
                -6.66e-08,
                -0.002735652,
                -0.000301667,
                0.005600544,
                4.10e-08,
                0.000713386,
                0.000329733,
                0,
                -7.31e-05,
                -0.00010089,
                0,
                4.94e-14,
            ],
            9: [
                0.000596217,
                -3.66e-09,
                0,
                -2.44e-05,
                1.93e-09,
                0,
                5.60e-05,
                0,
                -2.95e-07,
                5.30e-11,
                0,
                7.14e-06,
                0,
                0.064807392,
                2.06e-07,
                0.00974051,
                0,
                -1.16e-05,
                -3.96e-11,
                0,
                -5.27e-06,
                0,
                -0.05718687,
                -2.61e-07,
                -0.011936049,
                -0.000702529,
                0.013464849,
                1.65e-07,
                0.003686623,
                0.000759933,
                0,
                -0.00019293,
                -0.000182949,
                0,
                3.20e-14,
            ],
            12: [
                0.000000e00,
                3.304374e-08,
                -6.761157e-13,
                0.000000e00,
                0.000000e00,
                -5.496094e-09,
                1.958695e-13,
                0.000000e00,
                0.000000e00,
                6.105760e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -8.520444e-10,
                2.227182e-14,
                0.000000e00,
                0.000000e00,
                6.610470e-10,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                3.561099e-06,
                9.693068e-12,
                0.000000e00,
                1.148815e-06,
                0.000000e00,
                0.000000e00,
                -1.049434e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                1.734819e-10,
                -1.367980e-14,
                0.000000e00,
                0.000000e00,
                -5.044097e-10,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -2.717657e-06,
                -3.905605e-11,
                0.000000e00,
                -1.397796e-06,
                0.000000e00,
                0.000000e00,
                -4.132341e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                4.380618e-07,
                2.072263e-11,
                0.000000e00,
                4.398758e-07,
                0.000000e00,
                0.000000e00,
                5.991695e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -1.833849e-08,
                0.000000e00,
                -3.505028e-06,
                0.000000e00,
                1.713226e-06,
                0.000000e00,
                0.000000e00,
                6.273434e-18,
            ],
            14: [
                0.000000e00,
                4.368251e-08,
                2.260942e-13,
                0.000000e00,
                0.000000e00,
                -4.762111e-08,
                1.282504e-12,
                0.000000e00,
                0.000000e00,
                3.611544e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -3.197411e-09,
                8.656950e-14,
                0.000000e00,
                0.000000e00,
                3.617982e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                6.034818e-06,
                4.276140e-11,
                0.000000e00,
                4.908627e-06,
                0.000000e00,
                0.000000e00,
                -1.357859e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                1.131144e-10,
                -5.279738e-14,
                0.000000e00,
                0.000000e00,
                -2.836528e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -3.906973e-06,
                -1.642642e-10,
                0.000000e00,
                -6.562415e-06,
                0.000000e00,
                0.000000e00,
                -1.071227e-07,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                7.876487e-07,
                8.941839e-11,
                0.000000e00,
                2.276216e-06,
                0.000000e00,
                0.000000e00,
                1.772933e-07,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -6.455677e-08,
                0.000000e00,
                -1.423548e-05,
                0.000000e00,
                7.012498e-06,
                0.000000e00,
                0.000000e00,
                1.716278e-19,
            ],
        }

        @self.Constraint(doc="sA surrogate equation")
        def sA_cal(b):
            if b.Nef.value in [3, 6, 9]:
                return (
                    b.sA
                    == b.Xf * sA_coeffs[b.Nef.value][0]
                    + b.Xf**2 * sA_coeffs[b.Nef.value][1]
                    + b.RR * sA_coeffs[b.Nef.value][2]
                    + b.RR * b.Xf * sA_coeffs[b.Nef.value][3]
                    + b.RR * b.Xf**2 * sA_coeffs[b.Nef.value][4]
                    + b.RR**2 * sA_coeffs[b.Nef.value][5]
                    + b.RR**2 * b.Xf * sA_coeffs[b.Nef.value][6]
                    + b.TN * sA_coeffs[b.Nef.value][7]
                    + b.TN * b.Xf * sA_coeffs[b.Nef.value][8]
                    + b.TN * b.Xf**2 * sA_coeffs[b.Nef.value][9]
                    + b.TN * b.RR * sA_coeffs[b.Nef.value][10]
                    + b.TN * b.RR * b.Xf * sA_coeffs[b.Nef.value][11]
                    + b.TN * b.RR**2 * sA_coeffs[b.Nef.value][12]
                    + b.TN**2 * sA_coeffs[b.Nef.value][13]
                    + b.TN**2 * b.Xf * sA_coeffs[b.Nef.value][14]
                    + b.TN**2 * b.RR * sA_coeffs[b.Nef.value][15]
                    + b.Ts * sA_coeffs[b.Nef.value][16]
                    + b.Ts * b.Xf * sA_coeffs[b.Nef.value][17]
                    + b.Ts * b.Xf**2 * sA_coeffs[b.Nef.value][18]
                    + b.Ts * b.RR * sA_coeffs[b.Nef.value][19]
                    + b.Ts * b.RR * b.Xf * sA_coeffs[b.Nef.value][20]
                    + b.Ts * b.RR**2 * sA_coeffs[b.Nef.value][21]
                    + b.Ts * b.TN * sA_coeffs[b.Nef.value][22]
                    + b.Ts * b.TN * b.Xf * sA_coeffs[b.Nef.value][23]
                    + b.Ts * b.TN * b.RR * sA_coeffs[b.Nef.value][24]
                    + b.Ts * b.TN**2 * sA_coeffs[b.Nef.value][25]
                    + b.Ts**2 * sA_coeffs[b.Nef.value][26]
                    + b.Ts**2 * b.Xf * sA_coeffs[b.Nef.value][27]
                    + b.Ts**2 * b.RR * sA_coeffs[b.Nef.value][28]
                    + b.Ts**2 * b.TN * sA_coeffs[b.Nef.value][29]
                    + 1 * sA_coeffs[b.Nef.value][30]
                    + b.Ts**3 * sA_coeffs[b.Nef.value][31]
                    + b.TN**3 * sA_coeffs[b.Nef.value][32]
                    + b.RR**3 * sA_coeffs[b.Nef.value][33]
                    + b.Xf**3 * sA_coeffs[b.Nef.value][34]
                )

            else:  #  Nef in [12, 14]
                return (
                    b.sA
                    == b.Xf * sA_coeffs[b.Nef.value][0]
                    + b.Xf**2 * sA_coeffs[b.Nef.value][1]
                    + b.Xf**3 * sA_coeffs[b.Nef.value][2]
                    + b.RR * sA_coeffs[b.Nef.value][3]
                    + b.RR * b.Xf * sA_coeffs[b.Nef.value][4]
                    + b.RR * b.Xf**2 * sA_coeffs[b.Nef.value][5]
                    + b.RR * b.Xf**3 * sA_coeffs[b.Nef.value][6]
                    + b.RR**2 * sA_coeffs[b.Nef.value][7]
                    + b.RR**2 * b.Xf * sA_coeffs[b.Nef.value][8]
                    + b.RR**2 * b.Xf**2 * sA_coeffs[b.Nef.value][9]
                    + b.RR**3 * sA_coeffs[b.Nef.value][10]
                    + b.RR**3 * b.Xf * sA_coeffs[b.Nef.value][11]
                    + b.TN * sA_coeffs[b.Nef.value][12]
                    + b.TN * b.Xf * sA_coeffs[b.Nef.value][13]
                    + b.TN * b.Xf**2 * sA_coeffs[b.Nef.value][14]
                    + b.TN * b.Xf**3 * sA_coeffs[b.Nef.value][15]
                    + b.TN * b.RR * sA_coeffs[b.Nef.value][16]
                    + b.TN * b.RR * b.Xf * sA_coeffs[b.Nef.value][17]
                    + b.TN * b.RR * b.Xf**2 * sA_coeffs[b.Nef.value][18]
                    + b.TN * b.RR**2 * sA_coeffs[b.Nef.value][19]
                    + b.TN * b.RR**2 * b.Xf * sA_coeffs[b.Nef.value][20]
                    + b.TN * b.RR**3 * sA_coeffs[b.Nef.value][21]
                    + b.TN**2 * sA_coeffs[b.Nef.value][22]
                    + b.TN**2 * b.Xf * sA_coeffs[b.Nef.value][23]
                    + b.TN**2 * b.Xf**2 * sA_coeffs[b.Nef.value][24]
                    + b.TN**2 * b.RR * sA_coeffs[b.Nef.value][25]
                    + b.TN**2 * b.RR * b.Xf * sA_coeffs[b.Nef.value][26]
                    + b.TN**2 * b.RR**2 * sA_coeffs[b.Nef.value][27]
                    + b.TN**3 * sA_coeffs[b.Nef.value][28]
                    + b.TN**3 * b.Xf * sA_coeffs[b.Nef.value][29]
                    + b.TN**3 * b.RR * sA_coeffs[b.Nef.value][30]
                    + b.Ts * sA_coeffs[b.Nef.value][31]
                    + b.Ts * b.Xf * sA_coeffs[b.Nef.value][32]
                    + b.Ts * b.Xf**2 * sA_coeffs[b.Nef.value][33]
                    + b.Ts * b.Xf**3 * sA_coeffs[b.Nef.value][34]
                    + b.Ts * b.RR * sA_coeffs[b.Nef.value][35]
                    + b.Ts * b.RR * b.Xf * sA_coeffs[b.Nef.value][36]
                    + b.Ts * b.RR * b.Xf**2 * sA_coeffs[b.Nef.value][37]
                    + b.Ts * b.RR**2 * sA_coeffs[b.Nef.value][38]
                    + b.Ts * b.RR**2 * b.Xf * sA_coeffs[b.Nef.value][39]
                    + b.Ts * b.RR**3 * sA_coeffs[b.Nef.value][40]
                    + b.Ts * b.TN * sA_coeffs[b.Nef.value][41]
                    + b.Ts * b.TN * b.Xf * sA_coeffs[b.Nef.value][42]
                    + b.Ts * b.TN * b.Xf**2 * sA_coeffs[b.Nef.value][43]
                    + b.Ts * b.TN * b.RR * sA_coeffs[b.Nef.value][44]
                    + b.Ts * b.TN * b.RR * b.Xf * sA_coeffs[b.Nef.value][45]
                    + b.Ts * b.TN * b.RR**2 * sA_coeffs[b.Nef.value][46]
                    + b.Ts * b.TN**2 * sA_coeffs[b.Nef.value][47]
                    + b.Ts * b.TN**2 * b.Xf * sA_coeffs[b.Nef.value][48]
                    + b.Ts * b.TN**2 * b.RR * sA_coeffs[b.Nef.value][49]
                    + b.Ts * b.TN**3 * sA_coeffs[b.Nef.value][50]
                    + b.Ts**2 * sA_coeffs[b.Nef.value][51]
                    + b.Ts**2 * b.Xf * sA_coeffs[b.Nef.value][52]
                    + b.Ts**2 * b.Xf**2 * sA_coeffs[b.Nef.value][53]
                    + b.Ts**2 * b.RR * sA_coeffs[b.Nef.value][54]
                    + b.Ts**2 * b.RR * b.Xf * sA_coeffs[b.Nef.value][55]
                    + b.Ts**2 * b.RR**2 * sA_coeffs[b.Nef.value][56]
                    + b.Ts**2 * b.TN * sA_coeffs[b.Nef.value][57]
                    + b.Ts**2 * b.TN * b.Xf * sA_coeffs[b.Nef.value][58]
                    + b.Ts**2 * b.TN * b.RR * sA_coeffs[b.Nef.value][59]
                    + b.Ts**2 * b.TN**2 * sA_coeffs[b.Nef.value][60]
                    + b.Ts**3 * sA_coeffs[b.Nef.value][61]
                    + b.Ts**3 * b.Xf * sA_coeffs[b.Nef.value][62]
                    + b.Ts**3 * b.RR * sA_coeffs[b.Nef.value][63]
                    + b.Ts**3 * b.TN * sA_coeffs[b.Nef.value][64]
                    + 1 * sA_coeffs[b.Nef.value][65]
                    + b.Ts**4 * sA_coeffs[b.Nef.value][66]
                    + b.TN**4 * sA_coeffs[b.Nef.value][67]
                    + b.RR**4 * sA_coeffs[b.Nef.value][68]
                    + b.Xf**4 * sA_coeffs[b.Nef.value][69]
                )

        # Last effect vapor temperature is 10 degree higher than condenser inlet temperature
        @self.Constraint(doc="System configuration 1")
        def TN_Tin(b):
            return b.TN == b.Tin + 10

        # Steam flow rate calculation
        @self.Constraint(doc="Steam flow rate")
        def qs_cal(b):
            return b.steam_props[0].flow_mass_phase_comp[
                "Vap", "H2O"
            ] == pyunits.convert(
                sum(
                    b.distillate_props[0].flow_mass_phase_comp["Liq", j]
                    for j in b.distillate_props.component_list
                )
                / b.GOR,
                to_units=pyunits.kg / pyunits.s,
            )

        # Energy consumption
        @self.Constraint(doc="STEC calculation")
        def STEC_cal(b):
            return b.STEC == pyunits.convert(
                1
                / b.GOR
                * (b.steam_props[0].dh_vap_mass)
                * b.distillate_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(doc="Thermal power requirement calculation")
        def P_req_cal(b):
            return b.P_req == pyunits.convert(b.STEC * b.Capacity, to_units=pyunits.kW)

        # Mass flow rate
        @self.Constraint(doc="Feed and cooling water mass flow rate (kg/s)")
        def m_sw_cal(b):
            return b.m_sw * (b.h_cool - b.h_sw) == (
                (1 - b.Q_loss) * b.P_req
                - b.m_b * b.h_b
                - b.m_d * b.h_d
                + b.h_cool * b.m_f
            )

        # Volume flow rates
        @self.Constraint(doc="Feed and cooling water mass flow rate (m3/h)")
        def q_sw_cal(b):
            return b.q_sw == pyunits.convert(
                b.m_sw / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hour,
            )

        @self.Constraint(doc="Cooling water mass flow rate (m3/h)")
        def q_cooling_cal(b):
            return b.q_cooling == b.q_sw - b.qF

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.RR) is None:
            iscale.set_scaling_factor(self.RR, 1e1)

        if iscale.get_scaling_factor(self.Capacity) is None:
            # sf = iscale.get_scaling_factor(self.Xf, default=1, warning=True)
            iscale.set_scaling_factor(self.Capacity, 1e-3)

        if iscale.get_scaling_factor(self.Xf) is None:
            iscale.set_scaling_factor(self.Xf, 1e5)

        if iscale.get_scaling_factor(self.Tin) is None:
            iscale.set_scaling_factor(self.Tin, 1e2)

        if iscale.get_scaling_factor(self.h_sw) is None:
            iscale.set_scaling_factor(self.h_sw, 1e-3)

        if iscale.get_scaling_factor(self.m_f) is None:
            iscale.set_scaling_factor(self.m_f, 1e-2)

        if iscale.get_scaling_factor(self.P_req) is None:
            iscale.set_scaling_factor(
                self.P_req,
                iscale.get_scaling_factor(self.Capacity)
                * iscale.get_scaling_factor(self.STEC),
            )

        if iscale.get_scaling_factor(self.sA) is None:
            iscale.set_scaling_factor(self.sA, 1e-1)

        if iscale.get_scaling_factor(self.STEC) is None:
            iscale.set_scaling_factor(self.STEC, 1e-3)

        if iscale.get_scaling_factor(self.GOR) is None:
            iscale.set_scaling_factor(self.GOR, 1e-1)

        if iscale.get_scaling_factor(self.qs) is None:
            iscale.set_scaling_factor(self.qs, 1e-2)
