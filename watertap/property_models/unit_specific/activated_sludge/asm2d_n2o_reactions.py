#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
ASM2d-PSFe-GHG reaction package.

Reference:
[1] B. Solis, A. Guisasola, X. Flores-Alsina, U. Jeppsson, J.A. Baeza,
A plant-wide model describing GHG emissions and nutrient recovery options for water resource recovery facilities,
Water Research 215 (2022) https://www.sciencedirect.com/science/article/pii/S0043135422001865
"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import Bool, ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BurntToast

import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme


# Some more information about this module
__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASM2dN2OReactionParameterBlock")
class ASM2dN2OReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    CONFIG = ReactionParameterBlock.CONFIG()

    CONFIG.declare(
        "decay_switch",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Switching function for decay",
            doc="""Switching function for handling decay in reaction rate expressions,
**default** - True.
**Valid values:** {
**True** - the decay of heterotrophs and autotrophs is dependent on the electron acceptor present,
**False** - the decay of heterotrophs and autotrophs does not change""",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ASM2dGHGReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ASM2d paper
        # R1: Aerobic hydrolysis
        # R2: Anoxic hydrolysis (NO3-)
        # R3: Anoxic hydrolysis (NO2-)
        # R4: Anaerobic hydrolysis
        # R5: Aerobic growth on S_F
        # R6: Aerobic growth on S_A
        # R7: Anoxic growth on S_F (NO3- -> NO2-)
        # R8: Anoxic growth on S_F (NO2- -> NO)
        # R9: Anoxic growth on S_F (NO -> N2O)
        # R10: Anoxic growth on S_F (N2O- -> N2)
        # R11: Anoxic growth on S_A (NO3- -> NO2-)
        # R12: Anoxic growth on S_A (NO2- -> NO)
        # R13: Anoxic growth on S_A (NO -> N2O)
        # R14: Anoxic growth on S_A (N2O- -> N2)
        # R15: Fermentation
        # R16: Lysis of X_H
        # R17: Storage of X_PHA
        # R18: Aerobic storage of X_PP
        # R19: Anoxic storage of X_PP (NO3- -> NO2-)
        # R20: Anoxic storage of X_PP (NO2- -> NO)
        # R21: Anoxic storage of X_PP (NO -> N2O)
        # R22: Anoxic storage of X_PP (N2O- -> N2)
        # R23: Aerobic growth of X_PAO
        # R24: Anoxic growth of X_PAO (NO3- -> NO2-)
        # R25: Anoxic growth of X_PAO (NO2- -> NO)
        # R26: Anoxic growth of X_PAO (NO -> N2O)
        # R27: Anoxic growth of X_PAO (N2O- -> N2)
        # R28: Lysis of X_PAO
        # R29: Lysis of X_PP
        # R30: Lysis of X_PHA
        # R31: NH3 oxidation to NH2OH with O2 consumption
        # R32: NH2OH to NO coupled with O2 reduction (X_AOB growth)
        # R33: NO oxidation to NO2- coupled with O2 reduction
        # R34: NO to N2O coupled with NH2OH to No2- (N2O from NN pathway)
        # R35: HNO2 to N2O coupled with NH2OH to NO2- (N2O from ND pathway)
        # R36: Aerobic growth of XNOB
        # R37: Lysis of X_AOB
        # R38: Lysis of X_NOB
        self.rate_reaction_idx = pyo.Set(
            initialize=[
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
                "R13",
                "R14",
                "R15",
                "R16",
                "R17",
                "R18",
                "R19",
                "R20",
                "R21",
                "R22",
                "R23",
                "R24",
                "R25",
                "R26",
                "R27",
                "R28",
                "R29",
                "R30",
                "R31",
                "R32",
                "R33",
                "R34",
                "R35",
                "R36",
                "R37",
                "R38",
            ]
        )

        # Stoichiometric Parameters
        self.f_SI = pyo.Var(
            initialize=0.00,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Production of S_I in hydrolysis, [kg COD/kg COD]",
        )
        self.Y_H = pyo.Var(
            initialize=0.625,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield coefficient for heterotrophic biomass, [kg COD/kg COD]",
        )
        self.Y_PAO = pyo.Var(
            initialize=0.625,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield coefficient for P accumulating organisms (biomass/PHA), [kg COD/kg COD]",
        )
        self.Y_PO4 = pyo.Var(
            initialize=0.4,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="PP requirement (PO4 release) per PHA stored, [kg P/kg COD]",
        )
        self.Y_PHA = pyo.Var(
            initialize=0.20,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="PHA requirement for PP storage, [kg COD/kg P]",
        )
        self.Y_AOB = pyo.Var(
            initialize=0.18,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield of ammonia oxidizing bacteria, [kg COD/kg N]",
        )
        self.Y_NOB = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield of nitrite oxidizing bacteria, [kg COD/kg N]",
        )
        self.f_XIH = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD from lysis, [kg COD/kg COD]",
        )
        self.f_XIP = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD from lysis, [kg COD/kg COD]",
        )
        self.f_XIA = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD from lysis, [kg COD/kg COD]",
        )
        self.nG = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic growth factor",
        )
        self.i_CSI = pyo.Var(
            initialize=0.36178,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD S_I, [kg C/kg COD]",
        )
        self.i_CSF = pyo.Var(
            initialize=0.31843,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD S_F, [kg C/kg COD]",
        )
        self.i_CSA = pyo.Var(
            initialize=0.37500,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD S_A, [kg C/kg COD]",
        )
        self.i_CXI = pyo.Var(
            initialize=0.36178,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD X_I, [kg C/kg COD]",
        )
        self.i_CXS = pyo.Var(
            initialize=0.31843,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD X_S, [kg C/kg COD]",
        )
        self.i_CBM = pyo.Var(
            initialize=0.36612,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of biomass, [kg C/kg COD]",
        )
        self.i_CXPHA = pyo.Var(
            initialize=0.3,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of XPHA, [kg C/kg COD]",
        )
        self.i_NSI = pyo.Var(
            initialize=0.06003,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert soluble COD S_I, [kg N/kg COD]",
        )
        self.i_NSF = pyo.Var(
            initialize=0.03352,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of fermentable substrate, S_F, [kg N/kg COD]",
        )
        self.i_NXI = pyo.Var(
            initialize=0.06003,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert particulate COD X_I, [kg N/kg COD]",
        )
        self.i_NXS = pyo.Var(
            initialize=0.03352,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of slowly biodegradable substrate X_S, [kg N/kg COD]",
        )
        self.i_NBM = pyo.Var(
            initialize=0.08615,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of biomass, X_H, X_PAO, X_AUT, [kg N/kg COD]",
        )
        self.i_PSI = pyo.Var(
            initialize=0.00649,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert soluble COD, [kg P/kg COD]",
        )
        self.i_PSF = pyo.Var(
            initialize=0.00559,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of fermentable substrate, S_F, [kg P/kg COD]",
        )
        self.i_PXI = pyo.Var(
            initialize=0.00649,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert particulate COD X_I, [kg P/kg COD]",
        )
        self.i_PXS = pyo.Var(
            initialize=0.00559,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of slowly biodegradable substrate X_S, [kg P/kg COD]",
        )
        self.i_PBM = pyo.Var(
            initialize=0.02154,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of biomass, X_H, X_PAO, X_AUT, [kg P/kg COD]",
        )

        self.i_KXPP = pyo.Var(
            initialize=0.4204,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Potassium coefficient for polyphosphates",
        )
        self.i_MgXPP = pyo.Var(
            initialize=0.2614,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Magnesium coefficient for polyphosphates",
        )

        # Kinetic Parameters
        self.K_H = pyo.Var(
            initialize=2.46,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Hydrolysis rate constant",
        )
        self.hL_NO3 = pyo.Var(
            initialize=0.60,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic hydrolysis reduction factor for nitrate",
        )
        self.hL_fe = pyo.Var(
            initialize=0.40,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anaerobic hydrolysis reduction factor",
        )
        self.KL_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen, [kg O2/m^3]",
        )
        self.KL_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.KL_X = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for particulate COD, [kg X_S/kg X_H]",
        )
        self.hL_NO2 = pyo.Var(
            initialize=0.60,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic hydrolysis reduction factor for nitrite",
        )
        self.KL_NO2 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrite, [kg N/m^3]",
        )
        self.mu_H = pyo.Var(
            initialize=4.23,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate on substrate, [kg X_S/kg X_H/day]",
        )
        self.q_fe = pyo.Var(
            initialize=2.11,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum rate for fermentation, [kg S_F/kg X_H/day]",
        )
        self.hH_NO3 = pyo.Var(
            initialize=0.28,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitrate to nitrite)",
        )
        self.b_H = pyo.Var(
            initialize=0.28,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for lysis and decay",
        )
        self.hH_NO3_end = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for endogenous respiration",
        )
        self.KH_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen, [kg O2/m^3]",
        )
        self.K_F = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on SF, [kg COD/m^3]",
        )
        self.K_fe = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for fermentation of SF, [kg COD/m^3]",
        )
        self.KH_A = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA, [kg COD/m^3]",
        )
        self.KH_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.KH_NH4 = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for ammonium (nutrient), [kg N/m^3]",
        )
        self.KH_PO4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for SPO4 as nutrient, [kg P/m^3]",
        )
        self.hH_NO2 = pyo.Var(
            initialize=0.16,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitrite to nitric oxide)",
        )
        self.hH_NO = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitric oxide to nitrous oxide)",
        )
        self.hH_N2O = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitrous oxide to dinitrogen)",
        )
        self.KH2_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen (nitrate to nitrite), [kg O2/m^3]",
        )
        self.KH3_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen (nitrite to nitric oxide), [kg O2/m^3]",
        )
        self.KH4_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen (nitric oxide to nitrous oxide), [kg O2/m^3]",
        )
        self.KH5_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen (nitrous oxide to dinitrogen), [kg O2/m^3]",
        )
        self.K_F2 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on SF (nitrate to nitrite), [kg COD/m^3]",
        )
        self.K_F3 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on SF (nitrite to nitric oxide), [kg COD/m^3]",
        )
        self.K_F4 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on SF (nitric oxide to nitrous oxide, [kg COD/m^3]",
        )
        self.K_F5 = pyo.Var(
            initialize=40e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on SF (nitrous oxide to dinitrogen), [kg COD/m^3]",
        )
        self.KH_A2 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA (nitrate to nitrite), [kg COD/m^3]",
        )
        self.KH_A3 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA (nitrite to nitric oxide), [kg COD/m^3]",
        )
        self.KH_A4 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA (nitric oxide to nitrous oxide), [kg COD/m^3]",
        )
        self.KH_A5 = pyo.Var(
            initialize=40e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA (nitrous oxide to dinitrogen), [kg COD/m^3]",
        )
        self.KH_NO2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrite, [kg COD/m^3]",
        )
        self.KH_NO = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitric oxide, [kg COD/m^3]",
        )
        self.KH_N2O = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrous oxide, [kg COD/m^3]",
        )
        self.KI_NO = pyo.Var(
            initialize=3e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Inhibition coefficient for nitric oxide, [kg COD/m^3]",
        )
        self.q_PHA = pyo.Var(
            initialize=2.46,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for storage of X_PHA (base Xpp), [kg PHA/kg PAO/day]",
        )
        self.q_PP = pyo.Var(
            initialize=1.23,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for storage of X_PP, [kg PP/kg PAO/day]",
        )
        self.mu_PAO = pyo.Var(
            initialize=0.82,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of XPAO",
        )
        self.hP_NO3 = pyo.Var(
            initialize=0.28,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor under anoxic conditions",
        )
        self.b_PAO = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PAO",
        )
        self.hP_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PAOs",
        )
        self.b_PP = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PP",
        )
        self.hPP_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PP",
        )
        self.b_PHA = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PHA",
        )
        self.hPHA_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PHA",
        )
        self.KP_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen, [kg O2/m^3]",
        )
        self.KP_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for nitrate, [kg N/m^3]",
        )
        self.KP_A = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for acetate, [kg COD/m^3]",
        )
        self.KP_NH4 = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for ammonium, [kg N/m^3]",
        )
        self.KP_P = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphate for XPP formation, [kg P/m^3]",
        )
        self.KP_PO4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphate for growth, [kg P/m^3]",
        )
        self.KP_PP = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for poly-phosphate, [kg PP/kg PAO]",
        )
        self.K_MAX = pyo.Var(
            initialize=0.34,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Maximum ratio of X_PP/X_PAO, [kg PP/kg PAO]",
        )
        self.KI_PP = pyo.Var(
            initialize=0.02,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Inhibition coefficient for PP storage, [kg PP/kg PAO]",
        )
        self.KP_PHA = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for PHA, [kg PHA/kg PAO]",
        )
        self.hP_NO2 = pyo.Var(
            initialize=0.16,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitrite to nitric oxide)",
        )
        self.hP_NO = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitric oxide to nitrous oxide)",
        )
        self.hP_N2O = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification (nitrous oxide to dinitrogen)",
        )
        self.KP2_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen (nitrate to nitrite), [kg O2/m^3]",
        )
        self.KP3_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen (nitrite to nitric oxide), [kg O2/m^3]",
        )
        self.KP4_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen (nitric oxide to nitrous oxide), [kg O2/m^3]",
        )
        self.KP5_O2 = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen (nitrous oxide to dinitrogen), [kg O2/m^3]",
        )
        self.KP_NO2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for nitrite, [kg O2/m^3]",
        )
        self.KP_NO = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for nitric oxide), [kg O2/m^3]",
        )
        self.KP_N2O = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for nitrous oxide), [kg O2/m^3]",
        )
        self.q_AOB_AOM = pyo.Var(
            initialize=5.2,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum rate for the AMO reaction",
        )
        self.mu_AOB_HAO = pyo.Var(
            initialize=0.61,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of AOB",
        )
        self.q_AOB_HAO = pyo.Var(
            initialize=5.2,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum rate for the HAO reaction",
        )
        self.q_AOB_N2O_NN = pyo.Var(
            initialize=0.0078,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum N2O production rate by NH2OH oxidation pathway",
        )
        self.q_AOB_N2O_ND = pyo.Var(
            initialize=1.3008,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum N2O production rate by ND pathway",
        )
        self.b_AOB = pyo.Var(
            initialize=0.096,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Decay rate of AOB",
        )
        self.hAOB_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for AOB decay",
        )
        self.KAOB1_O2 = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for oxygen (AMO reaction)",
        )
        self.KAOB_NH4 = pyo.Var(
            initialize=0.004,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for NH4",
        )
        self.KAOB2_O2 = pyo.Var(
            initialize=0.6,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for oxygen (HAO reaction)",
        )
        self.KAOB_NH2OH = pyo.Var(
            initialize=0.3,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for NH2OH",
        )
        self.KAOB_P = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for phosphate",
        )
        self.KAOB_HAO_NO = pyo.Var(
            initialize=0.0003,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for NO (from HAO)",
        )
        self.KAOB_NN_NO = pyo.Var(
            initialize=0.008,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for NO (form NirK)",
        )
        self.KAOB_HNO2 = pyo.Var(
            initialize=0.0006,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for free nitrous acid",
        )
        self.KAOB_ND_O2 = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB constant for O2 effect on the ND pathway",
        )
        self.KAOB_I_O2 = pyo.Var(
            initialize=0.8,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N2O constant for production inhibition by O2",
        )
        self.KAOB_NO3 = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="AOB affinity constant for nitrate",
        )
        self.mu_NOB = pyo.Var(
            initialize=0.61,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of NOB",
        )
        self.b_NOB = pyo.Var(
            initialize=0.096,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Decay rate of NOB",
        )
        self.hNOB_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for NOB decay",
        )
        self.KNOB_O2 = pyo.Var(
            initialize=1.2,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="NOB affinity constant for oxygen",
        )
        self.KNOB_NO2 = pyo.Var(
            initialize=1e-6,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="NOB affinity constant for nitrite",
        )
        self.KNOB_P = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="NOB affinity constant for phosphate",
        )
        self.KNOB_NO3 = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="NOB affinity constant for nitrate",
        )
        self.D_N2 = pyo.Var(
            initialize=1.86e-9,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of N2",
        )
        self.D_N2O = pyo.Var(
            initialize=1.77e-9,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of N2O",
        )
        self.D_NO = pyo.Var(
            initialize=8.7e-11,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of N2",
        )
        self.D_O2 = pyo.Var(
            initialize=2.12e-9,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of O2",
        )
        self.D_CO2 = pyo.Var(
            initialize=1.92e-9,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of CO2",
        )
        self.D_H2S = pyo.Var(
            initialize=1.41e-9,
            units=pyo.units.cm / pyo.units.s,
            domain=pyo.NonNegativeReals,
            doc="Diffusivity of H2S",
        )
        self.H_N2 = pyo.Var(
            initialize=6.4e-6,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of N2 at reference temperature (288.15 K)",
        )
        self.H_N2O = pyo.Var(
            initialize=2.4e-4,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of N2O at reference temperature (288.15 K)",
        )
        self.H_NO = pyo.Var(
            initialize=1.9e-5,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of NO at reference temperature (288.15 K)",
        )
        self.H_O2 = pyo.Var(
            initialize=1.3e-5,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of O2 at reference temperature (288.15 K)",
        )
        self.H_CO2 = pyo.Var(
            initialize=3.3e-4,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of CO2 at reference temperature (288.15 K)",
        )
        self.H_H2S = pyo.Var(
            initialize=1e-3,
            units=pyo.units.mol / pyo.units.m**3 / pyo.units.Pa,
            domain=pyo.NonNegativeReals,
            doc="Henry constant of H2S at reference temperature (288.15 K)",
        )
        self.P_N2 = pyo.Var(
            initialize=0.78,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of N2 in air",
        )
        self.P_N2O = pyo.Var(
            initialize=0,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of N2O in air",
        )
        self.P_NO = pyo.Var(
            initialize=0,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of NO in air",
        )
        self.P_O2 = pyo.Var(
            initialize=0.21,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of O2 in air",
        )
        self.P_CO2 = pyo.Var(
            initialize=3.947e-4,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of CO2 in air",
        )
        self.P_H2S = pyo.Var(
            initialize=0,
            units=pyo.units.atm,
            domain=pyo.NonNegativeReals,
            doc="Partial pressure of H2S in air",
        )
        # ----------------------------------------------------------------------------------------------------

        # Reaction Stoichiometry - see https://app.box.com/file/1812341593943
        # This is the stoichiometric part the Peterson matrix in dict form
        # Note that reaction stoichiometry is on a mass basis.

        self.rate_reaction_stoichiometry = {
            # R1: Aerobic hydrolysis
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_O2"): 0,
            ("R1", "Liq", "S_F"): 1 - self.f_SI,
            ("R1", "Liq", "S_A"): 0,
            ("R1", "Liq", "S_I"): self.f_SI,
            ("R1", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R1", "Liq", "S_NH2OH"): 0,
            ("R1", "Liq", "S_N2O"): 0,
            ("R1", "Liq", "S_NO"): 0,
            ("R1", "Liq", "S_NO2"): 0,
            ("R1", "Liq", "S_N2"): 0,
            ("R1", "Liq", "S_NO3"): 0,
            ("R1", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R1", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_S"): -1,
            ("R1", "Liq", "X_H"): 0,
            ("R1", "Liq", "X_PAO"): 0,
            ("R1", "Liq", "X_PP"): 0,
            ("R1", "Liq", "X_PHA"): 0,
            ("R1", "Liq", "X_AOB"): 0,
            ("R1", "Liq", "X_NOB"): 0,
            ("R1", "Liq", "S_K"): 0,
            ("R1", "Liq", "S_Mg"): 0,
            # R2: Anoxic hydrolysis (NO3-)
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_O2"): 0,
            ("R2", "Liq", "S_F"): 1 - self.f_SI,
            ("R2", "Liq", "S_A"): 0,
            ("R2", "Liq", "S_I"): self.f_SI,
            ("R2", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R2", "Liq", "S_NH2OH"): 0,
            ("R2", "Liq", "S_N2O"): 0,
            ("R2", "Liq", "S_NO"): 0,
            ("R2", "Liq", "S_NO2"): 0,
            ("R2", "Liq", "S_N2"): 0,
            ("R2", "Liq", "S_NO3"): 0,
            ("R2", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R2", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_S"): -1,
            ("R2", "Liq", "X_H"): 0,
            ("R2", "Liq", "X_PAO"): 0,
            ("R2", "Liq", "X_PP"): 0,
            ("R2", "Liq", "X_PHA"): 0,
            ("R2", "Liq", "X_AOB"): 0,
            ("R2", "Liq", "X_NOB"): 0,
            ("R2", "Liq", "S_K"): 0,
            ("R2", "Liq", "S_Mg"): 0,
            # R3: Anoxic hydrolysis (NO2-)
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_O2"): 0,
            ("R3", "Liq", "S_F"): 1 - self.f_SI,
            ("R3", "Liq", "S_A"): 0,
            ("R3", "Liq", "S_I"): self.f_SI,
            ("R3", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R3", "Liq", "S_NH2OH"): 0,
            ("R3", "Liq", "S_N2O"): 0,
            ("R3", "Liq", "S_NO"): 0,
            ("R3", "Liq", "S_NO2"): 0,
            ("R3", "Liq", "S_N2"): 0,
            ("R3", "Liq", "S_NO3"): 0,
            ("R3", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R3", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_S"): -1,
            ("R3", "Liq", "X_H"): 0,
            ("R3", "Liq", "X_PAO"): 0,
            ("R3", "Liq", "X_PP"): 0,
            ("R3", "Liq", "X_PHA"): 0,
            ("R3", "Liq", "X_AOB"): 0,
            ("R3", "Liq", "X_NOB"): 0,
            ("R3", "Liq", "S_K"): 0,
            ("R3", "Liq", "S_Mg"): 0,
            # R4: Anaerobic hydrolysis
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_O2"): 0,
            ("R4", "Liq", "S_F"): 1 - self.f_SI,
            ("R4", "Liq", "S_A"): 0,
            ("R4", "Liq", "S_I"): self.f_SI,
            ("R4", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R4", "Liq", "S_NH2OH"): 0,
            ("R4", "Liq", "S_N2O"): 0,
            ("R4", "Liq", "S_NO"): 0,
            ("R4", "Liq", "S_NO2"): 0,
            ("R4", "Liq", "S_N2"): 0,
            ("R4", "Liq", "S_NO3"): 0,
            ("R4", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R4", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_S"): -1,
            ("R4", "Liq", "X_H"): 0,
            ("R4", "Liq", "X_PAO"): 0,
            ("R4", "Liq", "X_PP"): 0,
            ("R4", "Liq", "X_PHA"): 0,
            ("R4", "Liq", "X_AOB"): 0,
            ("R4", "Liq", "X_NOB"): 0,
            ("R4", "Liq", "S_K"): 0,
            ("R4", "Liq", "S_Mg"): 0,
            # R5: Aerobic growth on S_F
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R5", "Liq", "S_F"): -1 / self.Y_H,
            ("R5", "Liq", "S_A"): 0,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H),
            ("R5", "Liq", "S_NH2OH"): 0,
            ("R5", "Liq", "S_N2O"): 0,
            ("R5", "Liq", "S_NO"): 0,
            ("R5", "Liq", "S_NO2"): 0,
            ("R5", "Liq", "S_N2"): 0,
            ("R5", "Liq", "S_NO3"): 0,
            ("R5", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H),
            ("R5", "Liq", "S_IC"): -(self.i_CBM - self.i_CSF / self.Y_H),
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_S"): 0,
            ("R5", "Liq", "X_H"): 1,
            ("R5", "Liq", "X_PAO"): 0,
            ("R5", "Liq", "X_PP"): 0,
            ("R5", "Liq", "X_PHA"): 0,
            ("R5", "Liq", "X_AOB"): 0,
            ("R5", "Liq", "X_NOB"): 0,
            ("R5", "Liq", "S_K"): 0,
            ("R5", "Liq", "S_Mg"): 0,
            # R6: Aerobic growth on S_A
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R6", "Liq", "S_F"): 0,
            ("R6", "Liq", "S_A"): -1 / self.Y_H,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "S_NH4"): -self.i_NBM,
            ("R6", "Liq", "S_NH2OH"): 0,
            ("R6", "Liq", "S_N2O"): 0,
            ("R6", "Liq", "S_NO"): 0,
            ("R6", "Liq", "S_NO2"): 0,
            ("R6", "Liq", "S_N2"): 0,
            ("R6", "Liq", "S_NO3"): 0,
            ("R6", "Liq", "S_PO4"): -self.i_PBM,
            ("R6", "Liq", "S_IC"): -(self.i_CBM - self.i_CSA / self.Y_H),
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_H"): 1,
            ("R6", "Liq", "X_PAO"): 0,
            ("R6", "Liq", "X_PP"): 0,
            ("R6", "Liq", "X_PHA"): 0,
            ("R6", "Liq", "X_AOB"): 0,
            ("R6", "Liq", "X_NOB"): 0,
            ("R6", "Liq", "S_K"): 0,
            ("R6", "Liq", "S_Mg"): 0,
            # R7: Anoxic growth on S_F (NO3- to NO2-)
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_O2"): 0,
            ("R7", "Liq", "S_F"): -1 / self.Y_H * self.nG,
            ("R7", "Liq", "S_A"): 0,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H * self.nG),
            ("R7", "Liq", "S_NH2OH"): 0,
            ("R7", "Liq", "S_N2O"): 0,
            ("R7", "Liq", "S_NO"): 0,
            ("R7", "Liq", "S_NO2"): (1 - self.Y_H * self.nG)
            / (1.142857 * self.Y_H * self.nG),
            ("R7", "Liq", "S_N2"): 0,
            ("R7", "Liq", "S_NO3"): -(1 - self.Y_H * self.nG)
            / (1.142857 * self.Y_H * self.nG),
            ("R7", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H * self.nG),
            ("R7", "Liq", "S_IC"): -(self.i_CBM - self.i_CSF / self.Y_H * self.nG),
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "X_S"): 0,
            ("R7", "Liq", "X_H"): 1,
            ("R7", "Liq", "X_PAO"): 0,
            ("R7", "Liq", "X_PP"): 0,
            ("R7", "Liq", "X_PHA"): 0,
            ("R7", "Liq", "X_AOB"): 0,
            ("R7", "Liq", "X_NOB"): 0,
            ("R7", "Liq", "S_K"): 0,
            ("R7", "Liq", "S_Mg"): 0,
            # R8: Anoxic growth on S_F (NO2- to NO)
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_O2"): 0,
            ("R8", "Liq", "S_F"): -1 / self.Y_H * self.nG,
            ("R8", "Liq", "S_A"): 0,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H * self.nG),
            ("R8", "Liq", "S_NH2OH"): 0,
            ("R8", "Liq", "S_N2O"): 0,
            ("R8", "Liq", "S_NO"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R8", "Liq", "S_NO2"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R8", "Liq", "S_N2"): 0,
            ("R8", "Liq", "S_NO3"): 0,
            ("R8", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H * self.nG),
            ("R8", "Liq", "S_IC"): -(self.i_CBM - self.i_CSF / self.Y_H * self.nG),
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_S"): 0,
            ("R8", "Liq", "X_H"): 1,
            ("R8", "Liq", "X_PAO"): 0,
            ("R8", "Liq", "X_PP"): 0,
            ("R8", "Liq", "X_PHA"): 0,
            ("R8", "Liq", "X_AOB"): 0,
            ("R8", "Liq", "X_NOB"): 0,
            ("R8", "Liq", "S_K"): 0,
            ("R8", "Liq", "S_Mg"): 0,
            # R9: Anoxic growth on S_F (NO to N2O)
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_O2"): 0,
            ("R9", "Liq", "S_F"): -1 / self.Y_H * self.nG,
            ("R9", "Liq", "S_A"): 0,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H * self.nG),
            ("R9", "Liq", "S_NH2OH"): 0,
            ("R9", "Liq", "S_N2O"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R9", "Liq", "S_NO"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R9", "Liq", "S_NO2"): 0,
            ("R9", "Liq", "S_N2"): 0,
            ("R9", "Liq", "S_NO3"): 0,
            ("R9", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H * self.nG),
            ("R9", "Liq", "S_IC"): -(self.i_CBM - self.i_CSF / self.Y_H * self.nG),
            ("R9", "Liq", "X_I"): 0,
            ("R9", "Liq", "X_S"): 0,
            ("R9", "Liq", "X_H"): 1,
            ("R9", "Liq", "X_PAO"): 0,
            ("R9", "Liq", "X_PP"): 0,
            ("R9", "Liq", "X_PHA"): 0,
            ("R9", "Liq", "X_AOB"): 0,
            ("R9", "Liq", "X_NOB"): 0,
            ("R9", "Liq", "S_K"): 0,
            ("R9", "Liq", "S_Mg"): 0,
            # R10: Anoxic growth on S_F (N2O to N2)
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_O2"): 0,
            ("R10", "Liq", "S_F"): -1 / self.Y_H * self.nG,
            ("R10", "Liq", "S_A"): 0,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H * self.nG),
            ("R10", "Liq", "S_NH2OH"): 0,
            ("R10", "Liq", "S_N2O"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R10", "Liq", "S_NO"): 0,
            ("R10", "Liq", "S_NO2"): 0,
            ("R10", "Liq", "S_N2"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R10", "Liq", "S_NO3"): 0,
            ("R10", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H * self.nG),
            ("R10", "Liq", "S_IC"): -(self.i_CBM - self.i_CSF / self.Y_H * self.nG),
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "X_S"): 0,
            ("R10", "Liq", "X_H"): 1,
            ("R10", "Liq", "X_PAO"): 0,
            ("R10", "Liq", "X_PP"): 0,
            ("R10", "Liq", "X_PHA"): 0,
            ("R10", "Liq", "X_AOB"): 0,
            ("R10", "Liq", "X_NOB"): 0,
            ("R10", "Liq", "S_K"): 0,
            ("R10", "Liq", "S_Mg"): 0,
            # R11: Anoxic growth on S_A (NO3- to NO2-)
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_O2"): 0,
            ("R11", "Liq", "S_F"): 0,
            ("R11", "Liq", "S_A"): -1 / self.Y_H * self.nG,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "S_NH4"): -self.i_NBM,
            ("R11", "Liq", "S_NH2OH"): 0,
            ("R11", "Liq", "S_N2O"): 0,
            ("R11", "Liq", "S_NO"): 0,
            ("R11", "Liq", "S_NO2"): (1 - self.Y_H * self.nG)
            / (1.142857 * self.Y_H * self.nG),
            ("R11", "Liq", "S_N2"): 0,
            ("R11", "Liq", "S_NO3"): -(1 - self.Y_H * self.nG)
            / (1.142857 * self.Y_H * self.nG),
            ("R11", "Liq", "S_PO4"): -self.i_PBM,
            ("R11", "Liq", "S_IC"): -(self.i_CBM - self.i_CSA / self.Y_H * self.nG),
            ("R11", "Liq", "X_I"): 0,
            ("R11", "Liq", "X_S"): 0,
            ("R11", "Liq", "X_H"): 1,
            ("R11", "Liq", "X_PAO"): 0,
            ("R11", "Liq", "X_PP"): 0,
            ("R11", "Liq", "X_PHA"): 0,
            ("R11", "Liq", "X_AOB"): 0,
            ("R11", "Liq", "X_NOB"): 0,
            ("R11", "Liq", "S_K"): 0,
            ("R11", "Liq", "S_Mg"): 0,
            # R12: Anoxic growth on S_A (NO2- to NO)
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_O2"): 0,
            ("R12", "Liq", "S_F"): 0,
            ("R12", "Liq", "S_A"): -1 / self.Y_H * self.nG,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "S_NH4"): -self.i_NBM,
            ("R12", "Liq", "S_NH2OH"): 0,
            ("R12", "Liq", "S_N2O"): 0,
            ("R12", "Liq", "S_NO"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R12", "Liq", "S_NO2"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R12", "Liq", "S_N2"): 0,
            ("R12", "Liq", "S_NO3"): 0,
            ("R12", "Liq", "S_PO4"): -self.i_PBM,
            ("R12", "Liq", "S_IC"): -(self.i_CBM - self.i_CSA / self.Y_H * self.nG),
            ("R12", "Liq", "X_I"): 0,
            ("R12", "Liq", "X_S"): 0,
            ("R12", "Liq", "X_H"): 1,
            ("R12", "Liq", "X_PAO"): 0,
            ("R12", "Liq", "X_PP"): 0,
            ("R12", "Liq", "X_PHA"): 0,
            ("R12", "Liq", "X_AOB"): 0,
            ("R12", "Liq", "X_NOB"): 0,
            ("R12", "Liq", "S_K"): 0,
            ("R12", "Liq", "S_Mg"): 0,
            # R13: Anoxic growth on S_A (NO to N2O)
            ("R13", "Liq", "H2O"): 0,
            ("R13", "Liq", "S_O2"): 0,
            ("R13", "Liq", "S_F"): 0,
            ("R13", "Liq", "S_A"): -1 / self.Y_H * self.nG,
            ("R13", "Liq", "S_I"): 0,
            ("R13", "Liq", "S_NH4"): -self.i_NBM,
            ("R13", "Liq", "S_NH2OH"): 0,
            ("R13", "Liq", "S_N2O"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R13", "Liq", "S_NO"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R13", "Liq", "S_NO2"): 0,
            ("R13", "Liq", "S_N2"): 0,
            ("R13", "Liq", "S_NO3"): 0,
            ("R13", "Liq", "S_PO4"): -self.i_PBM,
            ("R13", "Liq", "S_IC"): -(self.i_CBM - self.i_CSA / self.Y_H * self.nG),
            ("R13", "Liq", "X_I"): 0,
            ("R13", "Liq", "X_S"): 0,
            ("R13", "Liq", "X_H"): 1,
            ("R13", "Liq", "X_PAO"): 0,
            ("R13", "Liq", "X_PP"): 0,
            ("R13", "Liq", "X_PHA"): 0,
            ("R13", "Liq", "X_AOB"): 0,
            ("R13", "Liq", "X_NOB"): 0,
            ("R13", "Liq", "S_K"): 0,
            ("R13", "Liq", "S_Mg"): 0,
            # R14: Anoxic growth on S_A (N2O to N2)
            ("R14", "Liq", "H2O"): 0,
            ("R14", "Liq", "S_O2"): 0,
            ("R14", "Liq", "S_F"): 0,
            ("R14", "Liq", "S_A"): -1 / self.Y_H * self.nG,
            ("R14", "Liq", "S_I"): 0,
            ("R14", "Liq", "S_NH4"): -self.i_NBM,
            ("R14", "Liq", "S_NH2OH"): 0,
            ("R14", "Liq", "S_N2O"): -(1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R14", "Liq", "S_NO"): 0,
            ("R14", "Liq", "S_NO2"): 0,
            ("R14", "Liq", "S_N2"): (1 - self.Y_H * self.nG)
            / (0.5714286 * self.Y_H * self.nG),
            ("R14", "Liq", "S_NO3"): 0,
            ("R14", "Liq", "S_PO4"): -self.i_PBM,
            ("R14", "Liq", "S_IC"): -(self.i_CBM - self.i_CSA / self.Y_H * self.nG),
            ("R14", "Liq", "X_I"): 0,
            ("R14", "Liq", "X_S"): 0,
            ("R14", "Liq", "X_H"): 1,
            ("R14", "Liq", "X_PAO"): 0,
            ("R14", "Liq", "X_PP"): 0,
            ("R14", "Liq", "X_PHA"): 0,
            ("R14", "Liq", "X_AOB"): 0,
            ("R14", "Liq", "X_NOB"): 0,
            ("R14", "Liq", "S_K"): 0,
            ("R14", "Liq", "S_Mg"): 0,
            # R15: Fermentation
            ("R15", "Liq", "H2O"): 0,
            ("R15", "Liq", "S_O2"): 0,
            ("R15", "Liq", "S_F"): -1,
            ("R15", "Liq", "S_A"): 1,
            ("R15", "Liq", "S_I"): 0,
            ("R15", "Liq", "S_NH4"): self.i_NSF,
            ("R15", "Liq", "S_NH2OH"): 0,
            ("R15", "Liq", "S_N2O"): 0,
            ("R15", "Liq", "S_NO"): 0,
            ("R15", "Liq", "S_NO2"): 0,
            ("R15", "Liq", "S_N2"): 0,
            ("R15", "Liq", "S_NO3"): 0,
            ("R15", "Liq", "S_PO4"): self.i_PSF,
            ("R15", "Liq", "S_IC"): -(self.i_CSA - self.i_CSF),
            ("R15", "Liq", "X_I"): 0,
            ("R15", "Liq", "X_S"): 0,
            ("R15", "Liq", "X_H"): 0,
            ("R15", "Liq", "X_PAO"): 0,
            ("R15", "Liq", "X_PP"): 0,
            ("R15", "Liq", "X_PHA"): 0,
            ("R15", "Liq", "X_AUT"): 0,
            ("R15", "Liq", "S_K"): 0,
            ("R15", "Liq", "S_Mg"): 0,
            # R16: Lysis of X_H
            ("R16", "Liq", "H2O"): 0,
            ("R16", "Liq", "S_O2"): 0,
            ("R16", "Liq", "S_F"): 0,
            ("R16", "Liq", "S_A"): 0,
            ("R16", "Liq", "S_I"): 0,
            ("R16", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R16", "Liq", "S_NH2OH"): 0,
            ("R16", "Liq", "S_N2O"): 0,
            ("R16", "Liq", "S_NO"): 0,
            ("R16", "Liq", "S_NO2"): 0,
            ("R16", "Liq", "S_N2"): 0,
            ("R16", "Liq", "S_NO3"): 0,
            ("R16", "Liq", "S_PO4"): -(
                self.f_XI * self.i_PXI + (1 - self.f_XI) * self.i_PXS - self.i_PBM
            ),
            ("R16", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CBM
            ),
            ("R16", "Liq", "X_I"): self.f_XI,
            ("R16", "Liq", "X_S"): 1 - self.f_XI,
            ("R16", "Liq", "X_H"): -1,
            ("R16", "Liq", "X_PAO"): 0,
            ("R16", "Liq", "X_PP"): 0,
            ("R16", "Liq", "X_PHA"): 0,
            ("R16", "Liq", "X_AOB"): 0,
            ("R16", "Liq", "X_NOB"): 0,
            ("R16", "Liq", "S_K"): 0,
            ("R16", "Liq", "S_Mg"): 0,
            # R17: Storage of X_PHA
            ("R17", "Liq", "H2O"): 0,
            ("R17", "Liq", "S_O2"): 0,
            ("R17", "Liq", "S_F"): 0,
            ("R17", "Liq", "S_A"): -1,
            ("R17", "Liq", "S_I"): 0,
            ("R17", "Liq", "S_NH4"): 0,
            ("R17", "Liq", "S_NH2OH"): 0,
            ("R17", "Liq", "S_N2O"): 0,
            ("R17", "Liq", "S_NO"): 0,
            ("R17", "Liq", "S_NO2"): 0,
            ("R17", "Liq", "S_N2"): 0,
            ("R17", "Liq", "S_NO3"): 0,
            ("R17", "Liq", "S_PO4"): self.Y_PO4,
            ("R17", "Liq", "S_IC"): self.i_CSA - self.i_CXPHA,
            ("R17", "Liq", "X_I"): 0,
            ("R17", "Liq", "X_S"): 0,
            ("R17", "Liq", "X_H"): 0,
            ("R17", "Liq", "X_PAO"): 0,
            ("R17", "Liq", "X_PP"): -self.Y_PO4,
            ("R17", "Liq", "X_PHA"): 1,
            ("R17", "Liq", "X_AOB"): 0,
            ("R17", "Liq", "X_NOB"): 0,
            ("R17", "Liq", "S_K"): self.Y_PO4 * self.i_KXPP,
            ("R17", "Liq", "S_Mg"): self.Y_PO4 * self.i_MgXPP,
            # R18: Aerobic Storage of X_PP
            ("R18", "Liq", "H2O"): 0,
            ("R18", "Liq", "S_O2"): -self.Y_PHA,
            ("R18", "Liq", "S_F"): 0,
            ("R18", "Liq", "S_A"): 0,
            ("R18", "Liq", "S_I"): 0,
            ("R18", "Liq", "S_NH4"): 0,
            ("R18", "Liq", "S_NH2OH"): 0,
            ("R18", "Liq", "S_N2O"): 0,
            ("R18", "Liq", "S_NO"): 0,
            ("R18", "Liq", "S_NO2"): 0,
            ("R18", "Liq", "S_N2"): 0,
            ("R18", "Liq", "S_NO3"): 0,
            ("R18", "Liq", "S_PO4"): -1,
            ("R18", "Liq", "S_IC"): self.Y_PHA * self.i_CXPHA,
            ("R18", "Liq", "X_I"): 0,
            ("R18", "Liq", "X_S"): 0,
            ("R18", "Liq", "X_H"): 0,
            ("R18", "Liq", "X_PAO"): 0,
            ("R18", "Liq", "X_PP"): 1,
            ("R18", "Liq", "X_PHA"): -self.Y_PHA,
            ("R18", "Liq", "X_AOB"): 0,
            ("R18", "Liq", "X_NOB"): 0,
            ("R18", "Liq", "S_K"): -self.i_KXPP,
            ("R18", "Liq", "S_Mg"): -self.i_MgXPP,
            # R19: Anoxic Storage of X_PP (NO3- to NO2-)
            ("R19", "Liq", "H2O"): 0,
            ("R19", "Liq", "S_O2"): 0,
            ("R19", "Liq", "S_F"): 0,
            ("R19", "Liq", "S_A"): 0,
            ("R19", "Liq", "S_I"): 0,
            ("R19", "Liq", "S_NH4"): 0,
            ("R19", "Liq", "S_NH2OH"): 0,
            ("R19", "Liq", "S_N2O"): 0,
            ("R19", "Liq", "S_NO"): 0,
            ("R19", "Liq", "S_NO2"): self.Y_PHA * self.nG / 1.142857,
            ("R19", "Liq", "S_N2"): 0,
            ("R19", "Liq", "S_NO3"): -self.Y_PHA * self.nG / 1.142857,
            ("R19", "Liq", "S_PO4"): -1,
            ("R19", "Liq", "S_IC"): self.Y_PHA * self.i_CXPHA * self.nG,
            ("R19", "Liq", "X_I"): 0,
            ("R19", "Liq", "X_S"): 0,
            ("R19", "Liq", "X_H"): 0,
            ("R19", "Liq", "X_PAO"): 0,
            ("R19", "Liq", "X_PP"): 1,
            ("R19", "Liq", "X_PHA"): -self.Y_PHA * self.nG,
            ("R19", "Liq", "X_AOB"): 0,
            ("R19", "Liq", "X_NOB"): 0,
            ("R19", "Liq", "S_K"): -self.i_KXPP,
            ("R19", "Liq", "S_Mg"): -self.i_MgXPP,
            # R20: Anoxic Storage of X_PP (NO2- to NO2)
            ("R20", "Liq", "H2O"): 0,
            ("R20", "Liq", "S_O2"): 0,
            ("R20", "Liq", "S_F"): 0,
            ("R20", "Liq", "S_A"): 0,
            ("R20", "Liq", "S_I"): 0,
            ("R20", "Liq", "S_NH4"): 0,
            ("R20", "Liq", "S_NH2OH"): 0,
            ("R20", "Liq", "S_N2O"): 0,
            ("R20", "Liq", "S_NO"): self.Y_PHA * self.nG / 0.5714286,
            ("R20", "Liq", "S_NO2"): -self.Y_PHA * self.nG / 0.5714286,
            ("R20", "Liq", "S_N2"): 0,
            ("R20", "Liq", "S_NO3"): 0,
            ("R20", "Liq", "S_PO4"): -1,
            ("R20", "Liq", "S_IC"): self.Y_PHA * self.i_CXPHA * self.nG,
            ("R20", "Liq", "X_I"): 0,
            ("R20", "Liq", "X_S"): 0,
            ("R20", "Liq", "X_H"): 0,
            ("R20", "Liq", "X_PAO"): 0,
            ("R20", "Liq", "X_PP"): 1,
            ("R20", "Liq", "X_PHA"): -self.Y_PHA * self.nG,
            ("R20", "Liq", "X_AOB"): 0,
            ("R20", "Liq", "X_NOB"): 0,
            ("R20", "Liq", "S_K"): -self.i_KXPP,
            ("R20", "Liq", "S_Mg"): -self.i_MgXPP,
            # R21: Anoxic Storage of X_PP (NO to N2O)
            ("R21", "Liq", "H2O"): 0,
            ("R21", "Liq", "S_O2"): 0,
            ("R21", "Liq", "S_F"): 0,
            ("R21", "Liq", "S_A"): 0,
            ("R21", "Liq", "S_I"): 0,
            ("R21", "Liq", "S_NH4"): 0,
            ("R21", "Liq", "S_NH2OH"): 0,
            ("R21", "Liq", "S_N2O"): self.Y_PHA * self.nG / 0.5714286,
            ("R21", "Liq", "S_NO"): -self.Y_PHA * self.nG / 0.5714286,
            ("R21", "Liq", "S_NO2"): 0,
            ("R21", "Liq", "S_N2"): 0,
            ("R21", "Liq", "S_NO3"): 0,
            ("R21", "Liq", "S_PO4"): -1,
            ("R21", "Liq", "S_IC"): self.Y_PHA * self.i_CXPHA * self.nG,
            ("R21", "Liq", "X_I"): 0,
            ("R21", "Liq", "X_S"): 0,
            ("R21", "Liq", "X_H"): 0,
            ("R21", "Liq", "X_PAO"): 0,
            ("R21", "Liq", "X_PP"): 1,
            ("R21", "Liq", "X_PHA"): -self.Y_PHA * self.nG,
            ("R21", "Liq", "X_AOB"): 0,
            ("R21", "Liq", "X_NOB"): 0,
            ("R21", "Liq", "S_K"): -self.i_KXPP,
            ("R21", "Liq", "S_Mg"): -self.i_MgXPP,
            # R22: Anoxic Storage of X_PP (N2O to N2)
            ("R22", "Liq", "H2O"): 0,
            ("R22", "Liq", "S_O2"): 0,
            ("R22", "Liq", "S_F"): 0,
            ("R22", "Liq", "S_A"): 0,
            ("R22", "Liq", "S_I"): 0,
            ("R22", "Liq", "S_NH4"): 0,
            ("R22", "Liq", "S_NH2OH"): 0,
            ("R22", "Liq", "S_N2O"): -self.Y_PHA * self.nG / 0.5714286,
            ("R22", "Liq", "S_NO"): 0,
            ("R22", "Liq", "S_NO2"): 0,
            ("R22", "Liq", "S_N2"): self.Y_PHA * self.nG / 0.5714286,
            ("R22", "Liq", "S_NO3"): 0,
            ("R22", "Liq", "S_PO4"): -1,
            ("R22", "Liq", "S_IC"): self.Y_PHA * self.i_CXPHA * self.nG,
            ("R22", "Liq", "X_I"): 0,
            ("R22", "Liq", "X_S"): 0,
            ("R22", "Liq", "X_H"): 0,
            ("R22", "Liq", "X_PAO"): 0,
            ("R22", "Liq", "X_PP"): 1,
            ("R22", "Liq", "X_PHA"): -self.Y_PHA * self.nG,
            ("R22", "Liq", "X_AOB"): 0,
            ("R22", "Liq", "X_NOB"): 0,
            ("R22", "Liq", "S_K"): -self.i_KXPP,
            ("R22", "Liq", "S_Mg"): -self.i_MgXPP,
            # R23: Aerobic growth of X_PAO
            ("R23", "Liq", "H2O"): 0,
            ("R23", "Liq", "S_O2"): 1 - 1 / self.Y_PAO,
            ("R23", "Liq", "S_F"): 0,
            ("R23", "Liq", "S_A"): 0,
            ("R23", "Liq", "S_I"): 0,
            ("R23", "Liq", "S_NH4"): -self.i_NBM,
            ("R23", "Liq", "S_NH2OH"): 0,
            ("R23", "Liq", "S_N2O"): 0,
            ("R23", "Liq", "S_NO"): 0,
            ("R23", "Liq", "S_NO2"): 0,
            ("R23", "Liq", "S_N2"): 0,
            ("R23", "Liq", "S_NO3"): 0,
            ("R23", "Liq", "S_PO4"): -self.i_PBM,
            ("R23", "Liq", "S_IC"): -self.i_CBM + self.i_CXPHA / self.Y_PAO,
            ("R23", "Liq", "X_I"): 0,
            ("R23", "Liq", "X_S"): 0,
            ("R23", "Liq", "X_H"): 0,
            ("R23", "Liq", "X_PAO"): 1,
            ("R23", "Liq", "X_PP"): 0,
            ("R23", "Liq", "X_PHA"): -1 / self.Y_PAO,
            ("R23", "Liq", "X_AOB"): 0,
            ("R23", "Liq", "X_NOB"): 0,
            ("R23", "Liq", "S_K"): 0,
            ("R23", "Liq", "S_Mg"): 0,
            # R24: Anoxic growth of X_PAO (NO3- to NO2-)
            ("R24", "Liq", "H2O"): 0,
            ("R24", "Liq", "S_O2"): 0,
            ("R24", "Liq", "S_F"): 0,
            ("R24", "Liq", "S_A"): 0,
            ("R24", "Liq", "S_I"): 0,
            ("R24", "Liq", "S_NH4"): -self.i_NBM,
            ("R24", "Liq", "S_NH2OH"): 0,
            ("R24", "Liq", "S_N2O"): 0,
            ("R24", "Liq", "S_NO"): 0,
            ("R24", "Liq", "S_NO2"): (1 - self.Y_PAO * self.nG)
            / (1.142857 * self.Y_PAO * self.nG),
            ("R24", "Liq", "S_N2"): 0,
            ("R24", "Liq", "S_NO3"): -(1 - self.Y_PAO * self.nG)
            / (1.142857 * self.Y_PAO * self.nG),
            ("R24", "Liq", "S_PO4"): -self.i_PBM,
            ("R24", "Liq", "S_IC"): -self.i_CBM + self.i_CXPHA / self.Y_PAO,
            ("R24", "Liq", "X_I"): 0,
            ("R24", "Liq", "X_S"): 0,
            ("R24", "Liq", "X_H"): 0,
            ("R24", "Liq", "X_PAO"): 1,
            ("R24", "Liq", "X_PP"): 0,
            ("R24", "Liq", "X_PHA"): -1 / (self.Y_PAO * self.nG),
            ("R24", "Liq", "X_AOB"): 0,
            ("R24", "Liq", "X_NOB"): 0,
            ("R24", "Liq", "S_K"): 0,
            ("R24", "Liq", "S_Mg"): 0,
            # R25: Anoxic growth of X_PAO (NO2- to NO)
            ("R25", "Liq", "H2O"): 0,
            ("R25", "Liq", "S_O2"): 0,
            ("R25", "Liq", "S_F"): 0,
            ("R25", "Liq", "S_A"): 0,
            ("R25", "Liq", "S_I"): 0,
            ("R25", "Liq", "S_NH4"): -self.i_NBM,
            ("R25", "Liq", "S_NH2OH"): 0,
            ("R25", "Liq", "S_N2O"): 0,
            ("R25", "Liq", "S_NO"): (1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R25", "Liq", "S_NO2"): -(1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R25", "Liq", "S_N2"): 0,
            ("R25", "Liq", "S_NO3"): 0,
            ("R25", "Liq", "S_PO4"): -self.i_PBM,
            ("R25", "Liq", "S_IC"): -self.i_CBM + self.i_CXPHA / self.Y_PAO,
            ("R25", "Liq", "X_I"): 0,
            ("R25", "Liq", "X_S"): 0,
            ("R25", "Liq", "X_H"): 0,
            ("R25", "Liq", "X_PAO"): 1,
            ("R25", "Liq", "X_PP"): 0,
            ("R25", "Liq", "X_PHA"): -1 / (self.Y_PAO * self.nG),
            ("R25", "Liq", "X_AOB"): 0,
            ("R25", "Liq", "X_NOB"): 0,
            ("R25", "Liq", "S_K"): 0,
            ("R25", "Liq", "S_Mg"): 0,
            # R26: Anoxic growth of X_PAO (NO to N2O)
            ("R26", "Liq", "H2O"): 0,
            ("R26", "Liq", "S_O2"): 0,
            ("R26", "Liq", "S_F"): 0,
            ("R26", "Liq", "S_A"): 0,
            ("R26", "Liq", "S_I"): 0,
            ("R26", "Liq", "S_NH4"): -self.i_NBM,
            ("R26", "Liq", "S_NH2OH"): 0,
            ("R26", "Liq", "S_N2O"): (1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R26", "Liq", "S_NO"): -(1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R26", "Liq", "S_NO2"): 0,
            ("R26", "Liq", "S_N2"): 0,
            ("R26", "Liq", "S_NO3"): 0,
            ("R26", "Liq", "S_PO4"): -self.i_PBM,
            ("R26", "Liq", "S_IC"): -self.i_CBM + self.i_CXPHA / self.Y_PAO,
            ("R26", "Liq", "X_I"): 0,
            ("R26", "Liq", "X_S"): 0,
            ("R26", "Liq", "X_H"): 0,
            ("R26", "Liq", "X_PAO"): 1,
            ("R26", "Liq", "X_PP"): 0,
            ("R26", "Liq", "X_PHA"): -1 / (self.Y_PAO * self.nG),
            ("R26", "Liq", "X_AOB"): 0,
            ("R26", "Liq", "X_NOB"): 0,
            ("R26", "Liq", "S_K"): 0,
            ("R26", "Liq", "S_Mg"): 0,
            # R27: Anoxic growth of X_PAO (N2O to N2)
            ("R27", "Liq", "H2O"): 0,
            ("R27", "Liq", "S_O2"): 0,
            ("R27", "Liq", "S_F"): 0,
            ("R27", "Liq", "S_A"): 0,
            ("R27", "Liq", "S_I"): 0,
            ("R27", "Liq", "S_NH4"): -self.i_NBM,
            ("R27", "Liq", "S_NH2OH"): 0,
            ("R27", "Liq", "S_N2O"): -(1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R27", "Liq", "S_NO"): 0,
            ("R27", "Liq", "S_NO2"): 0,
            ("R27", "Liq", "S_N2"): (1 - self.Y_PAO * self.nG)
            / (0.5714286 * self.Y_PAO * self.nG),
            ("R27", "Liq", "S_NO3"): 0,
            ("R27", "Liq", "S_PO4"): -self.i_PBM,
            ("R27", "Liq", "S_IC"): self.i_CBM - self.i_CXPHA / self.Y_PAO,
            ("R27", "Liq", "X_I"): 0,
            ("R27", "Liq", "X_S"): 0,
            ("R27", "Liq", "X_H"): 0,
            ("R27", "Liq", "X_PAO"): 1,
            ("R27", "Liq", "X_PP"): 0,
            ("R27", "Liq", "X_PHA"): -1 / (self.Y_PAO * self.nG),
            ("R27", "Liq", "X_AOB"): 0,
            ("R27", "Liq", "X_NOB"): 0,
            ("R27", "Liq", "S_K"): 0,
            ("R27", "Liq", "S_Mg"): 0,
            # R28: Lysis od X_PAO
            ("R28", "Liq", "H2O"): 0,
            ("R28", "Liq", "S_O2"): 0,
            ("R28", "Liq", "S_F"): 0,
            ("R28", "Liq", "S_A"): 0,
            ("R28", "Liq", "S_I"): 0,
            ("R28", "Liq", "S_NH4"): -(
                self.f_XI * self.I_NXI + (1 - self.f_XI) * self.I_NXS - self.I_NBM
            ),
            ("R28", "Liq", "S_NH2OH"): 0,
            ("R28", "Liq", "S_N2O"): 0,
            ("R28", "Liq", "S_NO"): 0,
            ("R28", "Liq", "S_NO2"): 0,
            ("R28", "Liq", "S_N2"): 0,
            ("R28", "Liq", "S_NO3"): 0,
            ("R28", "Liq", "S_PO4"): -(
                self.f_XI * self.I_PXI + (1 - self.f_XI) * self.I_PXS - self.I_PBM
            ),
            ("R28", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CBM
            ),
            ("R28", "Liq", "X_I"): self.f_XI,
            ("R28", "Liq", "X_S"): 1 - self.f_XI,
            ("R28", "Liq", "X_H"): 0,
            ("R28", "Liq", "X_PAO"): -1,
            ("R28", "Liq", "X_PP"): 0,
            ("R28", "Liq", "X_PHA"): 0,
            ("R28", "Liq", "X_AOB"): 0,
            ("R28", "Liq", "X_NOB"): 0,
            ("R28", "Liq", "S_K"): 0,
            ("R28", "Liq", "S_Mg"): 0,
            # R29: Lysis of X_PP
            ("R29", "Liq", "H2O"): 0,
            ("R29", "Liq", "S_O2"): 0,
            ("R29", "Liq", "S_F"): 0,
            ("R29", "Liq", "S_A"): 0,
            ("R29", "Liq", "S_I"): 0,
            ("R29", "Liq", "S_NH4"): 0,
            ("R29", "Liq", "S_NH2OH"): 0,
            ("R29", "Liq", "S_N2O"): 0,
            ("R29", "Liq", "S_NO"): 0,
            ("R29", "Liq", "S_NO2"): 0,
            ("R29", "Liq", "S_N2"): 0,
            ("R29", "Liq", "S_NO3"): 0,
            ("R29", "Liq", "S_PO4"): 1,
            ("R29", "Liq", "S_IC"): 0,
            ("R29", "Liq", "X_I"): 0,
            ("R29", "Liq", "X_S"): 0,
            ("R29", "Liq", "X_H"): 0,
            ("R29", "Liq", "X_PAO"): 0,
            ("R29", "Liq", "X_PP"): -1,
            ("R29", "Liq", "X_PHA"): 0,
            ("R29", "Liq", "X_AOB"): 0,
            ("R29", "Liq", "X_NOB"): 0,
            ("R29", "Liq", "S_K"): self.i_KXPP,
            ("R29", "Liq", "S_Mg"): self.i_MgXPP,
            # R30: Lysis of X_PHA
            ("R30", "Liq", "H2O"): 0,
            ("R30", "Liq", "S_O2"): 0,
            ("R30", "Liq", "S_F"): 0,
            ("R30", "Liq", "S_A"): 1,
            ("R30", "Liq", "S_I"): 0,
            ("R30", "Liq", "S_NH4"): 0,
            ("R30", "Liq", "S_NH2OH"): 0,
            ("R30", "Liq", "S_N2O"): 0,
            ("R30", "Liq", "S_NO"): 0,
            ("R30", "Liq", "S_NO2"): 0,
            ("R30", "Liq", "S_N2"): 0,
            ("R30", "Liq", "S_NO3"): 0,
            ("R30", "Liq", "S_PO4"): 0,
            ("R30", "Liq", "S_IC"): self.i_CXPHA - self.i_CSA,
            ("R30", "Liq", "X_I"): 0,
            ("R30", "Liq", "X_S"): 0,
            ("R30", "Liq", "X_H"): 0,
            ("R30", "Liq", "X_PAO"): 0,
            ("R30", "Liq", "X_PP"): 0,
            ("R30", "Liq", "X_PHA"): -1,
            ("R30", "Liq", "X_AOB"): 0,
            ("R30", "Liq", "X_NOB"): 0,
            ("R30", "Liq", "S_K"): 0,
            ("R30", "Liq", "S_Mg"): 0,
            # R31: NH3 oxidation to NH2OH with O2 consumption
            ("R31", "Liq", "H2O"): 0,
            ("R31", "Liq", "S_O2"): -1.42857,
            ("R31", "Liq", "S_F"): 0,
            ("R31", "Liq", "S_A"): 0,
            ("R31", "Liq", "S_I"): 0,
            ("R31", "Liq", "S_NH4"): -1,
            ("R31", "Liq", "S_NH2OH"): 1,
            ("R31", "Liq", "S_N2O"): 0,
            ("R31", "Liq", "S_NO"): 0,
            ("R31", "Liq", "S_NO2"): 0,
            ("R31", "Liq", "S_N2"): 0,
            ("R31", "Liq", "S_NO3"): 0,
            ("R31", "Liq", "S_PO4"): 0,
            ("R31", "Liq", "S_IC"): 0,
            ("R31", "Liq", "X_I"): 0,
            ("R31", "Liq", "X_S"): 0,
            ("R31", "Liq", "X_H"): 0,
            ("R31", "Liq", "X_PAO"): 0,
            ("R31", "Liq", "X_PP"): 0,
            ("R31", "Liq", "X_PHA"): 0,
            ("R31", "Liq", "X_AOB"): 0,
            ("R31", "Liq", "X_NOB"): 0,
            ("R31", "Liq", "S_K"): 0,
            ("R31", "Liq", "S_Mg"): 0,
            # R32: NH2OH to NO coupled with O2 reduction (X_AOB growth)
            ("R32", "Liq", "H2O"): 0,
            ("R32", "Liq", "S_O2"): -(1.714286 - self.Y_AOB) / self.Y_AOB,
            ("R32", "Liq", "S_F"): 0,
            ("R32", "Liq", "S_A"): 0,
            ("R32", "Liq", "S_I"): 0,
            ("R32", "Liq", "S_NH4"): -self.i_NBM,
            ("R32", "Liq", "S_NH2OH"): -1 / self.Y_AOB,
            ("R32", "Liq", "S_N2O"): 0,
            ("R32", "Liq", "S_NO"): 1 / self.Y_AOB,
            ("R32", "Liq", "S_NO2"): 0,
            ("R32", "Liq", "S_N2"): 0,
            ("R32", "Liq", "S_NO3"): 0,
            ("R32", "Liq", "S_PO4"): -self.i_PBM,
            ("R32", "Liq", "S_IC"): -self.i_CBM,
            ("R32", "Liq", "X_I"): 0,
            ("R32", "Liq", "X_S"): 0,
            ("R32", "Liq", "X_H"): 0,
            ("R32", "Liq", "X_PAO"): 0,
            ("R32", "Liq", "X_PP"): 0,
            ("R32", "Liq", "X_PHA"): 0,
            ("R32", "Liq", "X_AOB"): 1,
            ("R32", "Liq", "X_NOB"): 0,
            ("R32", "Liq", "S_K"): 0,
            ("R32", "Liq", "S_Mg"): 0,
            # R33: NO oxidation to NO2- coupled with O2 reduction
            ("R33", "Liq", "H2O"): 0,
            ("R33", "Liq", "S_O2"): -0.5714286,
            ("R33", "Liq", "S_F"): 0,
            ("R33", "Liq", "S_A"): 0,
            ("R33", "Liq", "S_I"): 0,
            ("R33", "Liq", "S_NH4"): 0,
            ("R33", "Liq", "S_NH2OH"): 0,
            ("R33", "Liq", "S_N2O"): 0,
            ("R33", "Liq", "S_NO"): -1,
            ("R33", "Liq", "S_NO2"): 1,
            ("R33", "Liq", "S_N2"): 0,
            ("R33", "Liq", "S_NO3"): 0,
            ("R33", "Liq", "S_PO4"): 0,
            ("R33", "Liq", "S_IC"): 0,
            ("R33", "Liq", "X_I"): 0,
            ("R33", "Liq", "X_S"): 0,
            ("R33", "Liq", "X_H"): 0,
            ("R33", "Liq", "X_PAO"): 0,
            ("R33", "Liq", "X_PP"): 0,
            ("R33", "Liq", "X_PHA"): 0,
            ("R33", "Liq", "X_AOB"): 0,
            ("R33", "Liq", "X_NOB"): 0,
            ("R33", "Liq", "S_K"): 0,
            ("R33", "Liq", "S_Mg"): 0,
            # R34: NO to N2O coupled with NH2OH to NO2- (N2O from NN pathway)
            ("R34", "Liq", "H2O"): 0,
            ("R34", "Liq", "S_O2"): 0,
            ("R34", "Liq", "S_F"): 0,
            ("R34", "Liq", "S_A"): 0,
            ("R34", "Liq", "S_I"): 0,
            ("R34", "Liq", "S_NH4"): 0,
            ("R34", "Liq", "S_NH2OH"): -1,
            ("R34", "Liq", "S_N2O"): 4,
            ("R34", "Liq", "S_NO"): -4,
            ("R34", "Liq", "S_NO2"): 1,
            ("R34", "Liq", "S_N2"): 0,
            ("R34", "Liq", "S_NO3"): 0,
            ("R34", "Liq", "S_PO4"): 0,
            ("R34", "Liq", "S_IC"): 0,
            ("R34", "Liq", "X_I"): 0,
            ("R34", "Liq", "X_S"): 0,
            ("R34", "Liq", "X_H"): 0,
            ("R34", "Liq", "X_PAO"): 0,
            ("R34", "Liq", "X_PP"): 0,
            ("R34", "Liq", "X_PHA"): 0,
            ("R34", "Liq", "X_AOB"): 0,
            ("R34", "Liq", "X_NOB"): 0,
            ("R34", "Liq", "S_K"): 0,
            ("R34", "Liq", "S_Mg"): 0,
            # R35: HNO2 to N2O coupled with NH2OH to NO2- (N2O from ND pathway)
            ("R35", "Liq", "H2O"): 0,
            ("R35", "Liq", "S_O2"): 0,
            ("R35", "Liq", "S_F"): 0,
            ("R35", "Liq", "S_A"): 0,
            ("R35", "Liq", "S_I"): 0,
            ("R35", "Liq", "S_NH4"): 0,
            ("R35", "Liq", "S_NH2OH"): -1,
            ("R35", "Liq", "S_N2O"): 2,
            ("R35", "Liq", "S_NO"): 0,
            ("R35", "Liq", "S_NO2"): -1,
            ("R35", "Liq", "S_N2"): 0,
            ("R35", "Liq", "S_NO3"): 0,
            ("R35", "Liq", "S_PO4"): 0,
            ("R35", "Liq", "S_IC"): 0,
            ("R35", "Liq", "X_I"): 0,
            ("R35", "Liq", "X_S"): 0,
            ("R35", "Liq", "X_H"): 0,
            ("R35", "Liq", "X_PAO"): 0,
            ("R35", "Liq", "X_PP"): 0,
            ("R35", "Liq", "X_PHA"): 0,
            ("R35", "Liq", "X_AOB"): 0,
            ("R35", "Liq", "X_NOB"): 0,
            ("R35", "Liq", "S_K"): 0,
            ("R35", "Liq", "S_Mg"): 0,
            # R36: Aerobic growth of X_NOB
            ("R36", "Liq", "H2O"): 0,
            ("R36", "Liq", "S_O2"): -(1.142857 - self.Y_NOB) / self.Y_NOB,
            ("R36", "Liq", "S_F"): 0,
            ("R36", "Liq", "S_A"): 0,
            ("R36", "Liq", "S_I"): 0,
            ("R36", "Liq", "S_NH4"): -self.i_NBM,
            ("R36", "Liq", "S_NH2OH"): 0,
            ("R36", "Liq", "S_N2O"): 0,
            ("R36", "Liq", "S_NO"): 0,
            ("R36", "Liq", "S_NO2"): -1 / self.Y_NOB,
            ("R36", "Liq", "S_N2"): 0,
            ("R36", "Liq", "S_NO3"): 1 / self.Y_NOB,
            ("R36", "Liq", "S_PO4"): -self.i_PBM,
            ("R36", "Liq", "S_IC"): -self.i_CBM,
            ("R36", "Liq", "X_I"): 0,
            ("R36", "Liq", "X_S"): 0,
            ("R36", "Liq", "X_H"): 0,
            ("R36", "Liq", "X_PAO"): 0,
            ("R36", "Liq", "X_PP"): 0,
            ("R36", "Liq", "X_PHA"): 0,
            ("R36", "Liq", "X_AOB"): 0,
            ("R36", "Liq", "X_NOB"): 1,
            ("R36", "Liq", "S_K"): 0,
            ("R36", "Liq", "S_Mg"): 0,
            # R37: Lysis of X_AOB
            ("R37", "Liq", "H2O"): 0,
            ("R37", "Liq", "S_O2"): 0,
            ("R37", "Liq", "S_F"): 0,
            ("R37", "Liq", "S_A"): 0,
            ("R37", "Liq", "S_I"): 0,
            ("R37", "Liq", "S_NH4"): -(
                self.f_XI * self.I_NXI + (1 - self.f_XI) * self.I_NXS - self.I_NBM
            ),
            ("R37", "Liq", "S_NH2OH"): 0,
            ("R37", "Liq", "S_N2O"): 0,
            ("R37", "Liq", "S_NO"): 0,
            ("R37", "Liq", "S_NO2"): 0,
            ("R37", "Liq", "S_N2"): 0,
            ("R37", "Liq", "S_NO3"): 0,
            ("R37", "Liq", "S_PO4"): -(
                self.f_XI * self.I_PXI + (1 - self.f_XI) * self.I_PXS - self.I_PBM
            ),
            ("R37", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CBM
            ),
            ("R37", "Liq", "X_I"): self.f_XI,
            ("R37", "Liq", "X_S"): 1 - self.f_XI,
            ("R37", "Liq", "X_H"): 0,
            ("R37", "Liq", "X_PAO"): 0,
            ("R37", "Liq", "X_PP"): 0,
            ("R37", "Liq", "X_PHA"): 0,
            ("R37", "Liq", "X_AOB"): -1,
            ("R37", "Liq", "X_NOB"): 0,
            ("R37", "Liq", "S_K"): 0,
            ("R37", "Liq", "S_Mg"): 0,
            # R38: Lysis of X_NOB
            ("R38", "Liq", "H2O"): 0,
            ("R38", "Liq", "S_O2"): 0,
            ("R38", "Liq", "S_F"): 0,
            ("R38", "Liq", "S_A"): 0,
            ("R38", "Liq", "S_I"): 0,
            ("R38", "Liq", "S_NH4"): -(
                self.f_XI * self.I_NXI + (1 - self.f_XI) * self.I_NXS - self.I_NBM
            ),
            ("R38", "Liq", "S_NH2OH"): 0,
            ("R38", "Liq", "S_N2O"): 0,
            ("R38", "Liq", "S_NO"): 0,
            ("R38", "Liq", "S_NO2"): 0,
            ("R38", "Liq", "S_N2"): 0,
            ("R38", "Liq", "S_NO3"): 0,
            ("R38", "Liq", "S_PO4"): -(
                self.f_XI * self.I_PXI + (1 - self.f_XI) * self.I_PXS - self.I_PBM
            ),
            ("R38", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CBM
            ),
            ("R38", "Liq", "X_I"): self.f_XI,
            ("R38", "Liq", "X_S"): 1 - self.f_XI,
            ("R38", "Liq", "X_H"): 0,
            ("R38", "Liq", "X_PAO"): 0,
            ("R38", "Liq", "X_PP"): 0,
            ("R38", "Liq", "X_PHA"): 0,
            ("R38", "Liq", "X_AOB"): 0,
            ("R38", "Liq", "X_NOB"): -1,
            ("R38", "Liq", "S_K"): 0,
            ("R38", "Liq", "S_Mg"): 0,
        }
        # Fix all the variables we just created
        for v in self.component_objects(pyo.Var, descend_into=False):
            v.fix()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "reaction_rate": {"method": "_rxn_rate"},
            }
        )
        obj.add_default_units(
            {
                "time": pyo.units.s,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.kmol,
                "temperature": pyo.units.K,
            }
        )


class ASM2dN2OReactionScaler(CustomScalerBase):
    """
    Scaler for the Activated Sludge Model No.2d reaction package.
    Variables are scaled by their default scaling factor (if no user input provided), and constraints
    are scaled using the inverse maximum scheme.
    """

    # TODO: Revisit this scaling factor
    DEFAULT_SCALING_FACTORS = {"reaction_rate": 1e2}

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):

        if model.is_property_constructed("reaction_rate"):
            for j in model.reaction_rate.values():
                self.scale_variable_by_default(j, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # TODO: Revisit this scaling methodology
        # Consider scale_constraint_by_default or scale_constraints_by_jacobian_norm
        if model.is_property_constructed("rate_expression"):
            for j in model.rate_expression.values():
                self.scale_constraint_by_nominal_value(
                    j,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )


class _ASM2dN2OReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ASM2dN2OReactionScaler

    def initialize(self, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        init_log.info("Initialization Complete.")


@declare_process_block_class(
    "ASM2dN2OReactionBlock", block_class=_ASM2dN2OReactionBlock
)
class ASM2dN2OReactionBlockData(ReactionBlockDataBase):
    """
    Reaction Block for ASM2d-PSFe-GHG.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super().build()

        # Create references to state vars
        # Concentration
        add_object_reference(self, "conc_mass_comp_ref", self.state_ref.conc_mass_comp)

    # Rate of reaction method - see asm2d_PSFe_GHG.c
    def _rxn_rate(self):
        self.reaction_rate = pyo.Var(
            self.params.rate_reaction_idx,
            initialize=0,
            doc="Rate of reaction",
            units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
        )
        self.eps = pyo.Param(
            initialize=1e-6,
            units=pyo.units.kg / pyo.units.m**3,
            mutable=True,
            doc="epsilon",
        )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KL_X")
        def MonodX(b):
            return (
                b.conc_mass_comp_ref["X_S"]
                / b.conc_mass_comp_ref["X_H"]
                / (
                    b.params.KL_X
                    + b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"]
                )
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KP_PHA")
        def MonodPHA(b):
            return (
                b.conc_mass_comp_ref["X_PHA"]
                / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
                / (
                    b.params.KP_PHA
                    + b.conc_mass_comp_ref["X_PHA"]
                    / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
                )
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KI_PP")
        def MonodIPP(b):
            return (
                b.params.K_MAX
                - b.conc_mass_comp_ref["X_PP"]
                / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
            ) / (
                b.params.KI_PP
                + b.params.K_MAX
                - b.conc_mass_comp_ref["X_PP"]
                / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KP_PP")
        def MonodPP(b):
            return (
                b.conc_mass_comp_ref["X_PP"]
                / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
                / (
                    b.params.KP_PP
                    + b.conc_mass_comp_ref["X_PHA"]
                    / (b.conc_mass_comp_ref["X_PAO"] + self.eps)
                )
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KH_PO4")
        def MonodPO4(b):
            return b.conc_mass_comp_ref["S_PO4"] / (
                b.params.KH_PO4 + b.conc_mass_comp_ref["S_PO4"]
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KH_NH4")
        def MonodNH4(b):
            return b.conc_mass_comp_ref["S_NH4"] / (
                b.params.KH_NH4 + b.conc_mass_comp_ref["S_NH4"]
            )

        @self.Expression(doc="Monod equation for KP_P")
        def MonodSPO4_P(b):
            return b.conc_mass_comp_ref["S_PO4"] / (
                b.params.KP_P + b.conc_mass_comp_ref["S_PO4"]
            )

        @self.Expression(self.flowsheet().time, doc="Monod equation for KP_NH4")
        def MonodSNH4_P(b):
            return b.conc_mass_comp_ref["S_NH4"] / (
                b.params.KP_NH4 + b.conc_mass_comp_ref["S_NH4"]
            )

        @self.Expression(self.flowsheet().time, doc="fSO2")
        def fSO2(b):
            return (
                b.conc_mass_comp_ref["S_O2"]
                / (
                    b.params.KA_AOB_ND_O2
                    + (1 - 2 * (b.params.KAOB_ND_O2 / b.params.KAOB_I_O2) ** 0.5)
                )
                * b.conc_mass_comp_ref["S_O2"]
                + b.conc_mass_comp_ref["S_O2"] ** 2 / b.params.KAOB_I_O2
            )

        try:

            def rate_expression_rule(b, r):
                if r == "R1":
                    # R1: Aerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * self.MonodX
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic hydrolysis (NO3-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.hL_NO3
                        * (
                            b.params.KL_O2
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KL_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * self.MonodX
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Anoxic hydrolysis (NO2-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.hL_NO2
                        * (
                            b.params.KL_O2
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.KL_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * self.MonodX
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Anaerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.hL_fe
                        * (
                            b.params.KL_O2
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.KL_NO3
                            / (
                                b.params.KL_NO3
                                + b.conc_mass_comp_ref["S_NO3"]
                                + b.conc_mass_comp_ref["S_NO2"]
                            )
                        )
                        * self.MonodX
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Aerobic growth on S_F
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Aerobic growth on S_A
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KH_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Anoxic growth on S_F (NO3- to NO2-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO3
                        * (
                            b.params.KH2_O2
                            / (b.params.KH2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F2 + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Anoxic growth on S_F (NO2- to NO)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO2
                        * (
                            b.params.KH3_O2
                            / (b.params.KH3_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F3 + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.KH_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Anoxic growth on S_F (NO to N2O)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO
                        * (
                            b.params.KH4_O2
                            / (b.params.KH4_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F4 + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (
                                b.params.KH_NO
                                + b.conc_mass_comp_ref["S_NO"]
                                + b.conc_mass_comp_ref["S_NO"] ** 2 / b.params.KI_NO
                            )
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R10":
                    # R10: Anoxic growth on S_F (N2O to N2)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_N2O
                        * (
                            b.params.KH5_O2
                            / (b.params.KH5_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F5 + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2O"]
                            / (b.params.KH_N2O + b.conc_mass_comp_ref["S_N2O"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Anoxic growth on S_A (NO3- to NO2-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO3
                        * (
                            b.params.KH2_O2
                            / (b.params.KH2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KH_A2 + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Anoxic growth on S_A (NO2- to NO)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO2
                        * (
                            b.params.KH3_O2
                            / (b.params.KH3_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KH_A3 + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.KH_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Anoxic growth on S_A (NO to N2O)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO
                        * (
                            b.params.KH4_O2
                            / (b.params.KH4_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KH_A4 + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (
                                b.params.KH_NO
                                + b.conc_mass_comp_ref["S_NO"]
                                + b.conc_mass_comp_ref["S_NO"] ** 2 / b.params.KI_NO
                            )
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Anoxic growth on S_A (N2O to N2)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_N2O
                        * (
                            b.params.KH5_O2
                            / (b.params.KH5_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KH_A5 + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                                + 1e-10 * pyo.units.kg / pyo.units.m**3
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2O"]
                            / (b.params.KH_N2O + b.conc_mass_comp_ref["S_N2O"])
                        )
                        * self.MonodNH4
                        * self.MonodPO4
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Fermentation
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_fe
                        * (
                            b.params.KH_O2
                            / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.KH_NO3
                            / (
                                b.params.KH_NO3
                                + b.conc_mass_comp_ref["S_NO3"]
                                + b.conc_mass_comp_ref["S_NO2"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_fe + b.conc_mass_comp_ref["S_F"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R16":
                    # R9: Lysis of X_H
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_H
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                                + b.params.hH_NO3_end
                                * b.params.KH_O2
                                / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                                * b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_H"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_H * b.conc_mass_comp_ref["X_H"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                elif r == "R17":
                    # R17: Storage of X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KP_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * self.MonodPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R18":
                    # R18: Aerobic storage of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * self.MonodIPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R19":
                    # R19: Anoxic storage of X_PP (NO3- to NO2-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * b.params.hP_NO3
                        * (
                            b.params.KP2_O2
                            / (b.params.KP2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * self.MonodIPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R20":
                    # R20: Anoxic storage of X_PP (NO2- to NO)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * b.params.hP_NO2
                        * (
                            b.params.KP3_O2
                            / (b.params.KP3_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.KP_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * self.MonodIPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R21":
                    # R21: Anoxic storage of X_PP (NO to N2O)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * b.params.hP_NO
                        * (
                            b.params.KP4_O2
                            / (b.params.KP4_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (
                                b.params.KP_NO
                                + b.conc_mass_comp_ref["S_NO"]
                                + b.conc_mass_comp_ref["S_NO"] ** 2 / b.params.KI_NO
                            )
                        )
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * self.MonodIPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R22":
                    # R22: Anoxic storage of X_PP (N2O to N2)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * b.params.hP_N2O
                        * (
                            b.params.KP5_O2
                            / (b.params.KP5_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2O"]
                            / (b.params.KP_N2O + b.conc_mass_comp_ref["S_N2O"])
                        )
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * self.MonodIPP
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R23":
                    # R23: Aerobic growth of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * self.MonodSNH4_P
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KP_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * self.MonodPHA
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R24":
                    # R24: Anoxic growth of X_PAO (NO3- to NO2-)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * b.params.hP_NO3
                        * (
                            b.params.KP2_O2
                            / (b.params.KP2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * self.MonodSNH4_P
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R25":
                    # R25: Anoxic growth of X_PAO (NO2- to NO)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * b.params.hP_NO2
                        * (
                            b.params.KP3_O2
                            / (b.params.KP3_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.KP_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * self.MonodSNH4_P
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R26":
                    # R26: Anoxic growth of X_PAO (NO to N2O)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * b.params.hP_NO
                        * (
                            b.params.KP4_O2
                            / (b.params.KP4_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (
                                b.params.KP_NO
                                + b.conc_mass_comp_ref["S_NO"]
                                + b.conc_mass_comp_ref["S_NO"] ** 2 / b.params.KI_NO
                            )
                        )
                        * self.MonodSNH4_P
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R27":
                    # R27: Anoxic growth of X_PAO (N2O to N2)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * b.params.hP_N2O
                        * (
                            b.params.KP5_O2
                            / (b.params.KP5_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2O"]
                            / (b.params.KP_N2O + b.conc_mass_comp_ref["S_N2O"])
                        )
                        * self.MonodSNH4_P
                        * self.MonodSPO4_P
                        * self.MonodPHA
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R28":
                    # R28: Lysis of X_PAO
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PAO
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                + b.params.hP_NO3_end
                                * b.params.KP_O2
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                * b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_PAO"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PAO * b.conc_mass_comp_ref["X_PAO"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                elif r == "R29":
                    # R29: Lysis of X_PP
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PP
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                + b.params.hPP_NO3_end
                                * b.params.KP_O2
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                * b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_PP"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PP * b.conc_mass_comp_ref["X_PP"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                elif r == "R30":
                    # R30: Lysis of X_PHA
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PHA
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                + b.params.hPHA_NO3_end
                                * b.params.KP_O2
                                / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                                * b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_PHA"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_PHA * b.conc_mass_comp_ref["X_PHA"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                elif r == "R31":
                    # R31: NH3 oxidation to NH2OH with O2 consumption
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_AOB_AMO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KAOB1_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (b.params.S_FA / (b.params.KAOB_NH4 + b.params.S_FA))
                        * b.conc_mass_comp_ref["X_AOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R32":
                    # R32: NH2OH to NO coupled with O2 reduction (X_AOB growth)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_AOB_AMO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KAOB2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH2OH"]
                            / (b.params.KAOB_NH2OH + b.conc_mass_comp_ref["S_NH2OH"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.conc_mass_comp_ref["S_NH4"] + self.eps * 1e-6)
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2"]
                            / (b.params.KAOB_P + b.conc_mass_comp_ref["S_N2"])
                        )
                        * b.conc_mass_comp_ref["X_AOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R33":
                    # R33: NO oxidation to NO2- with O2 reduction
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_AOB_HAO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KAOB2_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (b.params.KAOB_HAO_NO + b.conc_mass_comp_ref["S_NO"])
                        )
                        * b.conc_mass_comp_ref["X_AOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R34":
                    # R34: NO to N2O coupled with NH2OH to NO2- (N2O from NN pathway)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_AOB_N2O_NN
                        * (
                            b.conc_mass_comp_ref["S_NH2OH"]
                            / (b.params.KAOB_NH2OH + b.conc_mass_comp_ref["S_NH2OH"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (b.params.KAOB_NN_NO + b.conc_mass_comp_ref["S_NO"])
                        )
                        * b.conc_mass_comp_ref["X_AOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R35":
                    # R35: HNO2 to N2O coupled with NH2OH to NO2- (N2O from ND pathway)
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_AOB_N2O_ND
                        * (
                            b.conc_mass_comp_ref["S_NH2OH"]
                            / (b.params.KAOB_NH2OH + b.conc_mass_comp_ref["S_NH2OH"])
                        )
                        * (b.params.S_FNA / (b.params.KAOB_HNO2 + b.params.S_FNA))
                        * self.fSO2
                        * b.conc_mass_comp_ref["X_AOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R36":
                    # R36: Aerobic growth of X_NOB
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_NOB
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KNOB_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (b.params.S_FNA / (b.params.KNOB_NO2 + b.params.S_FNA))
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KNOB_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_NOB"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R37":
                    # R37: Lysis of X_AOB
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_AOB
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KAOB1_O2 + b.conc_mass_comp_ref["S_O2"])
                            )
                            + b.params.hAOB_NO3_end
                            * (
                                b.params.KAOB1_O2
                                / (b.params.KAOB1_O2 + b.conc_mass_comp_ref["S_O2"])
                            )
                            * (
                                b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KAOB_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_AOB"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_AOB * b.conc_mass_comp_ref["X_AOB"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                elif r == "R38":
                    # R38: Lysis of X_NOB
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_NOB
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KNOB_O2 + b.conc_mass_comp_ref["S_O2"])
                            )
                            + b.params.hNOB_NO3_end
                            * (
                                b.params.KNOB_O2
                                / (b.params.KNOB_O2 + b.conc_mass_comp_ref["S_O2"])
                            )
                            * (
                                b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KNOB_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_NOB"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_NOB * b.conc_mass_comp_ref["X_NOB"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                else:
                    raise BurntToast()

            self.rate_expression = pyo.Constraint(
                self.params.rate_reaction_idx,
                rule=rate_expression_rule,
                doc="ASM2d rate expressions",
            )

        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.rate_expression)
            raise

    def get_reaction_rate_basis(self):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.rate_expression.items():
            # TODO: Need to work out how to calculate good scaling factors
            # instead of a fixed 1e3.
            iscale.constraint_scaling_transform(c, 1e3, overwrite=True)
