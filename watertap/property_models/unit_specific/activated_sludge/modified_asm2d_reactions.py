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
Modified ASM2d reaction package.

Reference:

X. Flores-Alsina, K. Solon, C.K. Mbamba, S. Tait, K.V. Gernaey, U. Jeppsson, D.J. Batstone,
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions fordynamic simulations of anaerobic digestion processes,
Water Research. 95 (2016) 370-382. https://www.sciencedirect.com/science/article/pii/S0043135416301397
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
__author__ = "Marcus Holly, Adam Atia, Xinhong Liu"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ModifiedASM2dReactionParameterBlock")
class ModifiedASM2dReactionParameterData(ReactionParameterBlock):
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

        self._reaction_block_class = ModifiedASM2dReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ASM2d paper
        # R1: Aerobic hydrolysis
        # R2: Anoxic hydrolysis
        # R3: Anaerobic hydrolysis
        # R4: Aerobic growth on S_F
        # R5: Aerobic growth on S_A
        # R6: Anoxic growth on S_F
        # R7: Anoxic growth on S_A, denitrification
        # R8: Fermentation
        # R9: Lysis
        # R10: Storage of X_PHA
        # R11: Aerobic storage of X_PP
        # R12: Anoxic storage of X_PP
        # R13: Aerobic growth of X_PAO
        # R14: Anoxic growth of X_PAO
        # R15: Lysis of X_PAO
        # R16: Lysis of X_PP
        # R17: Lysis of X_PHA
        # R18: Aerobic growth of X_AUT
        # R19: Lysis of X_AUT
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
            ]
        )

        # Stoichiometric Parameters
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
        self.i_CXB = pyo.Var(
            initialize=0.36612,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="C content of inert soluble COD X_B, [kg C/kg COD]",
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
        self.f_XI = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD generated in lysis, [kg COD/kg COD]",
        )
        self.Y_PAO = pyo.Var(
            initialize=0.625,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield coefficient for P accumulating organisms (biomass/PHA), [kg COD/kg COD]",
        )
        self.Y_PO4 = pyo.Var(
            initialize=0.0129,
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
        self.Y_A = pyo.Var(
            initialize=0.24,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Yield of autotrophic biomass per NO3- N, [kg COD/kg N]",
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
        self.i_NOx_N2 = pyo.Var(
            initialize=2.8571,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Nitrogen oxide coefficient for N2",
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
            doc="Anoxic hydrolysis reduction factor",
        )
        self.hL_fe = pyo.Var(
            initialize=0.40,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anaerobic hydrolysis reduction factor",
        )
        self.KH_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen, [kg O2/m^3]",
        )
        self.KL_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen, [kg O2/m^3]",
        )
        self.KA_O2 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen, [kg O2/m^3]",
        )
        self.KP_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for oxygen, [kg O2/m^3]",
        )
        self.KH_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.KA_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.KL_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.KP_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for nitrate, [kg N/m^3]",
        )
        self.KL_X = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for particulate COD, [kg X_S/kg X_H]",
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
        self.b_H = pyo.Var(
            initialize=0.28,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for lysis and decay",
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
        self.KP_A = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for acetate SA, [kg COD/m^3]",
        )
        self.KH_NH4 = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for ammonium (nutrient), [kg N/m^3]",
        )
        self.KA_NH4 = pyo.Var(
            initialize=1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for SNH4, [kg N/m^3]",
        )
        self.KP_NH4 = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for ammonium, [kg N/m^3]",
        )
        self.KH_PO4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for SPO4 as nutrient, [kg P/m^3]",
        )
        self.KA_PO4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for SPO4, [kg P/m^3]",
        )
        self.KP_PO4 = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphate for growth, [kg P/m^3]",
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
            doc="Maximum growth rate of PAO",
        )
        self.b_PAO = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PAO",
        )
        self.b_PP = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PP",
        )
        self.b_PHA = pyo.Var(
            initialize=0.14,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PHA",
        )
        self.KP_P = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphate for XPP formation, [kg P/m^3]",
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

        self.mu_AUT = pyo.Var(
            initialize=0.61,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of X_AUT",
        )
        self.b_AUT = pyo.Var(
            initialize=0.09,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Decay rate of X_AUT",
        )
        self.hH_NO3 = pyo.Var(
            initialize=0.8,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor for denitrification",
        )
        self.hH_NO3_end = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for endogenous respiration",
        )
        self.hP_NO3 = pyo.Var(
            initialize=0.6,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Reduction factor under anoxic conditions",
        )
        self.hP_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PAOs",
        )
        self.hPP_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PP",
        )
        self.hPHA_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of PHA",
        )
        self.hAUT_NO3_end = pyo.Var(
            initialize=0.33,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic reduction factor for decay of autotrophs",
        )

        # Reaction Stoichiometry
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
            ("R1", "Liq", "S_N2"): 0,
            ("R1", "Liq", "S_NO3"): 0,
            ("R1", "Liq", "S_PO4"): -((1 - self.f_SI) * self.i_PSF - self.i_PXS),
            ("R1", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_S"): -1,
            ("R1", "Liq", "X_H"): 0,
            ("R1", "Liq", "X_PAO"): 0,
            ("R1", "Liq", "X_PP"): 0,
            ("R1", "Liq", "X_PHA"): 0,
            ("R1", "Liq", "X_AUT"): 0,
            ("R1", "Liq", "S_K"): 0,
            ("R1", "Liq", "S_Mg"): 0,
            # R2: Anoxic hydrolysis
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_O2"): 0,
            ("R2", "Liq", "S_F"): 1 - self.f_SI,
            ("R2", "Liq", "S_A"): 0,
            ("R2", "Liq", "S_I"): self.f_SI,
            ("R2", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R2", "Liq", "S_N2"): 0,
            ("R2", "Liq", "S_NO3"): 0,
            ("R2", "Liq", "S_PO4"): -((1 - self.f_SI) * self.i_PSF - self.i_PXS),
            ("R2", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_S"): -1,
            ("R2", "Liq", "X_H"): 0,
            ("R2", "Liq", "X_PAO"): 0,
            ("R2", "Liq", "X_PP"): 0,
            ("R2", "Liq", "X_PHA"): 0,
            ("R2", "Liq", "X_AUT"): 0,
            ("R2", "Liq", "S_K"): 0,
            ("R2", "Liq", "S_Mg"): 0,
            # R3: Anaerobic hydrolysis
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_O2"): 0,
            ("R3", "Liq", "S_F"): 1 - self.f_SI,
            ("R3", "Liq", "S_A"): 0,
            ("R3", "Liq", "S_I"): self.f_SI,
            ("R3", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R3", "Liq", "S_N2"): 0,
            ("R3", "Liq", "S_NO3"): 0,
            ("R3", "Liq", "S_PO4"): -((1 - self.f_SI) * self.i_PSF - self.i_PXS),
            ("R3", "Liq", "S_IC"): -(
                (1 - self.f_SI) * self.i_CSF + self.f_SI * self.i_CSI - self.i_CXS
            ),
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_S"): -1,
            ("R3", "Liq", "X_H"): 0,
            ("R3", "Liq", "X_PAO"): 0,
            ("R3", "Liq", "X_PP"): 0,
            ("R3", "Liq", "X_PHA"): 0,
            ("R3", "Liq", "X_AUT"): 0,
            ("R3", "Liq", "S_K"): 0,
            ("R3", "Liq", "S_Mg"): 0,
            # R4: Aerobic growth on S_F
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R4", "Liq", "S_F"): -1 / self.Y_H,
            ("R4", "Liq", "S_A"): 0,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H),
            ("R4", "Liq", "S_N2"): 0,
            ("R4", "Liq", "S_NO3"): 0,
            ("R4", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H),
            ("R4", "Liq", "S_IC"): -(self.i_CXB - self.i_CSF / self.Y_H),
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_S"): 0,
            ("R4", "Liq", "X_H"): 1,
            ("R4", "Liq", "X_PAO"): 0,
            ("R4", "Liq", "X_PP"): 0,
            ("R4", "Liq", "X_PHA"): 0,
            ("R4", "Liq", "X_AUT"): 0,
            ("R4", "Liq", "S_K"): 0,
            ("R4", "Liq", "S_Mg"): 0,
            # R5: Aerobic growth on S_A
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R5", "Liq", "S_F"): 0,
            ("R5", "Liq", "S_A"): -1 / self.Y_H,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "S_NH4"): -self.i_NBM,
            ("R5", "Liq", "S_N2"): 0,
            ("R5", "Liq", "S_NO3"): 0,
            ("R5", "Liq", "S_PO4"): -self.i_PBM,
            ("R5", "Liq", "S_IC"): -(self.i_CXB - self.i_CSA / self.Y_H),
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_S"): 0,
            ("R5", "Liq", "X_H"): 1,
            ("R5", "Liq", "X_PAO"): 0,
            ("R5", "Liq", "X_PP"): 0,
            ("R5", "Liq", "X_PHA"): 0,
            ("R5", "Liq", "X_AUT"): 0,
            ("R5", "Liq", "S_K"): 0,
            ("R5", "Liq", "S_Mg"): 0,
            # R6: Anoxic growth on S_F
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_O2"): 0,
            ("R6", "Liq", "S_F"): -1 / self.Y_H,
            ("R6", "Liq", "S_A"): 0,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H),
            ("R6", "Liq", "S_N2"): (1 - self.Y_H) / (2.86 * self.Y_H),
            ("R6", "Liq", "S_NO3"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R6", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H),
            ("R6", "Liq", "S_IC"): -(self.i_CXB - self.i_CSF / self.Y_H),
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_H"): 1,
            ("R6", "Liq", "X_PAO"): 0,
            ("R6", "Liq", "X_PP"): 0,
            ("R6", "Liq", "X_PHA"): 0,
            ("R6", "Liq", "X_AUT"): 0,
            ("R6", "Liq", "S_K"): 0,
            ("R6", "Liq", "S_Mg"): 0,
            # R7: Anoxic growth on S_A, denitrification
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_O2"): 0,
            ("R7", "Liq", "S_F"): 0,
            ("R7", "Liq", "S_A"): -1 / self.Y_H,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "S_NH4"): -self.i_NBM,
            ("R7", "Liq", "S_N2"): (1 - self.Y_H) / (2.86 * self.Y_H),
            ("R7", "Liq", "S_NO3"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R7", "Liq", "S_PO4"): -self.i_PBM,
            ("R7", "Liq", "S_IC"): -(self.i_CXB - self.i_CSA / self.Y_H),
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "X_S"): 0,
            ("R7", "Liq", "X_H"): 1,
            ("R7", "Liq", "X_PAO"): 0,
            ("R7", "Liq", "X_PP"): 0,
            ("R7", "Liq", "X_PHA"): 0,
            ("R7", "Liq", "X_AUT"): 0,
            ("R7", "Liq", "S_K"): 0,
            ("R7", "Liq", "S_Mg"): 0,
            # R8: Fermentation
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_O2"): 0,
            ("R8", "Liq", "S_F"): -1,
            ("R8", "Liq", "S_A"): 1,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "S_NH4"): -self.i_NSF,
            ("R8", "Liq", "S_N2"): 0,
            ("R8", "Liq", "S_NO3"): 0,
            ("R8", "Liq", "S_PO4"): -(-self.i_PSF),
            ("R8", "Liq", "S_IC"): -(self.i_CSA - self.i_CSF),
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_S"): 0,
            ("R8", "Liq", "X_H"): 0,
            ("R8", "Liq", "X_PAO"): 0,
            ("R8", "Liq", "X_PP"): 0,
            ("R8", "Liq", "X_PHA"): 0,
            ("R8", "Liq", "X_AUT"): 0,
            ("R8", "Liq", "S_K"): 0,
            ("R8", "Liq", "S_Mg"): 0,
            # R9: Lysis
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_O2"): 0,
            ("R9", "Liq", "S_F"): 0,
            ("R9", "Liq", "S_A"): 0,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R9", "Liq", "S_N2"): 0,
            ("R9", "Liq", "S_NO3"): 0,
            ("R9", "Liq", "S_PO4"): -(
                self.f_XI * self.i_PXI + (1 - self.f_XI) * self.i_PXS - self.i_PBM
            ),
            ("R9", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CXB
            ),
            ("R9", "Liq", "X_I"): self.f_XI,
            ("R9", "Liq", "X_S"): 1 - self.f_XI,
            ("R9", "Liq", "X_H"): -1,
            ("R9", "Liq", "X_PAO"): 0,
            ("R9", "Liq", "X_PP"): 0,
            ("R9", "Liq", "X_PHA"): 0,
            ("R9", "Liq", "X_AUT"): 0,
            ("R9", "Liq", "S_K"): 0,
            ("R9", "Liq", "S_Mg"): 0,
            # R10: Storage of X_PHA
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_O2"): 0,
            ("R10", "Liq", "S_F"): 0,
            ("R10", "Liq", "S_A"): -1,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "S_NH4"): 0,
            ("R10", "Liq", "S_N2"): 0,
            ("R10", "Liq", "S_NO3"): 0,
            ("R10", "Liq", "S_PO4"): -(-self.Y_PO4),
            ("R10", "Liq", "S_IC"): -(0.3 - self.i_CSA),
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "X_S"): 0,
            ("R10", "Liq", "X_H"): 0,
            ("R10", "Liq", "X_PAO"): 0,
            ("R10", "Liq", "X_PP"): -self.Y_PO4,
            ("R10", "Liq", "X_PHA"): 1,
            ("R10", "Liq", "X_AUT"): 0,
            ("R10", "Liq", "S_K"): self.Y_PO4 * self.i_KXPP,
            ("R10", "Liq", "S_Mg"): self.Y_PO4 * self.i_MgXPP,
            # R11: Aerobic storage of X_PP
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_O2"): -self.Y_PHA,
            ("R11", "Liq", "S_F"): 0,
            ("R11", "Liq", "S_A"): 0,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "S_NH4"): 0,
            ("R11", "Liq", "S_N2"): 0,
            ("R11", "Liq", "S_NO3"): 0,
            ("R11", "Liq", "S_PO4"): -1,
            ("R11", "Liq", "S_IC"): -(-self.Y_PHA * 0.3),
            ("R11", "Liq", "X_I"): 0,
            ("R11", "Liq", "X_S"): 0,
            ("R11", "Liq", "X_H"): 0,
            ("R11", "Liq", "X_PAO"): 0,
            ("R11", "Liq", "X_PP"): 1,
            ("R11", "Liq", "X_PHA"): -self.Y_PHA,
            ("R11", "Liq", "X_AUT"): 0,
            ("R11", "Liq", "S_K"): -self.i_KXPP,
            ("R11", "Liq", "S_Mg"): -self.i_MgXPP,
            # R12: Anoxic storage of X_PP
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_O2"): 0,
            ("R12", "Liq", "S_F"): 0,
            ("R12", "Liq", "S_A"): 0,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "S_NH4"): 0,
            ("R12", "Liq", "S_N2"): self.Y_PHA / 2.86,
            ("R12", "Liq", "S_NO3"): -0.07,
            ("R12", "Liq", "S_PO4"): -1,
            ("R12", "Liq", "S_IC"): -(-self.Y_PHA * 0.3),
            ("R12", "Liq", "X_I"): 0,
            ("R12", "Liq", "X_S"): 0,
            ("R12", "Liq", "X_H"): 0,
            ("R12", "Liq", "X_PAO"): 0,
            ("R12", "Liq", "X_PP"): 1,
            ("R12", "Liq", "X_PHA"): -self.Y_PHA,
            ("R12", "Liq", "X_AUT"): 0,
            ("R12", "Liq", "S_K"): -self.i_KXPP,
            ("R12", "Liq", "S_Mg"): -self.i_MgXPP,
            # R13: Aerobic growth of X_PAO
            ("R13", "Liq", "H2O"): 0,
            ("R13", "Liq", "S_O2"): -(1 - self.Y_PAO) / self.Y_PAO,
            ("R13", "Liq", "S_F"): 0,
            ("R13", "Liq", "S_A"): 0,
            ("R13", "Liq", "S_I"): 0,
            ("R13", "Liq", "S_NH4"): -self.i_NBM,
            ("R13", "Liq", "S_N2"): 0,
            ("R13", "Liq", "S_NO3"): 0,
            ("R13", "Liq", "S_PO4"): -self.i_PBM,
            ("R13", "Liq", "S_IC"): -(self.i_CXB - 0.3 / self.Y_PAO),
            ("R13", "Liq", "X_I"): 0,
            ("R13", "Liq", "X_S"): 0,
            ("R13", "Liq", "X_H"): 0,
            ("R13", "Liq", "X_PAO"): 1,
            ("R13", "Liq", "X_PP"): 0,
            ("R13", "Liq", "X_PHA"): -1 / self.Y_PAO,
            ("R13", "Liq", "X_AUT"): 0,
            ("R13", "Liq", "S_K"): 0,
            ("R13", "Liq", "S_Mg"): 0,
            # R14: Anoxic growth of X_PAO
            ("R14", "Liq", "H2O"): 0,
            ("R14", "Liq", "S_O2"): 0,
            ("R14", "Liq", "S_F"): 0,
            ("R14", "Liq", "S_A"): 0,
            ("R14", "Liq", "S_I"): 0,
            ("R14", "Liq", "S_NH4"): -self.i_NBM,
            ("R14", "Liq", "S_N2"): (1 - self.Y_H) / (2.86 * self.Y_H),
            ("R14", "Liq", "S_NO3"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R14", "Liq", "S_PO4"): -self.i_PBM,
            ("R14", "Liq", "S_IC"): -(self.i_CXB - 0.3 / self.Y_PAO),
            ("R14", "Liq", "X_I"): 0,
            ("R14", "Liq", "X_S"): 0,
            ("R14", "Liq", "X_H"): 0,
            ("R14", "Liq", "X_PAO"): 1,
            ("R14", "Liq", "X_PP"): 0,
            ("R14", "Liq", "X_PHA"): -1 / self.Y_PAO,
            ("R14", "Liq", "X_AUT"): 0,
            ("R14", "Liq", "S_K"): 0,
            ("R14", "Liq", "S_Mg"): 0,
            # R15: Lysis of X_PAO
            ("R15", "Liq", "H2O"): 0,
            ("R15", "Liq", "S_O2"): 0,
            ("R15", "Liq", "S_F"): 0,
            ("R15", "Liq", "S_A"): 0,
            ("R15", "Liq", "S_I"): 0,
            ("R15", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R15", "Liq", "S_N2"): 0,
            ("R15", "Liq", "S_NO3"): 0,
            ("R15", "Liq", "S_PO4"): -(
                self.f_XI * self.i_PXI + (1 - self.f_XI) * self.i_PXS - self.i_PBM
            ),
            ("R15", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CXB
            ),
            ("R15", "Liq", "X_I"): self.f_XI,
            ("R15", "Liq", "X_S"): 1 - self.f_XI,
            ("R15", "Liq", "X_H"): 0,
            ("R15", "Liq", "X_PAO"): -1,
            ("R15", "Liq", "X_PP"): 0,
            ("R15", "Liq", "X_PHA"): 0,
            ("R15", "Liq", "X_AUT"): 0,
            ("R15", "Liq", "S_K"): 0,
            ("R15", "Liq", "S_Mg"): 0,
            # R16: Lysis of X_PP
            ("R16", "Liq", "H2O"): 0,
            ("R16", "Liq", "S_O2"): 0,
            ("R16", "Liq", "S_F"): 0,
            ("R16", "Liq", "S_A"): 0,
            ("R16", "Liq", "S_I"): 0,
            ("R16", "Liq", "S_NH4"): 0,
            ("R16", "Liq", "S_N2"): 0,
            ("R16", "Liq", "S_NO3"): 0,
            ("R16", "Liq", "S_PO4"): 1,
            ("R16", "Liq", "S_IC"): 0,
            ("R16", "Liq", "X_I"): 0,
            ("R16", "Liq", "X_S"): 0,
            ("R16", "Liq", "X_H"): 0,
            ("R16", "Liq", "X_PAO"): 0,
            ("R16", "Liq", "X_PP"): -1,
            ("R16", "Liq", "X_PHA"): 0,
            ("R16", "Liq", "X_AUT"): 0,
            ("R16", "Liq", "S_K"): self.i_KXPP,
            ("R16", "Liq", "S_Mg"): self.i_MgXPP,
            # R17: Lysis of X_PHA
            ("R17", "Liq", "H2O"): 0,
            ("R17", "Liq", "S_O2"): 0,
            ("R17", "Liq", "S_F"): 0,
            ("R17", "Liq", "S_A"): 1,
            ("R17", "Liq", "S_I"): 0,
            ("R17", "Liq", "S_NH4"): 0,
            ("R17", "Liq", "S_N2"): 0,
            ("R17", "Liq", "S_NO3"): 0,
            ("R17", "Liq", "S_PO4"): 0,
            ("R17", "Liq", "S_IC"): -(self.i_CSA - 0.3),
            ("R17", "Liq", "X_I"): 0,
            ("R17", "Liq", "X_S"): 0,
            ("R17", "Liq", "X_H"): 0,
            ("R17", "Liq", "X_PAO"): 0,
            ("R17", "Liq", "X_PP"): 0,
            ("R17", "Liq", "X_PHA"): -1,
            ("R17", "Liq", "X_AUT"): 0,
            ("R17", "Liq", "S_K"): 0,
            ("R17", "Liq", "S_Mg"): 0,
            # R18: Aerobic growth of X_AUT
            ("R18", "Liq", "H2O"): 0,
            ("R18", "Liq", "S_O2"): -(4.5714 - self.Y_A) / self.Y_A,
            ("R18", "Liq", "S_F"): 0,
            ("R18", "Liq", "S_A"): 0,
            ("R18", "Liq", "S_I"): 0,
            ("R18", "Liq", "S_NH4"): -1 / self.Y_A - self.i_NBM,
            ("R18", "Liq", "S_N2"): 0,
            ("R18", "Liq", "S_NO3"): 1 / self.Y_A,
            ("R18", "Liq", "S_PO4"): -self.i_PBM,
            ("R18", "Liq", "S_IC"): -self.i_CXB,
            ("R18", "Liq", "X_I"): 0,
            ("R18", "Liq", "X_S"): 0,
            ("R18", "Liq", "X_H"): 0,
            ("R18", "Liq", "X_PAO"): 0,
            ("R18", "Liq", "X_PP"): 0,
            ("R18", "Liq", "X_PHA"): 0,
            ("R18", "Liq", "X_AUT"): 1,
            ("R18", "Liq", "S_K"): 0,
            ("R18", "Liq", "S_Mg"): 0,
            # R19: Lysis of X_AUT
            ("R19", "Liq", "H2O"): 0,
            ("R19", "Liq", "S_O2"): 0,
            ("R19", "Liq", "S_F"): 0,
            ("R19", "Liq", "S_A"): 0,
            ("R19", "Liq", "S_I"): 0,
            ("R19", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R19", "Liq", "S_N2"): 0,
            ("R19", "Liq", "S_NO3"): 0,
            ("R19", "Liq", "S_PO4"): -(
                self.f_XI * self.i_PXI + (1 - self.f_XI) * self.i_PXS - self.i_PBM
            ),
            ("R19", "Liq", "S_IC"): -(
                self.f_XI * self.i_CXI + (1 - self.f_XI) * self.i_CXS - self.i_CXB
            ),
            ("R19", "Liq", "X_I"): self.f_XI,
            ("R19", "Liq", "X_S"): 1 - self.f_XI,
            ("R19", "Liq", "X_H"): 0,
            ("R19", "Liq", "X_PAO"): 0,
            ("R19", "Liq", "X_PP"): 0,
            ("R19", "Liq", "X_PHA"): 0,
            ("R19", "Liq", "X_AUT"): -1,
            ("R19", "Liq", "S_K"): 0,
            ("R19", "Liq", "S_Mg"): 0,
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


class ModifiedASM2dReactionScaler(CustomScalerBase):
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


class _ModifiedASM2dReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ModifiedASM2dReactionScaler

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
    "ModifiedASM2dReactionBlock", block_class=_ModifiedASM2dReactionBlock
)
class ModifiedASM2dReactionBlockData(ReactionBlockDataBase):
    """
    Reaction Block for ASM2d.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super().build()

        # Create references to state vars
        # Concentration
        add_object_reference(self, "conc_mass_comp_ref", self.state_ref.conc_mass_comp)

    # Rate of reaction method
    def _rxn_rate(self):
        self.reaction_rate = pyo.Var(
            self.params.rate_reaction_idx,
            initialize=0,
            doc="Rate of reaction",
            units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
        )

        try:
            # TODO: Refer to Flores-Alsina c code to account for sulfur in rate expressions
            def rate_expression_rule(b, r):
                if r == "R1":
                    # R1: Aerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_S"]
                            / (
                                b.params.KL_X * b.conc_mass_comp_ref["X_H"]
                                + b.conc_mass_comp_ref["X_S"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic hydrolysis
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
                        * (
                            b.conc_mass_comp_ref["X_S"]
                            / (
                                b.params.KL_X * b.conc_mass_comp_ref["X_H"]
                                + b.conc_mass_comp_ref["X_S"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Anaerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.hL_fe
                        * (
                            b.params.KL_O2
                            / (b.params.KL_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.KL_NO3
                            / (b.params.KL_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_S"]
                            / (
                                b.params.KL_X * b.conc_mass_comp_ref["X_H"]
                                + b.conc_mass_comp_ref["X_S"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Aerobic growth on S_F
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
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KH_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KH_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Aerobic growth on S_A
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
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KH_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KH_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Anoxic growth on S_F
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO3
                        * (
                            b.params.KH_O2
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
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KH_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KH_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Anoxic growth on S_A, denitrification
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.hH_NO3
                        * (
                            b.params.KH_O2
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
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KH_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KH_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Fermentation
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_fe
                        * (
                            b.params.KH_O2
                            / (b.params.KH_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.KH_NO3
                            / (b.params.KH_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_fe + b.conc_mass_comp_ref["S_F"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Lysis
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
                elif r == "R10":
                    # R10: Storage of X_PHA
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PHA
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.KP_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_PP"]
                            / (
                                b.params.KP_PP * b.conc_mass_comp_ref["X_PAO"]
                                + b.conc_mass_comp_ref["X_PP"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Aerobic storage of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KP_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_PHA"]
                            / (
                                b.params.KP_PHA * b.conc_mass_comp_ref["X_PAO"]
                                + b.conc_mass_comp_ref["X_PHA"]
                            )
                        )
                        * (
                            (
                                b.params.K_MAX * b.conc_mass_comp_ref["X_PAO"]
                                - b.conc_mass_comp_ref["X_PP"]
                            )
                            / (
                                b.params.KI_PP * b.conc_mass_comp_ref["X_PAO"]
                                + b.params.K_MAX * b.conc_mass_comp_ref["X_PAO"]
                                - b.conc_mass_comp_ref["X_PP"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Anoxic storage of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_PP
                        * b.params.hP_NO3
                        * (
                            b.params.KP_O2
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KP_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_PHA"]
                            / (
                                b.params.KP_PHA * b.conc_mass_comp_ref["X_PAO"]
                                + b.conc_mass_comp_ref["X_PHA"]
                            )
                        )
                        * (
                            (
                                b.params.K_MAX * b.conc_mass_comp_ref["X_PAO"]
                                - b.conc_mass_comp_ref["X_PP"]
                            )
                            / (
                                b.params.KI_PP * b.conc_mass_comp_ref["X_PAO"]
                                + b.params.K_MAX * b.conc_mass_comp_ref["X_PAO"]
                                - b.conc_mass_comp_ref["X_PP"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Aerobic growth of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KP_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KP_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_PHA"]
                            / (
                                b.params.KP_PHA * b.conc_mass_comp_ref["X_PAO"]
                                + b.conc_mass_comp_ref["X_PHA"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Anoxic growth of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * b.params.hP_NO3
                        * (
                            b.params.KP_O2
                            / (b.params.KP_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.KP_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KP_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KP_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.conc_mass_comp_ref["X_PHA"]
                            / (
                                b.params.KP_PHA * b.conc_mass_comp_ref["X_PAO"]
                                + b.conc_mass_comp_ref["X_PHA"]
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Lysis of X_PAO
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
                elif r == "R16":
                    # R16: Lysis of X_PP
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
                elif r == "R17":
                    # R17: Lysis of X_PHA
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
                elif r == "R18":
                    # R18: Lysis of X_AUT
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_AUT
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.KA_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.KA_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.KA_PO4 + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_AUT"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R19":
                    # R19: Aerobic growth of X_AUT
                    if self.params.config.decay_switch:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_AUT
                            * (
                                b.conc_mass_comp_ref["S_O2"]
                                / (b.params.KA_O2 + b.conc_mass_comp_ref["S_O2"])
                                + b.params.hAUT_NO3_end
                                * b.params.KA_O2
                                / (b.params.KA_O2 + b.conc_mass_comp_ref["S_O2"])
                                * b.conc_mass_comp_ref["S_NO3"]
                                / (b.params.KA_NO3 + b.conc_mass_comp_ref["S_NO3"])
                            )
                            * b.conc_mass_comp_ref["X_AUT"],
                            to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                        )
                    else:
                        return b.reaction_rate[r] == pyo.units.convert(
                            b.params.b_AUT * b.conc_mass_comp_ref["X_AUT"],
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
