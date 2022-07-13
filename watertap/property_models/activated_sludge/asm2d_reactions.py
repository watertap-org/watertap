#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
ASM2d reaction package.

Important Note: ASM2d reactions depend on the presences of TSS in solution, however
TSS does not take part in the rate expressions. Thus, it is possible to have cases
where there is insufficient TSS present in the system for the reactions to occur
resulting in an infeasible solution.

Reference:

[1] Henze, M., Gujer, W., Mino, T., Matsuo, T., Wentzel, M.C., Marais, G.v.R.,
Van Loosdrecht, M.C.M., "Activated Sludge Model No.2D, ASM2D", 1999,
Wat. Sci. Tech. Vol. 39, No. 1, pp. 165-182
"""

# Import Pyomo libraries
import pyomo.environ as pyo

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


# Some more information about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASM2dReactionParameterBlock")
class ASM2dReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ASM2dReactionBlock

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
        # R20: Precipitation
        # R21: Re-dissolution
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
            ]
        )

        # Stoichiometric Parameters
        self.i_NSI = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert soluble COD S_I, [kg N/kg COD]",
        )
        self.i_NSF = pyo.Var(
            initialize=0.03,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of fermentable substrate, S_F, [kg N/kg COD]",
        )
        self.i_NXI = pyo.Var(
            initialize=0.02,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert particulate COD X_I, [kg N/kg COD]",
        )
        self.i_NXS = pyo.Var(
            initialize=0.04,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of slowly biodegradable substrate X_S, [kg N/kg COD]",
        )
        self.i_NBM = pyo.Var(
            initialize=0.07,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of biomass, X_H, X_PAO, X_AUT, [kg N/kg COD]",
        )
        self.i_PSI = pyo.Var(
            initialize=0.00,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert soluble COD S_I, [kg P/kg COD]",
        )
        self.i_PSF = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of fermentable substrate, S_F, [kg P/kg COD]",
        )
        self.i_PXI = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert particulate COD X_I, [kg P/kg COD]",
        )
        self.i_PXS = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of slowly biodegradable substrate X_S, [kg P/kg COD]",
        )
        self.i_PBM = pyo.Var(
            initialize=0.02,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of biomass, X_H, X_PAO, X_AUT, [kg P/kg COD]",
        )
        self.i_TSSXI = pyo.Var(
            initialize=0.75,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="TSS content of inert particulate COD X_I, [kg TSS/kg COD]",
        )
        self.i_TSSXS = pyo.Var(
            initialize=0.75,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="TSS content of slowly biodegradable substrate X_S, [kg TSS/kg COD]",
        )
        self.i_TSSBM = pyo.Var(
            initialize=0.90,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="TSS content of biomass, X_H, X_PAO, X_AUT, [kg TSS/kg COD]",
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
            initialize=0.40,
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

        # Kinetic Parameters
        self.K_H = pyo.Var(
            initialize=3,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Hydrolysis rate constant",
        )
        self.eta_NO3 = pyo.Var(
            initialize=0.60,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anoxic hydrolysis reduction factor",
        )
        self.eta_fe = pyo.Var(
            initialize=0.40,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Anaerobic hydrolysis reduction factor",
        )
        self.K_O2 = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for oxygen, [kg O2/m^3]",
        )
        self.K_NO3 = pyo.Var(
            initialize=5e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation/inhibition coefficient for nitrate, [kg N/m^3]",
        )
        self.K_X = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for particulate COD, [kg X_S/kg X_H]",
        )

        self.mu_H = pyo.Var(
            initialize=6,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate on substrate, [kg X_S/kg X_H/day]",
        )
        self.q_fe = pyo.Var(
            initialize=3,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum rate for fermentation, [kg S_F/kg X_H/day]",
        )
        self.b_H = pyo.Var(
            initialize=0.4,
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
        self.K_A = pyo.Var(
            initialize=4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for growth on acetate SA, [kg COD/m^3]",
        )
        self.K_NH4 = pyo.Var(
            initialize=5e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for ammonium (nutrient), [kg N/m^3]",
        )
        self.K_P = pyo.Var(
            initialize=1e-5,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphate (nutrient), [kg P/m^3]",
        )
        self.K_ALK = pyo.Var(
            initialize=1e-4,
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for alkalinity (HCO3-), [kmol HCO3-/m^3]",
        )

        self.q_PHA = pyo.Var(
            initialize=3,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for storage of X_PHA (base Xpp), [kg PHA/kg PAO/day]",
        )
        self.q_PP = pyo.Var(
            initialize=1.5,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for storage of X_PP, [kg PP/kg PAO/day]",
        )
        self.mu_PAO = pyo.Var(
            initialize=1,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of PAO",
        )
        self.b_PAO = pyo.Var(
            initialize=0.2,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PAO",
        )
        self.b_PP = pyo.Var(
            initialize=0.2,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PP",
        )
        self.b_PHA = pyo.Var(
            initialize=0.2,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate for Lysis of X_PHA",
        )
        self.K_PS = pyo.Var(
            initialize=2e-4,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for phosphorus in storage of PP, [kg P/m^3]",
        )
        self.K_PP = pyo.Var(
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
        self.K_IPP = pyo.Var(
            initialize=0.02,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Inhibition coefficient for PP storage, [kg PP/kg PAO]",
        )
        self.K_PHA = pyo.Var(
            initialize=0.01,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Saturation coefficient for PHA, [kg PHA/kg PAO]",
        )

        self.mu_AUT = pyo.Var(
            initialize=1,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Maximum growth rate of X_AUT",
        )
        self.b_AUT = pyo.Var(
            initialize=0.15,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Decay rate of X_AUT",
        )

        self.k_pre = pyo.Var(
            initialize=1e3,
            units=pyo.units.m**3 / pyo.units.kg / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for P precipitation, [m^3/kg Fe(OH)3/day",
        )
        self.k_red = pyo.Var(
            initialize=0.6,
            units=1 / pyo.units.day,
            domain=pyo.NonNegativeReals,
            doc="Rate constant for redissolution",
        )

        # Reaction Stoichiometry
        # This is the stoichiometric part the Peterson matrix in dict form
        # Note that reaction stoichiometry is on a mass basis.
        # For alkalinity, this requires converting the mass of species
        # reacted to mass of alkalinity converted using a charge balance
        mw_alk = 61 * pyo.units.kg / pyo.units.kmol
        mw_N = 14 * pyo.units.kg / pyo.units.kmol
        mw_P = 31 * pyo.units.kg / pyo.units.kmol
        self.rate_reaction_stoichiometry = {
            # R1: Aerobic hydrolysis
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_A"): 0,
            ("R1", "Liq", "S_F"): 1 - self.f_SI,
            ("R1", "Liq", "S_I"): self.f_SI,
            ("R1", "Liq", "S_N2"): 0,
            ("R1", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R1", "Liq", "S_NO3"): 0,
            ("R1", "Liq", "S_O2"): 0,
            ("R1", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R1", "Liq", "S_ALK"): (
                -((1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS)
                * (mw_alk / mw_N)
                - ((1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS)
                * (1.5 * mw_alk / mw_P)
            ),
            ("R1", "Liq", "X_AUT"): 0,
            ("R1", "Liq", "X_H"): 0,
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_MeOH"): 0,
            ("R1", "Liq", "X_MeP"): 0,
            ("R1", "Liq", "X_PAO"): 0,
            ("R1", "Liq", "X_PHA"): 0,
            ("R1", "Liq", "X_PP"): 0,
            ("R1", "Liq", "X_S"): -1,
            ("R1", "Liq", "X_TSS"): -self.i_TSSXS,
            # R2: Anoxic hydrolysis
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_A"): 0,
            ("R2", "Liq", "S_F"): 1 - self.f_SI,
            ("R2", "Liq", "S_I"): self.f_SI,
            ("R2", "Liq", "S_N2"): 0,
            ("R2", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R2", "Liq", "S_NO3"): 0,
            ("R2", "Liq", "S_O2"): 0,
            ("R2", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R2", "Liq", "S_ALK"): (
                -((1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS)
                * (mw_alk / mw_N)
                - ((1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS)
                * (1.5 * mw_alk / mw_P)
            ),
            ("R2", "Liq", "X_AUT"): 0,
            ("R2", "Liq", "X_H"): 0,
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_MeOH"): 0,
            ("R2", "Liq", "X_MeP"): 0,
            ("R2", "Liq", "X_PAO"): 0,
            ("R2", "Liq", "X_PHA"): 0,
            ("R2", "Liq", "X_PP"): 0,
            ("R2", "Liq", "X_S"): -1,
            ("R2", "Liq", "X_TSS"): -self.i_TSSXS,
            # R3: Anaerobic hydrolysis
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_A"): 0,
            ("R3", "Liq", "S_F"): 1 - self.f_SI,
            ("R3", "Liq", "S_I"): self.f_SI,
            ("R3", "Liq", "S_N2"): 0,
            ("R3", "Liq", "S_NH4"): -(
                (1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS
            ),
            ("R3", "Liq", "S_NO3"): 0,
            ("R3", "Liq", "S_O2"): 0,
            ("R3", "Liq", "S_PO4"): -(
                (1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS
            ),
            ("R3", "Liq", "S_ALK"): (
                -((1 - self.f_SI) * self.i_NSF + self.f_SI * self.i_NSI - self.i_NXS)
                * (mw_alk / mw_N)
                - ((1 - self.f_SI) * self.i_PSF + self.f_SI * self.i_PSI - self.i_PXS)
                * (1.5 * mw_alk / mw_P)
            ),
            ("R3", "Liq", "X_AUT"): 0,
            ("R3", "Liq", "X_H"): 0,
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_MeOH"): 0,
            ("R3", "Liq", "X_MeP"): 0,
            ("R3", "Liq", "X_PAO"): 0,
            ("R3", "Liq", "X_PHA"): 0,
            ("R3", "Liq", "X_PP"): 0,
            ("R3", "Liq", "X_S"): -1,
            ("R3", "Liq", "X_TSS"): -self.i_TSSXS,
            # R4: Aerobic growth on S_F
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_A"): 0,
            ("R4", "Liq", "S_F"): -1 / self.Y_H,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "S_N2"): 0,
            ("R4", "Liq", "S_NH4"): -(self.i_NBM - self.i_NSF / self.Y_H),
            ("R4", "Liq", "S_NO3"): 0,
            ("R4", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R4", "Liq", "S_PO4"): -(self.i_PBM - self.i_PSF / self.Y_H),
            ("R4", "Liq", "S_ALK"): -(self.i_NBM - self.i_NSF / self.Y_H)
            * mw_alk
            / mw_N
            + (self.i_PBM - self.i_PSF / self.Y_H) * 1.5 * mw_alk / mw_P,
            ("R4", "Liq", "X_AUT"): 0,
            ("R4", "Liq", "X_H"): 1,
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_MeOH"): 0,
            ("R4", "Liq", "X_MeP"): 0,
            ("R4", "Liq", "X_PAO"): 0,
            ("R4", "Liq", "X_PHA"): 0,
            ("R4", "Liq", "X_PP"): 0,
            ("R4", "Liq", "X_S"): 0,
            ("R4", "Liq", "X_TSS"): self.i_TSSBM,
            # R5: Aerobic growth on S_A
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_A"): -1 / self.Y_H,
            ("R5", "Liq", "S_F"): 0,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "S_N2"): 0,
            ("R5", "Liq", "S_NH4"): -self.i_NBM,
            ("R5", "Liq", "S_NO3"): 0,
            ("R5", "Liq", "S_O2"): 1 - 1 / self.Y_H,
            ("R5", "Liq", "S_PO4"): -self.i_PBM,
            ("R5", "Liq", "S_ALK"): -self.i_NBM * mw_alk / mw_N
            + self.i_PBM * 1.5 * mw_alk / mw_P
            + (1 / self.Y_H) * mw_alk / (64 * pyo.units.kg / pyo.units.kmol),
            ("R5", "Liq", "X_AUT"): 0,
            ("R5", "Liq", "X_H"): 1,
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_MeOH"): 0,
            ("R5", "Liq", "X_MeP"): 0,
            ("R5", "Liq", "X_PAO"): 0,
            ("R5", "Liq", "X_PHA"): 0,
            ("R5", "Liq", "X_PP"): 0,
            ("R5", "Liq", "X_S"): 0,
            ("R5", "Liq", "X_TSS"): self.i_TSSBM,
            # R6: Anoxic growth on S_F
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_A"): 0,
            ("R6", "Liq", "S_F"): -1 / self.Y_H,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "S_N2"): (1 - self.Y_H) / (2.86 * self.Y_H),
            ("R6", "Liq", "S_NH4"): -self.i_NBM + self.i_NSF / self.Y_H,
            ("R6", "Liq", "S_NO3"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R6", "Liq", "S_O2"): 0,
            ("R6", "Liq", "S_PO4"): -self.i_PBM + self.i_PSF / self.Y_H,
            ("R6", "Liq", "S_ALK"): (
                (-self.i_NBM + self.i_NSF / self.Y_H)
                + (1 - self.Y_H) / (2.86 * self.Y_H)
            )
            * mw_alk
            / mw_N
            - (-self.i_PBM + self.i_PSF / self.Y_H) * 1.5 * mw_alk / mw_P,
            ("R6", "Liq", "X_AUT"): 0,
            ("R6", "Liq", "X_H"): 1,
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "X_MeOH"): 0,
            ("R6", "Liq", "X_MeP"): 0,
            ("R6", "Liq", "X_PAO"): 0,
            ("R6", "Liq", "X_PHA"): 0,
            ("R6", "Liq", "X_PP"): 0,
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_TSS"): self.i_TSSBM,
            # R7: Anoxic growth on S_A, denitrification
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_A"): -1 / self.Y_H,
            ("R7", "Liq", "S_F"): 0,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "S_N2"): (1 - self.Y_H) / (2.86 * self.Y_H),
            ("R7", "Liq", "S_NH4"): -self.i_NBM,
            ("R7", "Liq", "S_NO3"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R7", "Liq", "S_O2"): 0,
            ("R7", "Liq", "S_PO4"): -self.i_PBM,
            ("R7", "Liq", "S_ALK"): (1 / self.Y_H)
            * mw_alk
            / (64 * pyo.units.kg / pyo.units.kmol)
            + (-self.i_NBM + (1 - self.Y_H) / (2.86 * self.Y_H)) * mw_alk / mw_N
            + self.i_PBM * 1.5 * mw_alk / mw_P,
            ("R7", "Liq", "X_AUT"): 0,
            ("R7", "Liq", "X_H"): 1,
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "X_MeOH"): 0,
            ("R7", "Liq", "X_MeP"): 0,
            ("R7", "Liq", "X_PAO"): 0,
            ("R7", "Liq", "X_PHA"): 0,
            ("R7", "Liq", "X_PP"): 0,
            ("R7", "Liq", "X_S"): 0,
            ("R7", "Liq", "X_TSS"): self.i_TSSBM,
            # R8: Fermentation
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_A"): 1,
            ("R8", "Liq", "S_F"): -1,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "S_N2"): 0,
            ("R8", "Liq", "S_NH4"): self.i_NSF,
            ("R8", "Liq", "S_NO3"): 0,
            ("R8", "Liq", "S_O2"): 0,
            ("R8", "Liq", "S_PO4"): self.i_PSF,
            ("R8", "Liq", "S_ALK"): -1 / 64
            + self.i_NSF * mw_alk / mw_N
            - self.i_PSF * 1.5 * mw_alk / mw_P,
            ("R8", "Liq", "X_AUT"): 0,
            ("R8", "Liq", "X_H"): 0,
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_MeOH"): 0,
            ("R8", "Liq", "X_MeP"): 0,
            ("R8", "Liq", "X_PAO"): 0,
            ("R8", "Liq", "X_PHA"): 0,
            ("R8", "Liq", "X_PP"): 0,
            ("R8", "Liq", "X_S"): 0,
            ("R8", "Liq", "X_TSS"): 0,
            # R9: Lysis
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_A"): 0,
            ("R9", "Liq", "S_F"): 0,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "S_N2"): 0,
            ("R9", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R9", "Liq", "S_NO3"): 0,
            ("R9", "Liq", "S_O2"): 0,
            ("R9", "Liq", "S_PO4"): self.i_PBM
            - self.f_XI * self.i_PXI
            - (1 - self.f_XI) * self.i_PXS,
            ("R9", "Liq", "S_ALK"): (
                self.i_NBM - self.f_XI * self.i_NXI - (1 - self.f_XI) * self.i_NXS
            )
            * mw_alk
            / mw_N
            - (self.i_PBM - self.f_XI * self.i_PXI - (1 - self.f_XI) * self.i_PXS)
            * 1.5
            * mw_alk
            / mw_P,
            ("R9", "Liq", "X_AUT"): 0,
            ("R9", "Liq", "X_H"): -1,
            ("R9", "Liq", "X_I"): self.f_XI,
            ("R9", "Liq", "X_MeOH"): 0,
            ("R9", "Liq", "X_MeP"): 0,
            ("R9", "Liq", "X_PAO"): 0,
            ("R9", "Liq", "X_PHA"): 0,
            ("R9", "Liq", "X_PP"): 0,
            ("R9", "Liq", "X_S"): 1 - self.f_XI,
            ("R9", "Liq", "X_TSS"): -(
                self.i_TSSBM - self.f_XI * self.i_TSSXI - (1 - self.f_XI) * self.i_TSSXS
            ),
            # R10: Storage of X_PHA
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_A"): -1,
            ("R10", "Liq", "S_F"): 0,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "S_N2"): 0,
            ("R10", "Liq", "S_NH4"): 0,
            ("R10", "Liq", "S_NO3"): 0,
            ("R10", "Liq", "S_O2"): 0,
            ("R10", "Liq", "S_PO4"): self.Y_PO4,
            ("R10", "Liq", "S_ALK"): mw_alk / (64 * pyo.units.kg / pyo.units.kmol)
            - self.Y_PO4 * 0.5 * mw_alk / mw_P,
            ("R10", "Liq", "X_AUT"): 0,
            ("R10", "Liq", "X_H"): 0,
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "X_MeOH"): 0,
            ("R10", "Liq", "X_MeP"): 0,
            ("R10", "Liq", "X_PAO"): 0,
            ("R10", "Liq", "X_PHA"): 1,
            ("R10", "Liq", "X_PP"): -self.Y_PO4,
            ("R10", "Liq", "X_S"): 0,
            ("R10", "Liq", "X_TSS"): 0.6 - self.Y_PO4 * 3.23,
            # R11: Aerobic storage of X_PP
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_A"): 0,
            ("R11", "Liq", "S_F"): 0,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "S_N2"): 0,
            ("R11", "Liq", "S_NH4"): 0,
            ("R11", "Liq", "S_NO3"): 0,
            ("R11", "Liq", "S_O2"): -self.Y_PHA,
            ("R11", "Liq", "S_PO4"): -1,
            ("R11", "Liq", "S_ALK"): 0.5 * mw_alk / mw_P,
            ("R11", "Liq", "X_AUT"): 0,
            ("R11", "Liq", "X_H"): 0,
            ("R11", "Liq", "X_I"): 0,
            ("R11", "Liq", "X_MeOH"): 0,
            ("R11", "Liq", "X_MeP"): 0,
            ("R11", "Liq", "X_PAO"): 0,
            ("R11", "Liq", "X_PHA"): -self.Y_PHA,
            ("R11", "Liq", "X_PP"): 1,
            ("R11", "Liq", "X_S"): 0,
            ("R11", "Liq", "X_TSS"): -(0.6 * self.Y_PHA - 3.23),
            # R12: Anoxic storage of X_PP
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_A"): 0,
            ("R12", "Liq", "S_F"): 0,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "S_N2"): self.Y_PHA
            * mw_N
            / (40 * pyo.units.kg / pyo.units.kmol),
            ("R12", "Liq", "S_NH4"): 0,
            ("R12", "Liq", "S_NO3"): -self.Y_PHA
            * mw_N
            / (40 * pyo.units.kg / pyo.units.kmol),
            ("R12", "Liq", "S_O2"): 0,
            ("R12", "Liq", "S_PO4"): -1,
            ("R12", "Liq", "S_ALK"): 0.5 * mw_alk / mw_P
            + self.Y_PHA / (40 * pyo.units.kg / pyo.units.kmol) * mw_alk,
            ("R12", "Liq", "X_AUT"): 0,
            ("R12", "Liq", "X_H"): 0,
            ("R12", "Liq", "X_I"): 0,
            ("R12", "Liq", "X_MeOH"): 0,
            ("R12", "Liq", "X_MeP"): 0,
            ("R12", "Liq", "X_PAO"): 0,
            ("R12", "Liq", "X_PHA"): -self.Y_PHA,
            ("R12", "Liq", "X_PP"): 1,
            ("R12", "Liq", "X_S"): 0,
            ("R12", "Liq", "X_TSS"): -(0.6 * self.Y_PHA - 3.23),
            # R13: Aerobic growth of X_PAO
            ("R13", "Liq", "H2O"): 0,
            ("R13", "Liq", "S_A"): 0,
            ("R13", "Liq", "S_F"): 0,
            ("R13", "Liq", "S_I"): 0,
            ("R13", "Liq", "S_N2"): 0,
            ("R13", "Liq", "S_NH4"): -self.i_NBM,
            ("R13", "Liq", "S_NO3"): 0,
            ("R13", "Liq", "S_O2"): -(1 / self.Y_H - 1),
            ("R13", "Liq", "S_PO4"): -self.i_PBM,
            ("R13", "Liq", "S_ALK"): -self.i_NBM * mw_alk / mw_N
            + self.i_PBM * 1.5 * mw_alk / mw_P,
            ("R13", "Liq", "X_AUT"): 0,
            ("R13", "Liq", "X_H"): 0,
            ("R13", "Liq", "X_I"): 0,
            ("R13", "Liq", "X_MeOH"): 0,
            ("R13", "Liq", "X_MeP"): 0,
            ("R13", "Liq", "X_PAO"): 1,
            ("R13", "Liq", "X_PHA"): -1 / self.Y_H,
            ("R13", "Liq", "X_PP"): 0,
            ("R13", "Liq", "X_S"): 0,
            ("R13", "Liq", "X_TSS"): (self.i_TSSBM - 0.6 / self.Y_H),
            # R14: Anoxic growth of X_PAO
            ("R14", "Liq", "H2O"): 0,
            ("R14", "Liq", "S_A"): 0,
            ("R14", "Liq", "S_F"): 0,
            ("R14", "Liq", "S_I"): 0,
            ("R14", "Liq", "S_N2"): -(1 - 1 / self.Y_H)
            * mw_N
            / (40 * pyo.units.kg / pyo.units.kmol),
            ("R14", "Liq", "S_NH4"): -self.i_NBM,
            ("R14", "Liq", "S_NO3"): (1 - 1 / self.Y_H)
            * mw_N
            / (40 * pyo.units.kg / pyo.units.kmol),
            ("R14", "Liq", "S_O2"): 0,
            ("R14", "Liq", "S_PO4"): -self.i_PBM,
            ("R14", "Liq", "S_ALK"): (
                -self.i_NBM
                - (1 - 1 / self.Y_H) * mw_N / (40 * pyo.units.kg / pyo.units.kmol)
            )
            * mw_alk
            / mw_N
            + self.i_PBM * mw_alk / mw_P,
            ("R14", "Liq", "X_AUT"): 0,
            ("R14", "Liq", "X_H"): 0,
            ("R14", "Liq", "X_I"): 0,
            ("R14", "Liq", "X_MeOH"): 0,
            ("R14", "Liq", "X_MeP"): 0,
            ("R14", "Liq", "X_PAO"): 1,
            ("R14", "Liq", "X_PHA"): -1 / self.Y_H,
            ("R14", "Liq", "X_PP"): 0,
            ("R14", "Liq", "X_S"): 0,
            ("R14", "Liq", "X_TSS"): self.i_TSSBM - 0.6 / self.Y_H,
            # R15: Lysis of X_PAO
            ("R15", "Liq", "H2O"): 0,
            ("R15", "Liq", "S_A"): 0,
            ("R15", "Liq", "S_F"): 0,
            ("R15", "Liq", "S_I"): 0,
            ("R15", "Liq", "S_N2"): 0,
            ("R15", "Liq", "S_NH4"): self.i_NBM
            - self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS,
            ("R15", "Liq", "S_NO3"): 0,
            ("R15", "Liq", "S_O2"): 0,
            ("R15", "Liq", "S_PO4"): self.i_PBM
            - self.f_XI * self.i_PXI
            - (1 - self.f_XI) * self.i_PXS,
            ("R15", "Liq", "S_ALK"): -(
                -self.i_NBM + self.f_XI * self.i_NXI + (1 - self.f_XI) * self.i_NXS
            )
            * mw_alk
            / mw_N
            + (-self.i_PBM + self.f_XI * self.i_PXI + (1 - self.f_XI) * self.i_PXS)
            * 1.5
            * mw_alk
            / mw_P,
            ("R15", "Liq", "X_AUT"): 0,
            ("R15", "Liq", "X_H"): 0,
            ("R15", "Liq", "X_I"): self.f_XI,
            ("R15", "Liq", "X_MeOH"): 0,
            ("R15", "Liq", "X_MeP"): 0,
            ("R15", "Liq", "X_PAO"): -1,
            ("R15", "Liq", "X_PHA"): 0,
            ("R15", "Liq", "X_PP"): 0,
            ("R15", "Liq", "X_S"): 1 - self.f_XI,
            ("R15", "Liq", "X_TSS"): -self.i_TSSBM
            + self.f_XI * self.i_TSSXI
            + (1 - self.f_XI) * self.i_TSSXS,
            # R16: Lysis of X_PP
            ("R16", "Liq", "H2O"): 0,
            ("R16", "Liq", "S_A"): 0,
            ("R16", "Liq", "S_F"): 0,
            ("R16", "Liq", "S_I"): 0,
            ("R16", "Liq", "S_N2"): 0,
            ("R16", "Liq", "S_NH4"): 0,
            ("R16", "Liq", "S_NO3"): 0,
            ("R16", "Liq", "S_O2"): 0,
            ("R16", "Liq", "S_PO4"): 1,
            ("R16", "Liq", "S_ALK"): -0.5 * mw_alk / mw_P,
            ("R16", "Liq", "X_AUT"): 0,
            ("R16", "Liq", "X_H"): 0,
            ("R16", "Liq", "X_I"): 0,
            ("R16", "Liq", "X_MeOH"): 0,
            ("R16", "Liq", "X_MeP"): 0,
            ("R16", "Liq", "X_PAO"): 0,
            ("R16", "Liq", "X_PHA"): 0,
            ("R16", "Liq", "X_PP"): -1,
            ("R16", "Liq", "X_S"): 0,
            ("R16", "Liq", "X_TSS"): -3.23,
            # R17: Lysis of X_PAH
            ("R17", "Liq", "H2O"): 0,
            ("R17", "Liq", "S_A"): 1,
            ("R17", "Liq", "S_F"): 0,
            ("R17", "Liq", "S_I"): 0,
            ("R17", "Liq", "S_N2"): 0,
            ("R17", "Liq", "S_NH4"): 0,
            ("R17", "Liq", "S_NO3"): 0,
            ("R17", "Liq", "S_O2"): 0,
            ("R17", "Liq", "S_PO4"): 0,
            ("R17", "Liq", "S_ALK"): mw_alk / (64 * pyo.units.kg / pyo.units.kmol)
            - mw_alk / (31 * pyo.units.kg / pyo.units.kmol),
            ("R17", "Liq", "X_AUT"): 0,
            ("R17", "Liq", "X_H"): 0,
            ("R17", "Liq", "X_I"): 0,
            ("R17", "Liq", "X_MeOH"): 0,
            ("R17", "Liq", "X_MeP"): 0,
            ("R17", "Liq", "X_PAO"): 0,
            ("R17", "Liq", "X_PHA"): -1,
            ("R17", "Liq", "X_PP"): 0,
            ("R17", "Liq", "X_S"): 0,
            ("R17", "Liq", "X_TSS"): -0.6,
            # R18: Aerobic growth of X_AUT
            ("R18", "Liq", "H2O"): 0,
            ("R18", "Liq", "S_A"): 0,
            ("R18", "Liq", "S_F"): 0,
            ("R18", "Liq", "S_I"): 0,
            ("R18", "Liq", "S_N2"): 0,
            ("R18", "Liq", "S_NH4"): -1 / self.Y_A - self.i_NBM,
            ("R18", "Liq", "S_NO3"): 1 / self.Y_A,
            ("R18", "Liq", "S_O2"): -(4.57 - self.Y_A) / self.Y_A,
            ("R18", "Liq", "S_PO4"): -self.i_PBM,
            ("R18", "Liq", "S_ALK"): (-1 / self.Y_A - self.i_NBM - 1 / self.Y_A)
            * mw_alk
            / mw_N
            + self.i_PBM * 1.5 * mw_alk / mw_P,
            ("R18", "Liq", "X_AUT"): 1,
            ("R18", "Liq", "X_H"): 0,
            ("R18", "Liq", "X_I"): 0,
            ("R18", "Liq", "X_MeOH"): 0,
            ("R18", "Liq", "X_MeP"): 0,
            ("R18", "Liq", "X_PAO"): 0,
            ("R18", "Liq", "X_PHA"): 0,
            ("R18", "Liq", "X_PP"): 0,
            ("R18", "Liq", "X_S"): 0,
            ("R18", "Liq", "X_TSS"): self.i_TSSBM,
            # R19: Lysis of X_AUT
            ("R19", "Liq", "H2O"): 0,
            ("R19", "Liq", "S_A"): 0,
            ("R19", "Liq", "S_F"): 0,
            ("R19", "Liq", "S_I"): 0,
            ("R19", "Liq", "S_N2"): 0,
            ("R19", "Liq", "S_NH4"): -self.f_XI * self.i_NXI
            - (1 - self.f_XI) * self.i_NXS
            + self.i_NBM,
            ("R19", "Liq", "S_NO3"): 0,
            ("R19", "Liq", "S_O2"): 0,
            ("R19", "Liq", "S_PO4"): -self.f_XI * self.i_PXI
            - (1 - self.f_XI) * self.i_PXS
            + self.i_PBM,
            ("R19", "Liq", "S_ALK"): (
                -self.f_XI * self.i_NXI - (1 - self.f_XI) * self.i_NXS + self.i_NBM
            )
            * mw_alk
            / mw_N
            - (-self.f_XI * self.i_PXI - (1 - self.f_XI) * self.i_PXS + self.i_PBM)
            * mw_alk
            / mw_P,
            ("R19", "Liq", "X_AUT"): -1,
            ("R19", "Liq", "X_H"): 0,
            ("R19", "Liq", "X_I"): self.f_XI,
            ("R19", "Liq", "X_MeOH"): 0,
            ("R19", "Liq", "X_MeP"): 0,
            ("R19", "Liq", "X_PAO"): 0,
            ("R19", "Liq", "X_PHA"): 0,
            ("R19", "Liq", "X_PP"): 0,
            ("R19", "Liq", "X_S"): 1 - self.f_XI,
            ("R19", "Liq", "X_TSS"): self.f_XI * self.i_TSSXI
            + (1 - self.f_XI) * self.i_TSSXS
            - self.i_TSSBM,
            # R20: Precipitation
            ("R20", "Liq", "H2O"): 0,
            ("R20", "Liq", "S_A"): 0,
            ("R20", "Liq", "S_F"): 0,
            ("R20", "Liq", "S_I"): 0,
            ("R20", "Liq", "S_N2"): 0,
            ("R20", "Liq", "S_NH4"): 0,
            ("R20", "Liq", "S_NO3"): 0,
            ("R20", "Liq", "S_O2"): 0,
            ("R20", "Liq", "S_PO4"): -1,
            ("R20", "Liq", "S_ALK"): 1.5 * mw_alk / mw_P,
            ("R20", "Liq", "X_AUT"): 0,
            ("R20", "Liq", "X_H"): 0,
            ("R20", "Liq", "X_I"): 0,
            ("R20", "Liq", "X_MeOH"): -3.45,
            ("R20", "Liq", "X_MeP"): 4.87,
            ("R20", "Liq", "X_PAO"): 0,
            ("R20", "Liq", "X_PHA"): 0,
            ("R20", "Liq", "X_PP"): 0,
            ("R20", "Liq", "X_S"): 0,
            ("R20", "Liq", "X_TSS"): 1.42,
            # R21: Re-dissolution
            ("R21", "Liq", "H2O"): 0,
            ("R21", "Liq", "S_A"): 0,
            ("R21", "Liq", "S_F"): 0,
            ("R21", "Liq", "S_I"): 0,
            ("R21", "Liq", "S_N2"): 0,
            ("R21", "Liq", "S_NH4"): 0,
            ("R21", "Liq", "S_NO3"): 0,
            ("R21", "Liq", "S_O2"): 0,
            ("R21", "Liq", "S_PO4"): 1,
            ("R21", "Liq", "S_ALK"): -1.5 * mw_alk / mw_P,
            ("R21", "Liq", "X_AUT"): 0,
            ("R21", "Liq", "X_H"): 0,
            ("R21", "Liq", "X_I"): 0,
            ("R21", "Liq", "X_MeOH"): 3.45,
            ("R21", "Liq", "X_MeP"): -4.87,
            ("R21", "Liq", "X_PAO"): 0,
            ("R21", "Liq", "X_PHA"): 0,
            ("R21", "Liq", "X_PP"): 0,
            ("R21", "Liq", "X_S"): 0,
            ("R21", "Liq", "X_TSS"): -1.42,
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


class _ASM2dReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    def initialize(blk, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info("Initialization Complete.")


@declare_process_block_class("ASM2dReactionBlock", block_class=_ASM2dReactionBlock)
class ASM2dReactionBlockData(ReactionBlockDataBase):
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

            def rate_expression_rule(b, r):
                if r == "R1":
                    # R1: Aerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            (b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"])
                            / (
                                b.params.K_X
                                + (
                                    b.conc_mass_comp_ref["X_S"]
                                    / b.conc_mass_comp_ref["X_H"]
                                )
                            )
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.eta_NO3
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            (b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"])
                            / (
                                b.params.K_X
                                + (
                                    b.conc_mass_comp_ref["X_S"]
                                    / b.conc_mass_comp_ref["X_H"]
                                )
                            )
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Anaerobic hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.K_H
                        * b.params.eta_fe
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.K_NO3
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            (b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"])
                            / (
                                b.params.K_X
                                + (
                                    b.conc_mass_comp_ref["X_S"]
                                    / b.conc_mass_comp_ref["X_H"]
                                )
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
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
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
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
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
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.K_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Anoxic growth on S_F
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.eta_NO3
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
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
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Anoxic growth on S_A, denitrification
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * b.params.eta_NO3
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (b.params.K_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_A"]
                            / (
                                b.conc_mass_comp_ref["S_F"]
                                + b.conc_mass_comp_ref["S_A"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Fermentation
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.q_fe
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.params.K_NO3
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_F"]
                            / (b.params.K_F + b.conc_mass_comp_ref["S_F"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Lysis
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
                            / (b.params.K_A + b.conc_mass_comp_ref["S_A"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["X_PP"]
                                / b.conc_mass_comp_ref["X_PAO"]
                            )
                            / (
                                b.params.K_PP
                                + (
                                    b.conc_mass_comp_ref["X_PP"]
                                    / b.conc_mass_comp_ref["X_PAO"]
                                )
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
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_PS + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["X_PHA"]
                                / b.conc_mass_comp_ref["X_PAO"]
                            )
                            / (
                                b.params.K_PHA
                                + (
                                    b.conc_mass_comp_ref["X_PHA"]
                                    / b.conc_mass_comp_ref["X_PAO"]
                                )
                            )
                        )
                        * (
                            (
                                b.params.K_MAX
                                - (
                                    b.conc_mass_comp_ref["X_PP"]
                                    / b.conc_mass_comp_ref["X_PAO"]
                                )
                            )
                            / (
                                b.params.K_IPP
                                + b.params.K_MAX
                                - (
                                    b.conc_mass_comp_ref["X_PP"]
                                    / b.conc_mass_comp_ref["X_PAO"]
                                )
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Anoxic storage of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.reaction_rate["R11"]
                        * b.params.eta_NO3
                        * (b.params.K_O2 / b.conc_mass_comp_ref["S_O2"])
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Aerobic growth of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_PAO
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["X_PHA"]
                                / b.conc_mass_comp_ref["X_PAO"]
                            )
                            / (
                                b.params.K_PHA
                                + (
                                    b.conc_mass_comp_ref["X_PHA"]
                                    / b.conc_mass_comp_ref["X_PAO"]
                                )
                            )
                        )
                        * b.conc_mass_comp_ref["X_PAO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Anoxic growth of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.reaction_rate["R13"]
                        * b.params.eta_NO3
                        * (b.params.K_O2 / b.conc_mass_comp_ref["S_O2"])
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Lysis of X_PAO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PAO
                        * b.conc_mass_comp_ref["X_PAO"]
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R16":
                    # R16: Lysis of X_PP
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PP
                        * b.conc_mass_comp_ref["X_PP"]
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R17":
                    # R17: Lysis of X_PAH
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_PHA
                        * b.conc_mass_comp_ref["X_PHA"]
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R18":
                    # R18: Aerobic growth of X_AUT
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_AUT * b.conc_mass_comp_ref["X_AUT"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R19":
                    # R19: Lysis of X_AUT
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_AUT
                        * (
                            b.conc_mass_comp_ref["S_O2"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O2"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_PO4"]
                            / (b.params.K_P + b.conc_mass_comp_ref["S_PO4"])
                        )
                        * b.conc_mass_comp_ref["X_AUT"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R20":
                    # R20: Precipitation
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_pre
                        * b.conc_mass_comp_ref["S_PO4"]
                        * b.conc_mass_comp_ref["X_MeOH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R21":
                    # R21: Re-dissolution
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_red
                        * b.conc_mass_comp_ref["X_MeP"]
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        ),
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

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.rate_expression.items():
            # TODO: Need to work out how to calculate good scaling factors
            # instead of a fixed 1e3.
            iscale.constraint_scaling_transform(c, 1e3, overwrite=True)
