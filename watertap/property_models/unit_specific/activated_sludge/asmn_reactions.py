#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
ASMN reaction package.

References:

[1] Hiatt, W.C., Grady, C.P.L. Jr., "An Updated Process Model for Carbon
Oxidation, Nitrification, and Denitrification", Water Environ. Res., 80, 2145
(2008).

[2] Snip, L.J.P., "Quantifying the greenhouse gas emissions of waste water
treatment plants", MSc thesis, Wageningen University (2010), Tables 5.1-5.2.

Note: S_N2 is not among the 20 components of [1]; it is added as the product
of N2O reduction (R5) following [2] so that N and COD balances close.
"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, Bool

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
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
import idaes.core.util.scaling as iscale

# Some more information about this module
__author__ = "Chenyu Wang"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASMNReactionParameterBlock")
class ASMNReactionParameterData(ReactionParameterBlock):
    """
    Reaction Parameter Block Class
    """

    CONFIG = ReactionParameterBlock.CONFIG()
    CONFIG.declare(
        "include_salt_inhibition",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include salt (S_Z) noncompetitive inhibition in mu_H, mu_A1, mu_A2",
        ),
    )
    CONFIG.declare(
        "include_aob_inhibitor",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include AOB-specific organic inhibitor (S_J) noncompetitive inhibition in mu_A1",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ASMNReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ASMN paper
        # R1: Aerobic growth of heterotrophs
        # R2: Anoxic growth of heterotrophs, reducing nitrate to nitrite
        # R3: Anoxic growth of heterotrophs, reducing nitrite to nitric oxide
        # R4: Anoxic growth of heterotrophs, reducing nitric oxide to nitrous oxide
        # R5: Anoxic growth of heterotrophs, reducing nitrous oxide to nitrogen
        # R6: ANRA, nitrate reduction to nitrite
        # R7: ANRA, nitrite reduction to ammonia
        # R8: Death and lysis of heterotrophs
        # R9: Autotrophic growth of AOB
        # R10: Autotrophic growth of NOB
        # R11: Mixotrophic growth of NOB
        # R12: Death and lysis of AOB
        # R13: Death and lysis of NOB
        # R14: Ammonification of soluble organic nitrogen
        # R15: Hydrolysis of particulate organics
        # R16: Hydrolysis of particulate organic nitrogen
        # R17: Biodegradation of an inhibitory compound
        # R18: Biodegradation of a special interest compound

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
            ]
        )

        # Stoichiometric Parameters
        self.Y_H = pyo.Var(
            initialize=0.6,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Heterotrophic yield, Y_H",
        )
        self.Y_A1 = pyo.Var(
            initialize=0.18,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Autotrophic yield, AOB, Y_A1",
        )
        self.Y_A2 = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Autotrophic yield, NOB, Y_A2",
        )
        self.eta_Y = pyo.Var(
            initialize=0.9,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic yield factor, eta_Y",
        )

        add_object_reference(self, "f_D", self.config.property_package.f_D)
        add_object_reference(self, "i_NXB", self.config.property_package.i_NXB)
        add_object_reference(self, "i_NXD", self.config.property_package.i_NXD)

        # Kinetic Parameters - Heterotrophs
        self.mu_H = pyo.Var(
            initialize=6.25,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for heterotrophic biomass, mu_H",
        )
        self.b_LH = pyo.Var(
            initialize=0.408,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient for heterotrophic biomass, b_L,H",
        )
        self.eta_g = pyo.Var(
            initialize=0.8,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic growth factor, R17, R18",
        )
        self.eta_g2 = pyo.Var(
            initialize=0.28,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic growth factor, R2",
        )
        self.eta_g3 = pyo.Var(
            initialize=0.16,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic growth factor, R3",
        )
        self.eta_g4 = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic growth factor, R4",
        )
        self.eta_g5 = pyo.Var(
            initialize=0.35,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic growth factor, R5",
        )
        self.eta_h = pyo.Var(
            initialize=0.4,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic hydrolysis factor",
        )

        # Half-saturation coefficients for substrate
        self.K_S1 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, R1, [mg COD/L]",
        )
        self.K_S2 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, R2, [mg COD/L]",
        )
        self.K_S3 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, R3, [mg COD/L]",
        )
        self.K_S4 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, R4, [mg COD/L]",
        )
        self.K_S5 = pyo.Var(
            initialize=40e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, R5, [mg COD/L]",
        )

        # Half-saturation coefficients for oxygen, heterotrophs
        self.K_OH1 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, heterotrophs, R1, [mg O2/L]",
        )
        self.K_OH2 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, heterotrophs, R2, [mg O2/L]",
        )
        self.K_OH3 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, heterotrophs, R3, [mg O2/L]",
        )
        self.K_OH4 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, heterotrophs, R4, [mg O2/L]",
        )
        self.K_OH5 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, heterotrophs, R5, [mg O2/L]",
        )

        # Half-saturation coefficients for nitrogen oxides, heterotrophs
        self.K_NO3 = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitrate-nitrogen, heterotrophs, [mg N/L]",
        )
        self.K_NO2 = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitrite-nitrogen, heterotrophs, [mg N/L]",
        )
        self.K_NO = pyo.Var(
            initialize=0.05e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitric oxide-nitrogen, heterotrophs, [mg N/L]",
        )
        self.K_N2O = pyo.Var(
            initialize=0.05e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitrous oxide-nitrogen, heterotrophs, [mg N/L]",
        )

        # NO inhibition coefficients
        self.K_I3NO = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitric oxide inhibition coefficient, R3, [mg N/L]",
        )
        self.K_I4NO = pyo.Var(
            initialize=0.3e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitric oxide inhibition coefficient, R4, [mg N/L]",
        )
        self.K_I5NO = pyo.Var(
            initialize=0.075e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitric oxide inhibition coefficient, R5, [mg N/L]",
        )

        # ANRA parameters
        self.K_6NO3 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitrate-nitrogen, ANRA R6, [mg N/L]",
        )
        self.K_I6NO2 = pyo.Var(
            initialize=0.05e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitrite inhibition coefficient, ANRA R6, [mg N/L]",
        )
        self.K_I6NH = pyo.Var(
            initialize=0.05e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Ammonia inhibition coefficient, ANRA R6, [mg N/L]",
        )
        self.K_7NO2 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for nitrite-nitrogen, ANRA R7, [mg N/L]",
        )
        self.K_I7NH = pyo.Var(
            initialize=0.05e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Ammonia inhibition coefficient, ANRA R7, [mg N/L]",
        )

        # Heterotroph ammonia limit
        self.K_N1 = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for NH3, heterotrophs, [mg N/L]",
        )

        # Hydrolysis
        self.k_h = pyo.Var(
            initialize=2.208,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Hydrolysis coefficient, k_h, [mg COD/ mg biomass COD * day]",
        )
        self.K_X = pyo.Var(
            initialize=0.15,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for hydrolysis of slowly biodegradable substrate, K_X",
        )

        # Ammonification
        self.k_a = pyo.Var(
            initialize=0.1608,
            units=pyo.units.L / pyo.units.mg / pyo.units.hr,
            domain=pyo.PositiveReals,
            doc="Ammonification rate coefficient, k_a",
        )

        # Kinetic Parameters - Autotrophs
        self.mu_A1 = pyo.Var(
            initialize=0.78,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for AOB, mu_A1",
        )
        self.mu_A2 = pyo.Var(
            initialize=0.78,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for NOB (autotrophic), mu_A2",
        )
        self.b_LA1 = pyo.Var(
            initialize=0.096,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient, AOB, b_L,A1",
        )
        self.b_LA2 = pyo.Var(
            initialize=0.096,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient, NOB, b_L,A2",
        )
        self.beta_g = pyo.Var(
            initialize=6,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="NOB mixotrophic growth factor, beta_g",
        )
        self.K_FA = pyo.Var(
            initialize=0.0075e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for free ammonia, [mg/L]",
        )
        self.K_FNA = pyo.Var(
            initialize=0.0001e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for free nitrous acid, [mg/L]",
        )
        self.K_S11 = pyo.Var(
            initialize=20e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for substrate, S_S, R11, [mg COD/L]",
        )
        self.K_OA1 = pyo.Var(
            initialize=0.6e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, AOB, [mg O2/L]",
        )
        self.K_OA2 = pyo.Var(
            initialize=1.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for O2, NOB, [mg O2/L]",
        )
        self.K_I9FA = pyo.Var(
            initialize=1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Free ammonia inhibition coefficient, R9, [mg N/L]",
        )
        self.K_I9FNA = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Free nitrous acid inhibition coefficient, R9, [mg N/L]",
        )
        self.K_I10FA = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Free ammonia inhibition coefficient, R10, [mg N/L]",
        )
        self.K_I10FNA = pyo.Var(
            initialize=0.04e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Free nitrous acid inhibition coefficient, R10, [mg N/L]",
        )

        # Inhibitor and priority pollutant kinetics
        self.K_SJ = pyo.Var(
            initialize=2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for AOB inhibitor, S_J",
        )
        self.K_SW = pyo.Var(
            initialize=2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for priority pollutant, S_W",
        )

        # Optional salt inhibition coefficients (only created if user opts in)
        # User must fix these to site-specific values; no defaults are given in Hiatt and Grady (2008).
        if self.config.include_salt_inhibition:
            self.K_IHZ = pyo.Var(
                units=pyo.units.kg / pyo.units.m**3,
                domain=pyo.PositiveReals,
                doc="Heterotrophic salt inhibition coefficient",
            )
            self.K_IA1Z = pyo.Var(
                units=pyo.units.kg / pyo.units.m**3,
                domain=pyo.PositiveReals,
                doc="AOB salt inhibition coefficient",
            )
            self.K_IA2Z = pyo.Var(
                units=pyo.units.kg / pyo.units.m**3,
                domain=pyo.PositiveReals,
                doc="NOB salt inhibition coefficient",
            )

        # Optional AOB-specific organic inhibitor coefficient
        if self.config.include_aob_inhibitor:
            self.K_ISJ = pyo.Var(
                units=pyo.units.kg / pyo.units.m**3,
                domain=pyo.PositiveReals,
                doc="AOB inhibitor inhibition coefficient",
            )

        # Reaction Stoichiometry
        # Note that reaction stoichiometry is on a mass basis.
        # For alkalinity, this requires converting the mass of nitrogen species
        # reacted to mass of alkalinity converted using a charge balance (effectively MW_C/MW_N)
        mw_alk = 12 * pyo.units.kg / pyo.units.kmol
        mw_n = 14 * pyo.units.kg / pyo.units.kmol

        # Electron equivalents [g COD/g N] as exact fractions
        # (Table 4 of [1] rounds to 1.143, 0.571, 1.14, 3.43);
        cod_2e = 8 / 7  # 2 e-/N: NO3- -> NO2-, NO2- -> NO3-
        cod_1e = 4 / 7  # 1 e-/N: NO2- -> NO, NO -> N2O, N2O -> N2
        cod_6e = 24 / 7  # 6 e-/N: NO2- -> NH4+, NH4+ -> NO2-

        # Helper denitrification stoichiometric coefficients
        A = (1 - self.Y_H * self.eta_Y) / (cod_2e * self.Y_H * self.eta_Y)
        B = (1 - self.Y_H * self.eta_Y) / (cod_1e * self.Y_H * self.eta_Y)

        # Component order follows Table 4 of Hiatt and Grady (2008):
        # X_S, X_BH, X_BA1, X_BA2, X_D, S_S, S_O, S_NO3, S_NO2, S_NO, S_N2O,
        # S_NH, S_NS, X_NS, S_ALK, S_J, S_W, then inert components (H2O, S_I, X_I, S_Z)
        # S_N2 (after S_N2O) is not in [1]; produced only by R5, per [2]
        self.rate_reaction_stoichiometry = {
            # R1: Aerobic growth of heterotrophs
            ("R1", "Liq", "X_S"): 0,
            ("R1", "Liq", "X_BH"): 1,
            ("R1", "Liq", "X_BA1"): 0,
            ("R1", "Liq", "X_BA2"): 0,
            ("R1", "Liq", "X_D"): 0,
            ("R1", "Liq", "S_S"): -1 / self.Y_H,
            ("R1", "Liq", "S_O"): -(1 - self.Y_H) / self.Y_H,
            ("R1", "Liq", "S_NO3"): 0,
            ("R1", "Liq", "S_NO2"): 0,
            ("R1", "Liq", "S_NO"): 0,
            ("R1", "Liq", "S_N2O"): 0,
            ("R1", "Liq", "S_N2"): 0,
            ("R1", "Liq", "S_NH"): -self.i_NXB,
            ("R1", "Liq", "S_NS"): 0,
            ("R1", "Liq", "X_NS"): 0,
            ("R1", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R1", "Liq", "S_J"): 0,
            ("R1", "Liq", "S_W"): 0,
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_I"): 0,
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "S_Z"): 0,
            # R2: Anoxic growth of heterotrophs, NO3 -> NO2
            ("R2", "Liq", "X_S"): 0,
            ("R2", "Liq", "X_BH"): 1,
            ("R2", "Liq", "X_BA1"): 0,
            ("R2", "Liq", "X_BA2"): 0,
            ("R2", "Liq", "X_D"): 0,
            ("R2", "Liq", "S_S"): -1 / (self.Y_H * self.eta_Y),
            ("R2", "Liq", "S_O"): 0,
            ("R2", "Liq", "S_NO3"): -A,
            ("R2", "Liq", "S_NO2"): A,
            ("R2", "Liq", "S_NO"): 0,
            ("R2", "Liq", "S_N2O"): 0,
            ("R2", "Liq", "S_N2"): 0,
            ("R2", "Liq", "S_NH"): -self.i_NXB,
            ("R2", "Liq", "S_NS"): 0,
            ("R2", "Liq", "X_NS"): 0,
            ("R2", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R2", "Liq", "S_J"): 0,
            ("R2", "Liq", "S_W"): 0,
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_I"): 0,
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "S_Z"): 0,
            # R3: Anoxic growth of heterotrophs, NO2 -> NO (alkalinity generated here)
            ("R3", "Liq", "X_S"): 0,
            ("R3", "Liq", "X_BH"): 1,
            ("R3", "Liq", "X_BA1"): 0,
            ("R3", "Liq", "X_BA2"): 0,
            ("R3", "Liq", "X_D"): 0,
            ("R3", "Liq", "S_S"): -1 / (self.Y_H * self.eta_Y),
            ("R3", "Liq", "S_O"): 0,
            ("R3", "Liq", "S_NO3"): 0,
            ("R3", "Liq", "S_NO2"): -B,
            ("R3", "Liq", "S_NO"): B,
            ("R3", "Liq", "S_N2O"): 0,
            ("R3", "Liq", "S_N2"): 0,
            ("R3", "Liq", "S_NH"): -self.i_NXB,
            ("R3", "Liq", "S_NS"): 0,
            ("R3", "Liq", "X_NS"): 0,
            ("R3", "Liq", "S_ALK"): (B - self.i_NXB) * (mw_alk / mw_n),
            ("R3", "Liq", "S_J"): 0,
            ("R3", "Liq", "S_W"): 0,
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_I"): 0,
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "S_Z"): 0,
            # R4: Anoxic growth of heterotrophs, NO -> N2O
            ("R4", "Liq", "X_S"): 0,
            ("R4", "Liq", "X_BH"): 1,
            ("R4", "Liq", "X_BA1"): 0,
            ("R4", "Liq", "X_BA2"): 0,
            ("R4", "Liq", "X_D"): 0,
            ("R4", "Liq", "S_S"): -1 / (self.Y_H * self.eta_Y),
            ("R4", "Liq", "S_O"): 0,
            ("R4", "Liq", "S_NO3"): 0,
            ("R4", "Liq", "S_NO2"): 0,
            ("R4", "Liq", "S_NO"): -B,
            ("R4", "Liq", "S_N2O"): B,
            ("R4", "Liq", "S_N2"): 0,
            ("R4", "Liq", "S_NH"): -self.i_NXB,
            ("R4", "Liq", "S_NS"): 0,
            ("R4", "Liq", "X_NS"): 0,
            ("R4", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R4", "Liq", "S_J"): 0,
            ("R4", "Liq", "S_W"): 0,
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "S_Z"): 0,
            # R5: Anoxic growth of heterotrophs, N2O -> N2
            ("R5", "Liq", "X_S"): 0,
            ("R5", "Liq", "X_BH"): 1,
            ("R5", "Liq", "X_BA1"): 0,
            ("R5", "Liq", "X_BA2"): 0,
            ("R5", "Liq", "X_D"): 0,
            ("R5", "Liq", "S_S"): -1 / (self.Y_H * self.eta_Y),
            ("R5", "Liq", "S_O"): 0,
            ("R5", "Liq", "S_NO3"): 0,
            ("R5", "Liq", "S_NO2"): 0,
            ("R5", "Liq", "S_NO"): 0,
            ("R5", "Liq", "S_N2O"): -B,
            ("R5", "Liq", "S_N2"): B,
            ("R5", "Liq", "S_NH"): -self.i_NXB,
            ("R5", "Liq", "S_NS"): 0,
            ("R5", "Liq", "X_NS"): 0,
            ("R5", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R5", "Liq", "S_J"): 0,
            ("R5", "Liq", "S_W"): 0,
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "S_Z"): 0,
            # R6: ANRA, nitrate reduction to nitrite
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_BH"): 0,
            ("R6", "Liq", "X_BA1"): 0,
            ("R6", "Liq", "X_BA2"): 0,
            ("R6", "Liq", "X_D"): 0,
            ("R6", "Liq", "S_S"): -cod_2e,
            ("R6", "Liq", "S_O"): 0,
            ("R6", "Liq", "S_NO3"): -1,
            ("R6", "Liq", "S_NO2"): 1,
            ("R6", "Liq", "S_NO"): 0,
            ("R6", "Liq", "S_N2O"): 0,
            ("R6", "Liq", "S_N2"): 0,
            ("R6", "Liq", "S_NH"): 0,
            ("R6", "Liq", "S_NS"): 0,
            ("R6", "Liq", "X_NS"): 0,
            ("R6", "Liq", "S_ALK"): 0,
            ("R6", "Liq", "S_J"): 0,
            ("R6", "Liq", "S_W"): 0,
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "S_Z"): 0,
            # R7: ANRA, nitrite reduction to ammonia
            ("R7", "Liq", "X_S"): 0,
            ("R7", "Liq", "X_BH"): 0,
            ("R7", "Liq", "X_BA1"): 0,
            ("R7", "Liq", "X_BA2"): 0,
            ("R7", "Liq", "X_D"): 0,
            ("R7", "Liq", "S_S"): -cod_6e,
            ("R7", "Liq", "S_O"): 0,
            ("R7", "Liq", "S_NO3"): 0,
            ("R7", "Liq", "S_NO2"): -1,
            ("R7", "Liq", "S_NO"): 0,
            ("R7", "Liq", "S_N2O"): 0,
            ("R7", "Liq", "S_N2"): 0,
            ("R7", "Liq", "S_NH"): 1,
            ("R7", "Liq", "S_NS"): 0,
            ("R7", "Liq", "X_NS"): 0,
            ("R7", "Liq", "S_ALK"): 2 * (mw_alk / mw_n),
            ("R7", "Liq", "S_J"): 0,
            ("R7", "Liq", "S_W"): 0,
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "S_Z"): 0,
            # R8: Death and lysis of heterotrophs
            ("R8", "Liq", "X_S"): 1 - self.f_D,
            ("R8", "Liq", "X_BH"): -1,
            ("R8", "Liq", "X_BA1"): 0,
            ("R8", "Liq", "X_BA2"): 0,
            ("R8", "Liq", "X_D"): self.f_D,
            ("R8", "Liq", "S_S"): 0,
            ("R8", "Liq", "S_O"): 0,
            ("R8", "Liq", "S_NO3"): 0,
            ("R8", "Liq", "S_NO2"): 0,
            ("R8", "Liq", "S_NO"): 0,
            ("R8", "Liq", "S_N2O"): 0,
            ("R8", "Liq", "S_N2"): 0,
            ("R8", "Liq", "S_NH"): 0,
            ("R8", "Liq", "S_NS"): 0,
            ("R8", "Liq", "X_NS"): self.i_NXB - self.f_D * self.i_NXD,
            ("R8", "Liq", "S_ALK"): 0,
            ("R8", "Liq", "S_J"): 0,
            ("R8", "Liq", "S_W"): 0,
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "S_Z"): 0,
            # R9: Autotrophic growth of AOB
            ("R9", "Liq", "X_S"): 0,
            ("R9", "Liq", "X_BH"): 0,
            ("R9", "Liq", "X_BA1"): 1,
            ("R9", "Liq", "X_BA2"): 0,
            ("R9", "Liq", "X_D"): 0,
            ("R9", "Liq", "S_S"): 0,
            ("R9", "Liq", "S_O"): -(cod_6e - self.Y_A1) / self.Y_A1,
            ("R9", "Liq", "S_NO3"): 0,
            ("R9", "Liq", "S_NO2"): 1 / self.Y_A1,
            ("R9", "Liq", "S_NO"): 0,
            ("R9", "Liq", "S_N2O"): 0,
            ("R9", "Liq", "S_N2"): 0,
            ("R9", "Liq", "S_NH"): -self.i_NXB - 1 / self.Y_A1,
            ("R9", "Liq", "S_NS"): 0,
            ("R9", "Liq", "X_NS"): 0,
            ("R9", "Liq", "S_ALK"): (-self.i_NXB - 2 / self.Y_A1) * (mw_alk / mw_n),
            ("R9", "Liq", "S_J"): 0,
            ("R9", "Liq", "S_W"): 0,
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "X_I"): 0,
            ("R9", "Liq", "S_Z"): 0,
            # R10: Autotrophic growth of NOB
            ("R10", "Liq", "X_S"): 0,
            ("R10", "Liq", "X_BH"): 0,
            ("R10", "Liq", "X_BA1"): 0,
            ("R10", "Liq", "X_BA2"): 1,
            ("R10", "Liq", "X_D"): 0,
            ("R10", "Liq", "S_S"): 0,
            ("R10", "Liq", "S_O"): -(cod_2e - self.Y_A2) / self.Y_A2,
            ("R10", "Liq", "S_NO3"): 1 / self.Y_A2,
            ("R10", "Liq", "S_NO2"): -1 / self.Y_A2,
            ("R10", "Liq", "S_NO"): 0,
            ("R10", "Liq", "S_N2O"): 0,
            ("R10", "Liq", "S_N2"): 0,
            ("R10", "Liq", "S_NH"): -self.i_NXB,
            ("R10", "Liq", "S_NS"): 0,
            ("R10", "Liq", "X_NS"): 0,
            ("R10", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R10", "Liq", "S_J"): 0,
            ("R10", "Liq", "S_W"): 0,
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "S_Z"): 0,
            # R11: Mixotrophic growth of NOB
            ("R11", "Liq", "X_S"): 0,
            ("R11", "Liq", "X_BH"): 0,
            ("R11", "Liq", "X_BA1"): 0,
            ("R11", "Liq", "X_BA2"): 1,
            ("R11", "Liq", "X_D"): 0,
            ("R11", "Liq", "S_S"): -1 / self.Y_H,
            ("R11", "Liq", "S_O"): -(1 - self.Y_H) / self.Y_H,
            ("R11", "Liq", "S_NO3"): 0,
            ("R11", "Liq", "S_NO2"): 0,
            ("R11", "Liq", "S_NO"): 0,
            ("R11", "Liq", "S_N2O"): 0,
            ("R11", "Liq", "S_N2"): 0,
            ("R11", "Liq", "S_NH"): -self.i_NXB,
            ("R11", "Liq", "S_NS"): 0,
            ("R11", "Liq", "X_NS"): 0,
            ("R11", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R11", "Liq", "S_J"): 0,
            ("R11", "Liq", "S_W"): 0,
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "X_I"): 0,
            ("R11", "Liq", "S_Z"): 0,
            # R12: Death and lysis of AOB
            ("R12", "Liq", "X_S"): 1 - self.f_D,
            ("R12", "Liq", "X_BH"): 0,
            ("R12", "Liq", "X_BA1"): -1,
            ("R12", "Liq", "X_BA2"): 0,
            ("R12", "Liq", "X_D"): self.f_D,
            ("R12", "Liq", "S_S"): 0,
            ("R12", "Liq", "S_O"): 0,
            ("R12", "Liq", "S_NO3"): 0,
            ("R12", "Liq", "S_NO2"): 0,
            ("R12", "Liq", "S_NO"): 0,
            ("R12", "Liq", "S_N2O"): 0,
            ("R12", "Liq", "S_N2"): 0,
            ("R12", "Liq", "S_NH"): 0,
            ("R12", "Liq", "S_NS"): 0,
            ("R12", "Liq", "X_NS"): self.i_NXB - self.f_D * self.i_NXD,
            ("R12", "Liq", "S_ALK"): 0,
            ("R12", "Liq", "S_J"): 0,
            ("R12", "Liq", "S_W"): 0,
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "X_I"): 0,
            ("R12", "Liq", "S_Z"): 0,
            # R13: Death and lysis of NOB
            ("R13", "Liq", "X_S"): 1 - self.f_D,
            ("R13", "Liq", "X_BH"): 0,
            ("R13", "Liq", "X_BA1"): 0,
            ("R13", "Liq", "X_BA2"): -1,
            ("R13", "Liq", "X_D"): self.f_D,
            ("R13", "Liq", "S_S"): 0,
            ("R13", "Liq", "S_O"): 0,
            ("R13", "Liq", "S_NO3"): 0,
            ("R13", "Liq", "S_NO2"): 0,
            ("R13", "Liq", "S_NO"): 0,
            ("R13", "Liq", "S_N2O"): 0,
            ("R13", "Liq", "S_N2"): 0,
            ("R13", "Liq", "S_NH"): 0,
            ("R13", "Liq", "S_NS"): 0,
            ("R13", "Liq", "X_NS"): self.i_NXB - self.f_D * self.i_NXD,
            ("R13", "Liq", "S_ALK"): 0,
            ("R13", "Liq", "S_J"): 0,
            ("R13", "Liq", "S_W"): 0,
            ("R13", "Liq", "H2O"): 0,
            ("R13", "Liq", "S_I"): 0,
            ("R13", "Liq", "X_I"): 0,
            ("R13", "Liq", "S_Z"): 0,
            # R14: Ammonification of soluble organic nitrogen
            ("R14", "Liq", "X_S"): 0,
            ("R14", "Liq", "X_BH"): 0,
            ("R14", "Liq", "X_BA1"): 0,
            ("R14", "Liq", "X_BA2"): 0,
            ("R14", "Liq", "X_D"): 0,
            ("R14", "Liq", "S_S"): 0,
            ("R14", "Liq", "S_O"): 0,
            ("R14", "Liq", "S_NO3"): 0,
            ("R14", "Liq", "S_NO2"): 0,
            ("R14", "Liq", "S_NO"): 0,
            ("R14", "Liq", "S_N2O"): 0,
            ("R14", "Liq", "S_N2"): 0,
            ("R14", "Liq", "S_NH"): 1,
            ("R14", "Liq", "S_NS"): -1,
            ("R14", "Liq", "X_NS"): 0,
            ("R14", "Liq", "S_ALK"): 1 * (mw_alk / mw_n),
            ("R14", "Liq", "S_J"): 0,
            ("R14", "Liq", "S_W"): 0,
            ("R14", "Liq", "H2O"): 0,
            ("R14", "Liq", "S_I"): 0,
            ("R14", "Liq", "X_I"): 0,
            ("R14", "Liq", "S_Z"): 0,
            # R15: Hydrolysis of particulate organics
            ("R15", "Liq", "X_S"): -1,
            ("R15", "Liq", "X_BH"): 0,
            ("R15", "Liq", "X_BA1"): 0,
            ("R15", "Liq", "X_BA2"): 0,
            ("R15", "Liq", "X_D"): 0,
            ("R15", "Liq", "S_S"): 1,
            ("R15", "Liq", "S_O"): 0,
            ("R15", "Liq", "S_NO3"): 0,
            ("R15", "Liq", "S_NO2"): 0,
            ("R15", "Liq", "S_NO"): 0,
            ("R15", "Liq", "S_N2O"): 0,
            ("R15", "Liq", "S_N2"): 0,
            ("R15", "Liq", "S_NH"): 0,
            ("R15", "Liq", "S_NS"): 0,
            ("R15", "Liq", "X_NS"): 0,
            ("R15", "Liq", "S_ALK"): 0,
            ("R15", "Liq", "S_J"): 0,
            ("R15", "Liq", "S_W"): 0,
            ("R15", "Liq", "H2O"): 0,
            ("R15", "Liq", "S_I"): 0,
            ("R15", "Liq", "X_I"): 0,
            ("R15", "Liq", "S_Z"): 0,
            # R16: Hydrolysis of particulate organic nitrogen
            ("R16", "Liq", "X_S"): 0,
            ("R16", "Liq", "X_BH"): 0,
            ("R16", "Liq", "X_BA1"): 0,
            ("R16", "Liq", "X_BA2"): 0,
            ("R16", "Liq", "X_D"): 0,
            ("R16", "Liq", "S_S"): 0,
            ("R16", "Liq", "S_O"): 0,
            ("R16", "Liq", "S_NO3"): 0,
            ("R16", "Liq", "S_NO2"): 0,
            ("R16", "Liq", "S_NO"): 0,
            ("R16", "Liq", "S_N2O"): 0,
            ("R16", "Liq", "S_N2"): 0,
            ("R16", "Liq", "S_NH"): 0,
            ("R16", "Liq", "S_NS"): 1,
            ("R16", "Liq", "X_NS"): -1,
            ("R16", "Liq", "S_ALK"): 0,
            ("R16", "Liq", "S_J"): 0,
            ("R16", "Liq", "S_W"): 0,
            ("R16", "Liq", "H2O"): 0,
            ("R16", "Liq", "S_I"): 0,
            ("R16", "Liq", "X_I"): 0,
            ("R16", "Liq", "S_Z"): 0,
            # R17: Biodegradation of an inhibitory compound
            ("R17", "Liq", "X_S"): 0,
            ("R17", "Liq", "X_BH"): 1,
            ("R17", "Liq", "X_BA1"): 0,
            ("R17", "Liq", "X_BA2"): 0,
            ("R17", "Liq", "X_D"): 0,
            ("R17", "Liq", "S_S"): 0,
            ("R17", "Liq", "S_O"): -(1 - self.Y_H) / self.Y_H,
            ("R17", "Liq", "S_NO3"): 0,
            ("R17", "Liq", "S_NO2"): 0,
            ("R17", "Liq", "S_NO"): 0,
            ("R17", "Liq", "S_N2O"): 0,
            ("R17", "Liq", "S_N2"): 0,
            ("R17", "Liq", "S_NH"): -self.i_NXB,
            ("R17", "Liq", "S_NS"): 0,
            ("R17", "Liq", "X_NS"): 0,
            ("R17", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R17", "Liq", "S_J"): -1 / self.Y_H,
            ("R17", "Liq", "S_W"): 0,
            ("R17", "Liq", "H2O"): 0,
            ("R17", "Liq", "S_I"): 0,
            ("R17", "Liq", "X_I"): 0,
            ("R17", "Liq", "S_Z"): 0,
            # R18: Biodegradation of a special interest compound
            ("R18", "Liq", "X_S"): 0,
            ("R18", "Liq", "X_BH"): 1,
            ("R18", "Liq", "X_BA1"): 0,
            ("R18", "Liq", "X_BA2"): 0,
            ("R18", "Liq", "X_D"): 0,
            ("R18", "Liq", "S_S"): 0,
            ("R18", "Liq", "S_O"): -(1 - self.Y_H) / self.Y_H,
            ("R18", "Liq", "S_NO3"): 0,
            ("R18", "Liq", "S_NO2"): 0,
            ("R18", "Liq", "S_NO"): 0,
            ("R18", "Liq", "S_N2O"): 0,
            ("R18", "Liq", "S_N2"): 0,
            ("R18", "Liq", "S_NH"): -self.i_NXB,
            ("R18", "Liq", "S_NS"): 0,
            ("R18", "Liq", "X_NS"): 0,
            ("R18", "Liq", "S_ALK"): -self.i_NXB * (mw_alk / mw_n),
            ("R18", "Liq", "S_J"): 0,
            ("R18", "Liq", "S_W"): -1 / self.Y_H,
            ("R18", "Liq", "H2O"): 0,
            ("R18", "Liq", "S_I"): 0,
            ("R18", "Liq", "X_I"): 0,
            ("R18", "Liq", "S_Z"): 0,
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


class ASMNReactionScaler(CustomScalerBase):
    """
    Scaler for the Activated Sludge Model for Nitrogen reaction package.

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


class _ASMNReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ASMNReactionScaler

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


@declare_process_block_class("ASMNReactionBlock", block_class=_ASMNReactionBlock)
class ASMNReactionBlockData(ReactionBlockDataBase):
    """
    ReactionBlock for ASMN.
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

            # Free ammonia and free nitrous acid (Anthonisen et al., 1976)
            T_K = pyo.units.convert(
                self.state_ref.temperature / pyo.units.K,
                to_units=pyo.units.dimensionless,
            )
            pH = self.state_ref.pH
            K_B_over_K_W = pyo.exp(6344 / T_K)
            K_A = pyo.exp(-2300 / T_K)
            FA = self.conc_mass_comp_ref["S_NH"] * 10**pH / (K_B_over_K_W + 10**pH)
            FNA = self.conc_mass_comp_ref["S_NO2"] / (K_A * 10**pH)

            # Composite maximum specific growth rates with supplementary functions
            # Heterotrophs: ammonia limit is always present; salt inhibition is optional
            mu_H_eff = self.params.mu_H * (
                self.conc_mass_comp_ref["S_NH"]
                / (self.params.K_N1 + self.conc_mass_comp_ref["S_NH"])
            )
            if self.params.config.include_salt_inhibition:
                mu_H_eff = mu_H_eff * (
                    self.params.K_IHZ
                    / (self.params.K_IHZ + self.conc_mass_comp_ref["S_Z"])
                )

            # AOB: salt inhibition and AOB-specific organic inhibitor are both optional
            mu_A1_eff = self.params.mu_A1
            if self.params.config.include_aob_inhibitor:
                mu_A1_eff = mu_A1_eff * (
                    self.params.K_ISJ
                    / (self.params.K_ISJ + self.conc_mass_comp_ref["S_J"])
                )
            if self.params.config.include_salt_inhibition:
                mu_A1_eff = mu_A1_eff * (
                    self.params.K_IA1Z
                    / (self.params.K_IA1Z + self.conc_mass_comp_ref["S_Z"])
                )

            # NOB: salt inhibition is optional
            mu_A2_eff = self.params.mu_A2
            if self.params.config.include_salt_inhibition:
                mu_A2_eff = mu_A2_eff * (
                    self.params.K_IA2Z
                    / (self.params.K_IA2Z + self.conc_mass_comp_ref["S_Z"])
                )

            # Sum of nitrogen oxides for hydrolysis and R17/R18 anoxic switching
            sum_NOX = (
                self.conc_mass_comp_ref["S_NO3"]
                + self.conc_mass_comp_ref["S_NO2"]
                + self.conc_mass_comp_ref["S_NO"]
                + self.conc_mass_comp_ref["S_N2O"]
            )

            def rate_expression_rule(b, r):
                if r == "R1":
                    # R1: Aerobic growth of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S1 + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic growth of heterotrophs, NO3 -> NO2
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * b.params.eta_g2
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S2 + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.params.K_OH2
                            / (b.params.K_OH2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Anoxic growth of heterotrophs, NO2 -> NO
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * b.params.eta_g3
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S3 + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.K_NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * (
                            b.params.K_OH3
                            / (b.params.K_OH3 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.params.K_I3NO
                            / (b.params.K_I3NO + b.conc_mass_comp_ref["S_NO"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Anoxic growth of heterotrophs, NO -> N2O
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * b.params.eta_g4
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S4 + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (
                                b.params.K_NO
                                + b.conc_mass_comp_ref["S_NO"]
                                + b.conc_mass_comp_ref["S_NO"] ** 2 / b.params.K_I4NO
                            )
                        )
                        * (
                            b.params.K_OH4
                            / (b.params.K_OH4 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Anoxic growth of heterotrophs, N2O -> N2
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * b.params.eta_g5
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S5 + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_N2O"]
                            / (b.params.K_N2O + b.conc_mass_comp_ref["S_N2O"])
                        )
                        * (
                            b.params.K_OH5
                            / (b.params.K_OH5 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.params.K_I5NO
                            / (b.params.K_I5NO + b.conc_mass_comp_ref["S_NO"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: ANRA, nitrate reduction to nitrite
                    return b.reaction_rate[r] == pyo.units.convert(
                        1.2
                        * b.params.i_NXB
                        * (
                            b.conc_mass_comp_ref["S_NO3"]
                            / (b.params.K_6NO3 + b.conc_mass_comp_ref["S_NO3"])
                        )
                        * (
                            b.params.K_I6NH
                            / (b.params.K_I6NH + b.conc_mass_comp_ref["S_NH"])
                        )
                        * (
                            b.params.K_I6NO2
                            / (b.params.K_I6NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * (
                            b.reaction_rate["R1"]
                            + b.reaction_rate["R2"]
                            + b.reaction_rate["R3"]
                            + b.reaction_rate["R4"]
                            + b.reaction_rate["R5"]
                            + b.reaction_rate["R9"]
                            + b.reaction_rate["R10"]
                            + b.reaction_rate["R11"]
                            + b.reaction_rate["R17"]
                            + b.reaction_rate["R18"]
                            - b.reaction_rate["R14"]
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: ANRA, nitrite reduction to ammonia
                    return b.reaction_rate[r] == pyo.units.convert(
                        1.2
                        * b.params.i_NXB
                        * (
                            b.conc_mass_comp_ref["S_NO2"]
                            / (b.params.K_7NO2 + b.conc_mass_comp_ref["S_NO2"])
                        )
                        * (
                            b.params.K_I7NH
                            / (b.params.K_I7NH + b.conc_mass_comp_ref["S_NH"])
                        )
                        * (
                            b.reaction_rate["R1"]
                            + b.reaction_rate["R2"]
                            + b.reaction_rate["R3"]
                            + b.reaction_rate["R4"]
                            + b.reaction_rate["R5"]
                            + b.reaction_rate["R9"]
                            + b.reaction_rate["R10"]
                            + b.reaction_rate["R11"]
                            + b.reaction_rate["R17"]
                            + b.reaction_rate["R18"]
                            - b.reaction_rate["R14"]
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Death and lysis of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_LH * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Autotrophic growth of AOB
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_A1_eff
                        * (FA / (b.params.K_FA + FA + FA**2 / b.params.K_I9FA))
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_OA1 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (b.params.K_I9FNA / (b.params.K_I9FNA + FNA))
                        * b.conc_mass_comp_ref["X_BA1"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R10":
                    # R10: Autotrophic growth of NOB
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_A2_eff
                        * (FNA / (b.params.K_FNA + FNA + FNA**2 / b.params.K_I10FNA))
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_OA2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (b.params.K_I10FA / (b.params.K_I10FA + FA))
                        * b.conc_mass_comp_ref["X_BA2"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Mixotrophic growth of NOB
                    return b.reaction_rate[r] == (
                        b.params.beta_g
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S11 + b.conc_mass_comp_ref["S_S"])
                        )
                        * b.reaction_rate["R10"]
                    )
                elif r == "R12":
                    # R12: Death and lysis of AOB
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_LA1 * b.conc_mass_comp_ref["X_BA1"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R13":
                    # R13: Death and lysis of NOB
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_LA2 * b.conc_mass_comp_ref["X_BA2"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R14":
                    # R14: Ammonification of soluble organic nitrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_a
                        * b.conc_mass_comp_ref["S_NS"]
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R15":
                    # R15: Hydrolysis of particulate organics
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_h
                        * b.conc_mass_comp_ref["X_S"]
                        / (
                            b.params.K_X * b.conc_mass_comp_ref["X_BH"]
                            + b.conc_mass_comp_ref["X_S"]
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["S_O"]
                                / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            )
                            + b.params.eta_h
                            * b.params.K_OH1
                            / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            * (sum_NOX / (b.params.K_NO3 + sum_NOX))
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R16":
                    # R16: Hydrolysis of particulate organic nitrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.reaction_rate["R15"]
                        * (b.conc_mass_comp_ref["X_NS"] / b.conc_mass_comp_ref["X_S"]),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R17":
                    # R17: Biodegradation of an inhibitory compound
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * (b.conc_mass_comp_ref["S_J"] / b.conc_mass_comp_ref["S_S"])
                        * b.conc_mass_comp_ref["X_BH"]
                        * (
                            b.conc_mass_comp_ref["S_J"]
                            / (b.params.K_SJ + b.conc_mass_comp_ref["S_J"])
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["S_O"]
                                / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            )
                            + b.params.eta_g
                            * b.params.K_OH1
                            / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            * (sum_NOX / (b.params.K_NO3 + sum_NOX))
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R18":
                    # R18: Biodegradation of a special interest compound
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H_eff
                        * (b.conc_mass_comp_ref["S_W"] / b.conc_mass_comp_ref["S_S"])
                        * b.conc_mass_comp_ref["X_BH"]
                        * (
                            b.conc_mass_comp_ref["S_W"]
                            / (b.params.K_SW + b.conc_mass_comp_ref["S_W"])
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["S_O"]
                                / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            )
                            + b.params.eta_g
                            * b.params.K_OH1
                            / (b.params.K_OH1 + b.conc_mass_comp_ref["S_O"])
                            * (sum_NOX / (b.params.K_NO3 + sum_NOX))
                        ),
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                else:
                    raise BurntToast()

            self.rate_expression = pyo.Constraint(
                self.params.rate_reaction_idx,
                rule=rate_expression_rule,
                doc="ASMN rate expressions",
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
