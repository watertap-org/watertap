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
ASM1 reaction package.

References:

[1] Henze, M., Grady, C.P.L., Gujer, W., Marais, G.v.R., Matsuo, T.,
"Activated Sludge Model No. 1", 1987, IAWPRC Task Group on Mathematical Modeling
for Design and Operation of Biological Wastewater Treatment
[2] J. Alex, L. Benedetti, J. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons,
J.P. Steyer and P. Vanrolleghem, "Benchmark Simulation Model no. 1 (BSM1)", 2018
"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import BurntToast, ConfigurationError
import idaes.logger as idaeslog
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme

# Some more information about this module
__author__ = "Chenyu Wang, Adam Atia"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASM3ReactionParameterBlock")
class ASM3ReactionParameterData(ReactionParameterBlock):
    """
    Reaction Parameter Block Class
    """

    CONFIG = ReactionParameterBlock.CONFIG()

    CONFIG.declare(
        "reference_temperature",
        ConfigValue(
            default="20C",
            domain=In(["10C", "20C"]),
            description="Reference temperature for kinetic parameters",
            doc="""Indicates reference temperature for kinetic parameters
        **default** - "20C".
        **Valid values:** {
        **"10C"** - 10 Celsius degree,
        **"20C"** - 20 Celsius degree}""",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ASM3ReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ASM3 paper

        # R1: Hydrolysis

        # Heterotrophic organisms, aerobic and denitrifying activity
        # R2: Aerobic storage of S_S
        # R3: Anoxic storage of S_S
        # R4: Aerobic growth of X_H
        # R5: Anoxic growth (denitrific)
        # R6: Aerobic endog. respiration
        # R7: Anoxic endog. respiration
        # R8: Aerobic respiration of X_STO
        # R9: Anoxic respiration of X_STO

        # Autotrophic organisms, nitrifying activity
        # R10: Aerobic growth of X_A
        # R11: Aerobic endog. respiration
        # R12: Anoxic endog. respiration

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
                "R10",
                "R11",
                "R12",
            ]
        )

        # Stoichiometric Parameters
        self.Y_A = pyo.Var(
            initialize=0.24,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of autotrophic biomass per N03-N (g-COD-X_A / g-N-S_NOX)",
        )
        self.Y_H_O2 = pyo.Var(
            initialize=0.63,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Aerobic yield of heterotrophic biomass (g-COD-X_H / g-COD-X_STO)",
        )
        self.Y_H_NOX = pyo.Var(
            initialize=0.54,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic yield of heterotrophic biomass (g-COD-X_H / g-COD-X_STO)",
        )
        self.Y_STO_O2 = pyo.Var(
            initialize=0.85,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Aerobic yield of stored product per S_S (g-COD-X_STO / g-COD-S_S)",
        )
        self.Y_STO_NOX = pyo.Var(
            initialize=0.80,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic yield of stored product per per S_S (g-COD-X_STO / g-COD-S_S)",
        )

        add_object_reference(self, "f_SI", self.config.property_package.f_SI)
        add_object_reference(self, "f_XI", self.config.property_package.f_XI)
        add_object_reference(self, "i_NSI", self.config.property_package.i_NSI)
        add_object_reference(self, "i_NSS", self.config.property_package.i_NSS)
        add_object_reference(self, "i_NXI", self.config.property_package.i_NXI)
        add_object_reference(self, "i_NXS", self.config.property_package.i_NXS)
        add_object_reference(self, "i_NBM", self.config.property_package.i_NBM)
        add_object_reference(self, "i_SSXI", self.config.property_package.i_SSXI)
        add_object_reference(self, "i_SSXS", self.config.property_package.i_SSXS)
        add_object_reference(self, "i_SSBM", self.config.property_package.i_SSBM)
        add_object_reference(self, "i_SSSTO", self.config.property_package.i_SSSTO)

        # Kinetic Parameters
        k_H_dict = {"10C": 2, "20C": 3}
        self.k_H = pyo.Var(
            k_H_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=k_H_dict,
            units=pyo.units.day**-1,
            doc="Hydrolysis rate constant (g-COD-X_S / g-COD-X_H / day)",
        )
        self.K_X = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Hydrolysis saturation constant (g-COD-X_S / g-COD-X_H)",
        )

        # Heterotrophic organisms X_H, aerobic and denitrifying activity
        k_STO_dict = {"10C": 2.5, "20C": 5}
        self.k_STO = pyo.Var(
            k_STO_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=k_STO_dict,
            units=pyo.units.day**-1,
            doc="Storage rate constant (g-COD-S_S / g-COD-X_H / day)",
        )
        self.eta_NOX = pyo.Var(
            initialize=0.6,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Anoxic reduction factor",
        )
        self.K_O2 = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Saturation constant for S_NO2 (kg-O2 / m3)",
        )
        self.K_NOX = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Saturation constant for S_NOX (kg-N-NO3 / m3)",
        )
        self.K_S = pyo.Var(
            initialize=2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Saturation constant for substrate S_S (kg-COD-S_S / m3)",
        )
        self.K_STO = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Saturation constant for for X_STO (g-COD-X_STO / g-COD-X_H)",
        )
        mu_H_dict = {"10C": 1, "20C": 2}
        self.mu_H = pyo.Var(
            mu_H_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=mu_H_dict,
            units=pyo.units.day**-1,
            doc="Heterotrophic max. growth rate of X_H (day^-1)",
        )
        self.K_NH4 = pyo.Var(
            initialize=0.01e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Saturation constant for ammonium, S_NH4 (kg-N / m3)",
        )
        self.K_ALK = pyo.Var(
            initialize=0.1e-3,
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Saturation constant for alkalinity for X_H (kmol-HCO3- / m3)",
        )
        b_H_O2_dict = {"10C": 0.1, "20C": 0.2}
        self.b_H_O2 = pyo.Var(
            b_H_O2_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_H_O2_dict,
            units=pyo.units.day**-1,
            doc="Aerobic endogenous respiration rate of X_H (day^-1)",
        )
        b_H_NOX_dict = {"10C": 0.05, "20C": 0.1}
        self.b_H_NOX = pyo.Var(
            b_H_NOX_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_H_NOX_dict,
            units=pyo.units.day**-1,
            doc="Anoxic endogenous respiration rate of X_H (day^-1)",
        )
        b_STO_O2_dict = {"10C": 0.1, "20C": 0.2}
        self.b_STO_O2 = pyo.Var(
            b_STO_O2_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_STO_O2_dict,
            units=pyo.units.day**-1,
            doc="Aerobic respiration rate for X_STO (day^-1)",
        )
        b_STO_NOX_dict = {"10C": 0.05, "20C": 0.1}
        self.b_STO_NOX = pyo.Var(
            b_STO_NOX_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_STO_NOX_dict,
            units=pyo.units.day**-1,
            doc="Anoxic respiration rate for X_STO (day^-1)",
        )

        # Autotrophic organisms X_A, nitrifying activity
        mu_A_dict = {"10C": 0.35, "20C": 1}
        self.mu_A = pyo.Var(
            mu_A_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=mu_A_dict,
            units=pyo.units.day**-1,
            doc="Autotrophic max. growth rate of X_A (day^-1)",
        )
        self.K_A_NH4 = pyo.Var(
            initialize=1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Ammonium substrate saturation for X_A (kg-N / m3)",
        )
        self.K_A_O2 = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Oxygen saturation for nitrifiers (kg-O2 / m3)",
        )
        self.K_A_ALK = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kmol / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Bicarbonate saturation for nitrifiers (kmol-HCO3- / m3)",
        )
        b_A_O2_dict = {"10C": 0.05, "20C": 0.15}
        self.b_A_O2 = pyo.Var(
            b_A_O2_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_A_O2_dict,
            units=pyo.units.day**-1,
            doc="Aerobic endogenous respiration rate of X_A (day^-1)",
        )
        b_A_NOX_dict = {"10C": 0.02, "20C": 0.05}
        self.b_A_NOX = pyo.Var(
            b_A_NOX_dict.keys(),
            domain=pyo.PositiveReals,
            initialize=b_A_NOX_dict,
            units=pyo.units.day**-1,
            doc="Anoxic endogenous respiration rate of X_A (day^-1)",
        )

        # Stoichiometric numbers from Table 1
        # obtained by \sum_i^12 νji*ikI
        x1 = 1.0 - self.f_SI
        x2 = -1.0 + self.Y_STO_O2
        x3 = (-1.0 + self.Y_STO_NOX) / (64.0 / 14.0 - 24.0 / 14.0)
        x4 = 1.0 - 1.0 / self.Y_H_O2
        x5 = (+1.0 - 1.0 / self.Y_H_NOX) / (64.0 / 14.0 - 24.0 / 14.0)
        x6 = -1.0 + self.f_XI
        x7 = (self.f_XI - 1.0) / (64.0 / 14.0 - 24.0 / 14.0)
        x8 = -1.0
        x9 = -1.0 / (64.0 / 14.0 - 24.0 / 14.0)
        x10 = -(64.0 / 14.0) / self.Y_A + 1.0
        x11 = self.f_XI - 1.0
        x12 = (self.f_XI - 1.0) / (64.0 / 14.0 - 24.0 / 14.0)

        y1 = -self.f_SI * self.i_NSI - (1.0 - self.f_SI) * self.i_NSS + self.i_NXS
        y2 = self.i_NSS
        y3 = self.i_NSS
        y4 = -self.i_NBM
        y5 = -self.i_NBM
        y6 = -self.f_XI * self.i_NXI + self.i_NBM
        y7 = -self.f_XI * self.i_NXI + self.i_NBM
        y10 = -1.0 / self.Y_A - self.i_NBM
        y11 = -self.f_XI * self.i_NXI + self.i_NBM
        y12 = -self.f_XI * self.i_NXI + self.i_NBM

        z1 = y1 / 14.0
        z2 = y2 / 14.0
        z3 = y3 / 14.0 - x3 / 14.0
        z4 = y4 / 14.0
        z5 = y5 / 14.0 - x5 / 14.0
        z6 = y6 / 14.0
        z7 = y7 / 14.0 - x7 / 14.0
        z9 = -x9 / 14.0
        z10 = y10 / 14.0 - 1.0 / (self.Y_A * 14.0)
        z11 = y11 / 14.0
        z12 = y12 / 14.0 - x12 / 14.0

        t1 = -self.i_SSXS
        t2 = self.Y_STO_O2 * self.i_SSSTO
        t3 = self.Y_STO_NOX * self.i_SSSTO
        t4 = self.i_SSBM - 1.0 / self.Y_H_O2 * self.i_SSSTO
        t5 = self.i_SSBM - 1.0 / self.Y_H_NOX * self.i_SSSTO
        t6 = self.f_XI * self.i_SSXI - self.i_SSBM
        t7 = self.f_XI * self.i_SSXI - self.i_SSBM
        t8 = -self.i_SSSTO
        t9 = -self.i_SSSTO
        t10 = self.i_SSBM
        t11 = self.f_XI * self.i_SSXI - self.i_SSBM
        t12 = self.f_XI * self.i_SSXI - self.i_SSBM

        # Reaction Stoichiometry
        # This is the stoichiometric part the Peterson matrix in dict form
        # Note that reaction stoichiometry is on a mass basis.
        # For alkalinity, this requires converting the mass of nitrogen species
        # reacted to mass of alkalinity converted using a charge balance (effectively MW_C/MW_N)
        mw_alk = 12 * pyo.units.kg / pyo.units.kmol
        mw_n = 14 * pyo.units.kg / pyo.units.kmol
        self.rate_reaction_stoichiometry = {
            # R1: Hydrolysis
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_O"): 0,
            ("R1", "Liq", "S_I"): self.f_SI,
            ("R1", "Liq", "S_S"): x1,
            ("R1", "Liq", "S_NH4"): y1,
            ("R1", "Liq", "S_N2"): 0,
            ("R1", "Liq", "S_NOX"): 0,
            ("R1", "Liq", "S_ALK"): z1,
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_S"): -1,
            ("R1", "Liq", "X_H"): 0,
            ("R1", "Liq", "X_STO"): 0,
            ("R1", "Liq", "X_A"): 0,
            ("R1", "Liq", "X_TSS"): t1,
            # Heterotrophic organisms, aerobic and denitrifying activity
            # R2: Aerobic storage of S_S
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_O"): x2,
            ("R2", "Liq", "S_I"): 0,
            ("R2", "Liq", "S_S"): -1,
            ("R2", "Liq", "S_NH4"): y2,
            ("R2", "Liq", "S_N2"): 0,
            ("R2", "Liq", "S_NOX"): 0,
            ("R2", "Liq", "S_ALK"): z2,
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_S"): 0,
            ("R2", "Liq", "X_H"): 0,
            ("R2", "Liq", "X_STO"): self.Y_STO_O2,
            ("R2", "Liq", "X_A"): 0,
            ("R2", "Liq", "X_TSS"): t2,
            # R3: Anoxic storage of S_S
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_O"): 0,
            ("R3", "Liq", "S_I"): 0,
            ("R3", "Liq", "S_S"): -1,
            ("R3", "Liq", "S_NH4"): y3,
            ("R3", "Liq", "S_N2"): -x3,
            ("R3", "Liq", "S_NOX"): x3,
            ("R3", "Liq", "S_ALK"): z3,
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_S"): 0,
            ("R3", "Liq", "X_H"): 0,
            ("R3", "Liq", "X_STO"): self.Y_STO_NOX,
            ("R3", "Liq", "X_A"): 0,
            ("R3", "Liq", "X_TSS"): t3,
            # R4: Aerobic growth of X_H
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_O"): x4,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "S_S"): 0,
            ("R4", "Liq", "S_NH4"): y4,
            ("R4", "Liq", "S_N2"): 0,
            ("R4", "Liq", "S_NOX"): 0,
            ("R4", "Liq", "S_ALK"): z4,
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_S"): 0,
            ("R4", "Liq", "X_H"): 1,
            ("R4", "Liq", "X_STO"): -1 / self.Y_H_O2,
            ("R4", "Liq", "X_A"): 0,
            ("R4", "Liq", "X_TSS"): t4,
            # R5: Anoxic growth (denitrific.)
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_O"): 0,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "S_S"): 0,
            ("R5", "Liq", "S_NH4"): y4,
            ("R5", "Liq", "S_N2"): -x5,
            ("R5", "Liq", "S_NOX"): x5,
            ("R5", "Liq", "S_ALK"): z5,
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_S"): 0,
            ("R5", "Liq", "X_H"): 1,
            ("R5", "Liq", "X_STO"): -1 / self.Y_H_NOX,
            ("R5", "Liq", "X_A"): 0,
            ("R5", "Liq", "X_TSS"): t5,
            # R6: Aerobic endog. respiration
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_O"): x6,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "S_S"): 0,
            ("R6", "Liq", "S_NH4"): y6,
            ("R6", "Liq", "S_N2"): 0,
            ("R6", "Liq", "S_NOX"): 0,
            ("R6", "Liq", "S_ALK"): z6,
            ("R6", "Liq", "X_I"): self.f_XI,
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_H"): -1,
            ("R6", "Liq", "X_STO"): 0,
            ("R6", "Liq", "X_A"): 0,
            ("R6", "Liq", "X_TSS"): t6,
            # R7: Anoxic endog. respiration
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_O"): 0,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "S_S"): 0,
            ("R7", "Liq", "S_NH4"): y7,
            ("R7", "Liq", "S_N2"): -x7,
            ("R7", "Liq", "S_NOX"): x7,
            ("R7", "Liq", "S_ALK"): z7,
            ("R7", "Liq", "X_I"): self.f_XI,
            ("R7", "Liq", "X_S"): 0,
            ("R7", "Liq", "X_H"): -1,
            ("R7", "Liq", "X_STO"): 0,
            ("R7", "Liq", "X_A"): 0,
            ("R7", "Liq", "X_TSS"): t7,
            # R8: Aerobic respiration of X_STO
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_O"): x8,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "S_S"): 0,
            ("R8", "Liq", "S_NH4"): 0,
            ("R8", "Liq", "S_N2"): 0,
            ("R8", "Liq", "S_NOX"): 0,
            ("R8", "Liq", "S_ALK"): 0,
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_S"): 0,
            ("R8", "Liq", "X_H"): 0,
            ("R8", "Liq", "X_STO"): -1,
            ("R8", "Liq", "X_A"): 0,
            ("R8", "Liq", "X_TSS"): t8,
            # R9: Anoxic respiration of X_STO
            ("R9", "Liq", "H2O"): 0,
            ("R9", "Liq", "S_O"): 0,
            ("R9", "Liq", "S_I"): 0,
            ("R9", "Liq", "S_S"): 0,
            ("R9", "Liq", "S_NH4"): 0,
            ("R9", "Liq", "S_N2"): -x9,
            ("R9", "Liq", "S_NOX"): x9,
            ("R9", "Liq", "S_ALK"): z9,
            ("R9", "Liq", "X_I"): 0,
            ("R9", "Liq", "X_S"): 0,
            ("R9", "Liq", "X_H"): 0,
            ("R9", "Liq", "X_STO"): -1,
            ("R9", "Liq", "X_A"): 0,
            ("R9", "Liq", "X_TSS"): t9,
            # Autotrophic organisms, nitrifying activity
            # R10: Aerobic growth of X_A
            ("R10", "Liq", "H2O"): 0,
            ("R10", "Liq", "S_O"): x10,
            ("R10", "Liq", "S_I"): 0,
            ("R10", "Liq", "S_S"): 0,
            ("R10", "Liq", "S_NH4"): y10,
            ("R10", "Liq", "S_N2"): 0,
            ("R10", "Liq", "S_NOX"): 1 / self.Y_A,
            ("R10", "Liq", "S_ALK"): z10,
            ("R10", "Liq", "X_I"): 0,
            ("R10", "Liq", "X_S"): 0,
            ("R10", "Liq", "X_H"): 0,
            ("R10", "Liq", "X_STO"): 0,
            ("R10", "Liq", "X_A"): 1,
            ("R10", "Liq", "X_TSS"): t10,
            # R11: Aerobic endog. respiration
            ("R11", "Liq", "H2O"): 0,
            ("R11", "Liq", "S_O"): x11,
            ("R11", "Liq", "S_I"): 0,
            ("R11", "Liq", "S_S"): 0,
            ("R11", "Liq", "S_NH4"): y11,
            ("R11", "Liq", "S_N2"): 0,
            ("R11", "Liq", "S_NOX"): 0,
            ("R11", "Liq", "S_ALK"): z11,
            ("R11", "Liq", "X_I"): self.f_XI,
            ("R11", "Liq", "X_S"): 0,
            ("R11", "Liq", "X_H"): 0,
            ("R11", "Liq", "X_STO"): 0,
            ("R11", "Liq", "X_A"): -1,
            ("R11", "Liq", "X_TSS"): t11,
            # R12: Anoxic endog. respiration
            ("R12", "Liq", "H2O"): 0,
            ("R12", "Liq", "S_O"): 0,
            ("R12", "Liq", "S_I"): 0,
            ("R12", "Liq", "S_S"): 0,
            ("R12", "Liq", "S_NH4"): y12,
            ("R12", "Liq", "S_N2"): -x12,
            ("R12", "Liq", "S_NOX"): x12,
            ("R12", "Liq", "S_ALK"): z12,
            ("R12", "Liq", "X_I"): self.f_XI,
            ("R12", "Liq", "X_S"): 0,
            ("R12", "Liq", "X_H"): 0,
            ("R12", "Liq", "X_STO"): 0,
            ("R12", "Liq", "X_A"): -1,
            ("R12", "Liq", "X_TSS"): t12,
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


class ASM3ReactionScaler(CustomScalerBase):
    """
    Scaler for the Activated Sludge Model No.1 reaction package.

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


class _ASM3ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    default_scaler = ASM3ReactionScaler

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


@declare_process_block_class("ASM3ReactionBlock", block_class=_ASM3ReactionBlock)
class ASM3ReactionBlockData(ReactionBlockDataBase):
    """
    ReactionBlock for ASM3.
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
                if self.params.config.reference_temperature == "20C":
                    k_H = b.params.k_H["20C"]
                    k_STO = b.params.k_STO["20C"]
                    mu_H = b.params.mu_H["20C"]
                    b_H_O2 = b.params.b_H_O2["20C"]
                    b_H_NOX = b.params.b_H_NOX["20C"]
                    b_STO_O2 = b.params.b_STO_O2["20C"]
                    b_STO_NOX = b.params.b_STO_NOX["20C"]
                    mu_A = b.params.mu_A["20C"]
                    b_A_O2 = b.params.b_A_O2["20C"]
                    b_A_NOX = b.params.b_A_NOX["20C"]
                elif self.params.config.reference_temperature == "10C":
                    k_H = b.params.k_H["10C"]
                    k_STO = b.params.k_STO["10C"]
                    mu_H = b.params.mu_H["10C"]
                    b_H_O2 = b.params.b_H_O2["10C"]
                    b_H_NOX = b.params.b_H_NOX["10C"]
                    b_STO_O2 = b.params.b_STO_O2["10C"]
                    b_STO_NOX = b.params.b_STO_NOX["10C"]
                    mu_A = b.params.mu_A["10C"]
                    b_A_O2 = b.params.b_A_O2["10C"]
                    b_A_NOX = b.params.b_A_NOX["10C"]
                else:
                    raise ConfigurationError(
                        "Reference temperature only supports '10C' and '20C'"
                    )

                if r == "R1":
                    # R1: Hydrolysis
                    return b.reaction_rate[r] == pyo.units.convert(
                        k_H
                        * (b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"])
                        / (
                            b.params.K_X
                            + b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_H"]
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                # Heterotrophic organisms, aerobic and denitrifying activity
                elif r == "R2":
                    # R2: Aerobic storage of S_S
                    return b.reaction_rate[r] == pyo.units.convert(
                        k_STO
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S + b.conc_mass_comp_ref["S_S"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Anoxic storage of S_S
                    return b.reaction_rate[r] == pyo.units.convert(
                        k_STO
                        * b.params.eta_NOX
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NOX"]
                            / (b.params.K_NOX + b.conc_mass_comp_ref["S_NOX"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S + b.conc_mass_comp_ref["S_S"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Aerobic growth
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * (b.conc_mass_comp_ref["X_STO"] / b.conc_mass_comp_ref["X_H"])
                        / (
                            b.params.K_STO
                            + b.conc_mass_comp_ref["X_STO"]
                            / b.conc_mass_comp_ref["X_H"]
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Anoxic growth
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_H
                        * b.params.eta_NOX
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NOX"]
                            / (b.params.K_NOX + b.conc_mass_comp_ref["S_NOX"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_ALK + b.state_ref.alkalinity)
                        )
                        * (b.conc_mass_comp_ref["X_STO"] / b.conc_mass_comp_ref["X_H"])
                        / (
                            b.params.K_STO
                            + b.conc_mass_comp_ref["X_STO"]
                            / b.conc_mass_comp_ref["X_H"]
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Aerobic endogenous respiration
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_H_O2
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Anoxic endogenous respiration
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_H_NOX
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NOX"]
                            / (b.params.K_NOX + b.conc_mass_comp_ref["S_NOX"])
                        )
                        * b.conc_mass_comp_ref["X_H"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Aerobic respiration of X_STO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_STO_O2
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_STO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R9":
                    # R9: Anoxic respiration of X_STO
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_STO_NOX
                        * (
                            b.params.K_O2
                            / (b.params.K_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NOX"]
                            / (b.params.K_NOX + b.conc_mass_comp_ref["S_NOX"])
                        )
                        * b.conc_mass_comp_ref["X_STO"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                # Autotrophic organisms, nitrifying activity
                elif r == "R10":
                    # R10: Aerobic growth of X_A, nitrification
                    return b.reaction_rate[r] == pyo.units.convert(
                        mu_A
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_A_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NH4"]
                            / (b.params.K_A_NH4 + b.conc_mass_comp_ref["S_NH4"])
                        )
                        * (
                            b.state_ref.alkalinity
                            / (b.params.K_A_ALK + b.state_ref.alkalinity)
                        )
                        * b.conc_mass_comp_ref["X_A"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R11":
                    # R11: Aerobic endogenous respiration
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_A_O2
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_A_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_A"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R12":
                    # R12: Anoxic endogenous respiration
                    return b.reaction_rate[r] == pyo.units.convert(
                        b_A_NOX
                        * (
                            b.params.K_A_O2
                            / (b.params.K_A_O2 + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NOX"]
                            / (b.params.K_NOX + b.conc_mass_comp_ref["S_NOX"])
                        )
                        * b.conc_mass_comp_ref["X_A"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                else:
                    raise BurntToast()

            self.rate_expression = pyo.Constraint(
                self.params.rate_reaction_idx,
                rule=rate_expression_rule,
                doc="ASM3 rate expressions",
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
        # iscale.constraint_scaling_transform(self.rate_expression["R5"], 1e3)
        # iscale.constraint_scaling_transform(self.rate_expression["R3"], 1e3)
        # iscale.constraint_scaling_transform(self.rate_expression["R4"], 1e3)
