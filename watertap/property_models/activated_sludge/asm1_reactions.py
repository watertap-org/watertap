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


# Some more information about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASM1ReactionParameterBlock")
class ASM1ReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = ASM1ReactionBlock

        # Reaction Index
        # Reaction names based on standard numbering in ASM1 paper
        # R1: Aerobic growth of heterotrophs
        # R2: Anoxic growth of heterotrophs
        # R3: Aerobic growth of autotrophs
        # R4: Decay of heterotrophs
        # R5: Decay of autotrophs
        # R6: Ammonification of soluble organic nitrogen
        # R7: Hydrolysis of entrapped organics
        # R8: Hydrolysis of entrapped organic nitrogen
        self.rate_reaction_idx = pyo.Set(
            initialize=["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
        )

        # Stoichiometric Parameters
        self.Y_A = pyo.Var(
            initialize=0.24,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of cell COD formed per g N consumed, Y_A",
        )
        self.Y_H = pyo.Var(
            initialize=0.67,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of cell COD formed per g COD oxidized, Y_H",
        )
        self.f_p = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of biomass yielding particulate products, f_p",
        )
        self.i_xb = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Mass fraction of N per COD in biomass, i_xb",
        )
        self.i_xp = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Mass fraction of N per COD in particulates, i_xp",
        )

        # Kinetic Parameters
        self.mu_A = pyo.Var(
            initialize=0.5,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for autotrophic biomass, mu_A",
        )
        self.mu_H = pyo.Var(
            initialize=4,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for heterotrophic biomass, mu_H",
        )
        self.K_S = pyo.Var(
            initialize=10e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for heterotrophic biomass, K_S",
        )
        self.K_OH = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Oxygen half-saturation coefficient for heterotrophic biomass, K_O,H",
        )
        self.K_OA = pyo.Var(
            initialize=0.4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Oxygen half-saturation coefficient for autotrophic biomass, K_O,A",
        )
        self.K_NO = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitrate half-saturation coefficient for denitrifying heterotrophic biomass, K_NO",
        )
        self.b_H = pyo.Var(
            initialize=0.3,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient for heterotrophic biomass, b_H",
        )
        self.b_A = pyo.Var(
            initialize=0.05,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient for autotrophic biomass, b_A",
        )
        self.eta_g = pyo.Var(
            initialize=0.8,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Correction factor for mu_H under anoxic conditions, eta_g",
        )
        self.eta_h = pyo.Var(
            initialize=0.8,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Correction factor for hydrolysis under anoxic conditions, eta_h",
        )
        self.k_h = pyo.Var(
            initialize=3,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific hydrolysis rate, k_h",
        )
        self.K_X = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for hydrolysis of slowly biodegradable substrate, K_X",
        )
        self.K_NH = pyo.Var(
            initialize=1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Ammonia Half-saturation coefficient for autotrophic biomass, K_NH",
        )
        self.k_a = pyo.Var(
            initialize=0.05e3,
            units=pyo.units.m**3 / pyo.units.kg / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Ammonification rate, k_a",
        )

        # Reaction Stoichiometry
        # This is the stoichiometric part the Peterson matrix in dict form
        # Note that reaction stoichiometry is on a mass basis.
        # For alkalinity, this requires converting the mass of nitrogen species
        # reacted to mass of alkalinity converted using a charge balance (effectively MW_C/MW_N)
        mw_alk = 12 * pyo.units.kg / pyo.units.kmol
        mw_n = 14 * pyo.units.kg / pyo.units.kmol
        self.rate_reaction_stoichiometry = {
            # R1: Aerobic growth of heterotrophs
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "S_I"): 0,
            ("R1", "Liq", "S_S"): -1 / self.Y_H,
            ("R1", "Liq", "X_I"): 0,
            ("R1", "Liq", "X_S"): 0,
            ("R1", "Liq", "X_BH"): 1,
            ("R1", "Liq", "X_BA"): 0,
            ("R1", "Liq", "X_P"): 0,
            ("R1", "Liq", "S_O"): -(1 - self.Y_H) / self.Y_H,
            ("R1", "Liq", "S_NO"): 0,
            ("R1", "Liq", "S_NH"): -self.i_xb,
            ("R1", "Liq", "S_ND"): 0,
            ("R1", "Liq", "X_ND"): 0,
            ("R1", "Liq", "S_ALK"): -self.i_xb * (mw_alk / mw_n),
            # R2: Anoxic growth of heterotrophs
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "S_I"): 0,
            ("R2", "Liq", "S_S"): -1 / self.Y_H,
            ("R2", "Liq", "X_I"): 0,
            ("R2", "Liq", "X_S"): 0,
            ("R2", "Liq", "X_BH"): 1,
            ("R2", "Liq", "X_BA"): 0,
            ("R2", "Liq", "X_P"): 0,
            ("R2", "Liq", "S_O"): 0,
            ("R2", "Liq", "S_NO"): -(1 - self.Y_H) / (2.86 * self.Y_H),
            ("R2", "Liq", "S_NH"): -self.i_xb,
            ("R2", "Liq", "S_ND"): 0,
            ("R2", "Liq", "X_ND"): 0,
            ("R2", "Liq", "S_ALK"): ((1 - self.Y_H) / (2.86 * self.Y_H) - self.i_xb)
            * (mw_alk / mw_n),
            # R3: Aerobic growth of autotrophs
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "S_I"): 0,
            ("R3", "Liq", "S_S"): 0,
            ("R3", "Liq", "X_I"): 0,
            ("R3", "Liq", "X_S"): 0,
            ("R3", "Liq", "X_BH"): 0,
            ("R3", "Liq", "X_BA"): 1,
            ("R3", "Liq", "X_P"): 0,
            ("R3", "Liq", "S_O"): -(4.57 - self.Y_A) / self.Y_A,
            ("R3", "Liq", "S_NO"): 1 / self.Y_A,
            ("R3", "Liq", "S_NH"): -self.i_xb - 1 / self.Y_A,
            ("R3", "Liq", "S_ND"): 0,
            ("R3", "Liq", "X_ND"): 0,
            ("R3", "Liq", "S_ALK"): (-self.i_xb - 2 / (self.Y_A)) * (mw_alk / mw_n),
            # R4: Decay of heterotrophs
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "S_I"): 0,
            ("R4", "Liq", "S_S"): 0,
            ("R4", "Liq", "X_I"): 0,
            ("R4", "Liq", "X_S"): 1 - self.f_p,
            ("R4", "Liq", "X_BH"): -1,
            ("R4", "Liq", "X_BA"): 0,
            ("R4", "Liq", "X_P"): self.f_p,
            ("R4", "Liq", "S_O"): 0,
            ("R4", "Liq", "S_NO"): 0,
            ("R4", "Liq", "S_NH"): 0,
            ("R4", "Liq", "S_ND"): 0,
            ("R4", "Liq", "X_ND"): self.i_xb - self.f_p * self.i_xp,
            ("R4", "Liq", "S_ALK"): 0,
            # R5: Decay of autotrophs
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "S_I"): 0,
            ("R5", "Liq", "S_S"): 0,
            ("R5", "Liq", "X_I"): 0,
            ("R5", "Liq", "X_S"): 1 - self.f_p,
            ("R5", "Liq", "X_BH"): 0,
            ("R5", "Liq", "X_BA"): -1,
            ("R5", "Liq", "X_P"): self.f_p,
            ("R5", "Liq", "S_O"): 0,
            ("R5", "Liq", "S_NO"): 0,
            ("R5", "Liq", "S_NH"): 0,
            ("R5", "Liq", "S_ND"): 0,
            ("R5", "Liq", "X_ND"): self.i_xb - self.f_p * self.i_xp,
            ("R5", "Liq", "S_ALK"): 0,
            # R6: Ammonification of soluble organic nitrogen
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "S_I"): 0,
            ("R6", "Liq", "S_S"): 0,
            ("R6", "Liq", "X_I"): 0,
            ("R6", "Liq", "X_S"): 0,
            ("R6", "Liq", "X_BH"): 0,
            ("R6", "Liq", "X_BA"): 0,
            ("R6", "Liq", "X_P"): 0,
            ("R6", "Liq", "S_O"): 0,
            ("R6", "Liq", "S_NO"): 0,
            ("R6", "Liq", "S_NH"): 1,
            ("R6", "Liq", "S_ND"): -1,
            ("R6", "Liq", "X_ND"): 0,
            ("R6", "Liq", "S_ALK"): 1 * (mw_alk / mw_n),
            # R7: Hydrolysis of entrapped organics
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "S_I"): 0,
            ("R7", "Liq", "S_S"): 1,
            ("R7", "Liq", "X_I"): 0,
            ("R7", "Liq", "X_S"): -1,
            ("R7", "Liq", "X_BH"): 0,
            ("R7", "Liq", "X_BA"): 0,
            ("R7", "Liq", "X_P"): 0,
            ("R7", "Liq", "S_O"): 0,
            ("R7", "Liq", "S_NO"): 0,
            ("R7", "Liq", "S_NH"): 0,
            ("R7", "Liq", "S_ND"): 0,
            ("R7", "Liq", "X_ND"): 0,
            ("R7", "Liq", "S_ALK"): 0,
            # R8: Hydrolysis of entrapped organic nitrogen
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "S_I"): 0,
            ("R8", "Liq", "S_S"): 0,
            ("R8", "Liq", "X_I"): 0,
            ("R8", "Liq", "X_S"): 0,
            ("R8", "Liq", "X_BH"): 0,
            ("R8", "Liq", "X_BA"): 0,
            ("R8", "Liq", "X_P"): 0,
            ("R8", "Liq", "S_O"): 0,
            ("R8", "Liq", "S_NO"): 0,
            ("R8", "Liq", "S_NH"): 0,
            ("R8", "Liq", "S_ND"): 1,
            ("R8", "Liq", "X_ND"): -1,
            ("R8", "Liq", "S_ALK"): 0,
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


class _ASM1ReactionBlock(ReactionBlockBase):
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


@declare_process_block_class("ASM1ReactionBlock", block_class=_ASM1ReactionBlock)
class ASM1ReactionBlockData(ReactionBlockDataBase):
    """
    ReactionBlcok for ASM1.
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
                    # R1: Aerobic growth of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_OH + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic growth of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_H
                        * (
                            b.conc_mass_comp_ref["S_S"]
                            / (b.params.K_S + b.conc_mass_comp_ref["S_S"])
                        )
                        * (
                            b.params.K_OH
                            / (b.params.K_OH + b.conc_mass_comp_ref["S_O"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_NO"]
                            / (b.params.K_NO + b.conc_mass_comp_ref["S_NO"])
                        )
                        * b.params.eta_g
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Aerobic growth of autotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.mu_A
                        * (
                            b.conc_mass_comp_ref["S_NH"]
                            / (b.params.K_NH + b.conc_mass_comp_ref["S_NH"])
                        )
                        * (
                            b.conc_mass_comp_ref["S_O"]
                            / (b.params.K_OA + b.conc_mass_comp_ref["S_O"])
                        )
                        * b.conc_mass_comp_ref["X_BA"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Decay of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_H * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Decay of autotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.b_A * b.conc_mass_comp_ref["X_BA"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Ammonification of soluble organic nitrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_a
                        * b.conc_mass_comp_ref["S_ND"]
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Hydrolysis of entrapped organics
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.k_h
                        * (b.conc_mass_comp_ref["X_S"] / b.conc_mass_comp_ref["X_BH"])
                        / (
                            b.params.K_X
                            + (
                                b.conc_mass_comp_ref["X_S"]
                                / b.conc_mass_comp_ref["X_BH"]
                            )
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["S_O"]
                                / (b.params.K_OH + b.conc_mass_comp_ref["S_O"])
                            )
                            + b.params.eta_h
                            * b.params.K_OH
                            / (b.params.K_OH + b.conc_mass_comp_ref["S_O"])
                            * (
                                b.conc_mass_comp_ref["S_NO"]
                                / (b.params.K_NO + b.conc_mass_comp_ref["S_NO"])
                            )
                        )
                        * b.conc_mass_comp_ref["X_BH"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Hydrolysis of entrapped organic nitrogen
                    return b.reaction_rate[r] == (
                        b.reaction_rate["R7"]
                        * (b.conc_mass_comp_ref["X_ND"] / b.conc_mass_comp_ref["X_S"])
                    )
                else:
                    raise BurntToast()

            self.rate_expression = pyo.Constraint(
                self.params.rate_reaction_idx,
                rule=rate_expression_rule,
                doc="ASM1 rate expressions",
            )

        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.rate_expression)
            raise

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.mass
