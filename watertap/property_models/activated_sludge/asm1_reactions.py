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


# Some more inforation about this module
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
        self.yield_A = pyo.Var(
            initialize=0.24,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of cell COD formed per g N consumed, Y_A",
        )
        self.yield_H = pyo.Var(
            initialize=0.67,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Yield of cell COD formed per g COD oxidized, Y_H",
        )
        self.f_p = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of biomass yielding particulate products, y_p",
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
        self.max_growth_rate_autotrophic = pyo.Var(
            initialize=0.5,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for autotrophic biomass, mu_A",
        )
        self.max_growth_rate_heterotrophic = pyo.Var(
            initialize=4,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific growth rate for heterotrophic biomass, mu_H",
        )
        self.half_saturation_coeff_heterotrophic = pyo.Var(
            initialize=10e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for heterotrophic biomass, K_S",
        )
        self.half_saturation_coeff_oxygen_heterotrophic = pyo.Var(
            initialize=0.2e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Oxygen half-saturation coefficient for heterotrophic biomass, K_O,H",
        )
        self.half_saturation_coeff_oxygen_autotrophic = pyo.Var(
            initialize=0.4e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Oxygen half-saturation coefficient for autotrophic biomass, K_O,A",
        )
        self.half_saturation_coeff_nitrate_heterotrophic = pyo.Var(
            initialize=0.5e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Nitrate half-saturation coefficient for denitrifying heterotrophic biomass, K_NO",
        )
        self.decay_coeff_heterotrophic = pyo.Var(
            initialize=0.3,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Decay coefficient for heterotrophic biomass, b_H",
        )
        self.decay_coeff_autotrophic = pyo.Var(
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
            doc="Correction factor for hydolysis under anoxic conditions, eta_h",
        )
        self.max_hydrolysis_rate = pyo.Var(
            initialize=3,
            units=1 / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Maximum specific hydrolysis rate, k_h",
        )
        self.half_saturation_coeff_hydrolysis = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Half-saturation coefficient for hydolysis of slowly biodegradable substrate, K_X",
        )
        self.half_saturation_coeff_ammonium = pyo.Var(
            initialize=1e-3,
            units=pyo.units.kg / pyo.units.m**3,
            domain=pyo.PositiveReals,
            doc="Ammonia Half-saturation coefficient for autotrophic biomass, K_NH",
        )
        self.ammonification_rate = pyo.Var(
            initialize=0.05e3,
            units=pyo.units.m**3 / pyo.units.kg / pyo.units.day,
            domain=pyo.PositiveReals,
            doc="Ammonification rate, k_a",
        )

        # Reaction Stoichiometry
        # This is the stoichiometric part the Perterson matrix in dict form
        self.rate_reaction_stoichiometry = {
            # R1: Aerobic growth of heterotrophs
            ("R1", "Liq", "H2O"): 0,
            ("R1", "Liq", "inert_soluble"): 0,
            ("R1", "Liq", "substrate"): -1 / self.yield_H,
            ("R1", "Liq", "inert_particulate"): 0,
            ("R1", "Liq", "biodegradable"): 0,
            ("R1", "Liq", "heterotrophic"): 1,
            ("R1", "Liq", "autotrophic"): 0,
            ("R1", "Liq", "decay_particulate"): 0,
            ("R1", "Liq", "oxygen"): -(1 - self.yield_H) / self.yield_H,
            ("R1", "Liq", "nitrates"): 0,
            ("R1", "Liq", "ammonium"): -self.i_xb,
            ("R1", "Liq", "nitrogen_soluble"): 0,
            ("R1", "Liq", "nitrogen_particulate"): 0,
            ("R1", "Liq", "alkalinity"): -self.i_xb / 14,
            # R2: Anoxic growth of heterotrophs
            ("R2", "Liq", "H2O"): 0,
            ("R2", "Liq", "inert_soluble"): 0,
            ("R2", "Liq", "substrate"): -1 / self.yield_H,
            ("R2", "Liq", "inert_particulate"): 0,
            ("R2", "Liq", "biodegradable"): 0,
            ("R2", "Liq", "heterotrophic"): 1,
            ("R2", "Liq", "autotrophic"): 0,
            ("R2", "Liq", "decay_particulate"): 0,
            ("R2", "Liq", "oxygen"): 0,
            ("R2", "Liq", "nitrates"): -(1 - self.yield_H) / (2.86 * self.yield_H),
            ("R2", "Liq", "ammonium"): -self.i_xb,
            ("R2", "Liq", "nitrogen_soluble"): 0,
            ("R2", "Liq", "nitrogen_particulate"): 0,
            ("R2", "Liq", "alkalinity"): (1 - self.yield_H) / (14 * 2.86 * self.yield_H)
            - self.i_xb / 14,
            # R3: Aerobic growth of autotrophs
            ("R3", "Liq", "H2O"): 0,
            ("R3", "Liq", "inert_soluble"): 0,
            ("R3", "Liq", "substrate"): 0,
            ("R3", "Liq", "inert_particulate"): 0,
            ("R3", "Liq", "biodegradable"): 0,
            ("R3", "Liq", "heterotrophic"): 0,
            ("R3", "Liq", "autotrophic"): 1,
            ("R3", "Liq", "decay_particulate"): 0,
            ("R3", "Liq", "oxygen"): -(4.57 - self.yield_A) / self.yield_A,
            ("R3", "Liq", "nitrates"): 1 / self.yield_A,
            ("R3", "Liq", "ammonium"): -self.i_xb - 1 / self.yield_A,
            ("R3", "Liq", "nitrogen_soluble"): 0,
            ("R3", "Liq", "nitrogen_particulate"): 0,
            ("R3", "Liq", "alkalinity"): -self.i_xb / 14 - 1 / (7 * self.yield_A),
            # R4: Decay of heterotrophs
            ("R4", "Liq", "H2O"): 0,
            ("R4", "Liq", "inert_soluble"): 0,
            ("R4", "Liq", "substrate"): 0,
            ("R4", "Liq", "inert_particulate"): 0,
            ("R4", "Liq", "biodegradable"): 1 - self.f_p,
            ("R4", "Liq", "heterotrophic"): -1,
            ("R4", "Liq", "autotrophic"): 0,
            ("R4", "Liq", "decay_particulate"): self.f_p,
            ("R4", "Liq", "oxygen"): 0,
            ("R4", "Liq", "nitrates"): 0,
            ("R4", "Liq", "ammonium"): 0,
            ("R4", "Liq", "nitrogen_soluble"): 0,
            ("R4", "Liq", "nitrogen_particulate"): self.i_xb - self.f_p * self.i_xp,
            ("R4", "Liq", "alkalinity"): 0,
            # R5: Decay of autotrophs
            ("R5", "Liq", "H2O"): 0,
            ("R5", "Liq", "inert_soluble"): 0,
            ("R5", "Liq", "substrate"): 0,
            ("R5", "Liq", "inert_particulate"): 0,
            ("R5", "Liq", "biodegradable"): 1 - self.f_p,
            ("R5", "Liq", "heterotrophic"): 0,
            ("R5", "Liq", "autotrophic"): -1,
            ("R5", "Liq", "decay_particulate"): self.f_p,
            ("R5", "Liq", "oxygen"): 0,
            ("R5", "Liq", "nitrates"): 0,
            ("R5", "Liq", "ammonium"): 0,
            ("R5", "Liq", "nitrogen_soluble"): 0,
            ("R5", "Liq", "nitrogen_particulate"): self.i_xb - self.f_p * self.i_xp,
            ("R5", "Liq", "alkalinity"): 0,
            # R6: Ammonification of soluble organic nitrogen
            ("R6", "Liq", "H2O"): 0,
            ("R6", "Liq", "inert_soluble"): 0,
            ("R6", "Liq", "substrate"): 0,
            ("R6", "Liq", "inert_particulate"): 0,
            ("R6", "Liq", "biodegradable"): 0,
            ("R6", "Liq", "heterotrophic"): 0,
            ("R6", "Liq", "autotrophic"): 0,
            ("R6", "Liq", "decay_particulate"): 0,
            ("R6", "Liq", "oxygen"): 0,
            ("R6", "Liq", "nitrates"): 0,
            ("R6", "Liq", "ammonium"): 1,
            ("R6", "Liq", "nitrogen_soluble"): -1,
            ("R6", "Liq", "nitrogen_particulate"): 0,
            ("R6", "Liq", "alkalinity"): 1 / 14,
            # R7: Hydrolysis of entrapped organics
            ("R7", "Liq", "H2O"): 0,
            ("R7", "Liq", "inert_soluble"): 0,
            ("R7", "Liq", "substrate"): 1,
            ("R7", "Liq", "inert_particulate"): 0,
            ("R7", "Liq", "biodegradable"): -1,
            ("R7", "Liq", "heterotrophic"): 0,
            ("R7", "Liq", "autotrophic"): 0,
            ("R7", "Liq", "decay_particulate"): 0,
            ("R7", "Liq", "oxygen"): 0,
            ("R7", "Liq", "nitrates"): 0,
            ("R7", "Liq", "ammonium"): 0,
            ("R7", "Liq", "nitrogen_soluble"): 0,
            ("R7", "Liq", "nitrogen_particulate"): 0,
            ("R7", "Liq", "alkalinity"): 0,
            # R8: Hydrolysis of entrapped organic nitrogen
            ("R8", "Liq", "H2O"): 0,
            ("R8", "Liq", "inert_soluble"): 0,
            ("R8", "Liq", "substrate"): 0,
            ("R8", "Liq", "inert_particulate"): 0,
            ("R8", "Liq", "biodegradable"): 0,
            ("R8", "Liq", "heterotrophic"): 0,
            ("R8", "Liq", "autotrophic"): 0,
            ("R8", "Liq", "decay_particulate"): 0,
            ("R8", "Liq", "oxygen"): 0,
            ("R8", "Liq", "nitrates"): 0,
            ("R8", "Liq", "ammonium"): 0,
            ("R8", "Liq", "nitrogen_soluble"): 1,
            ("R8", "Liq", "nitrogen_particulate"): -1,
            ("R8", "Liq", "alkalinity"): 0,
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
                "amount": pyo.units.mol,
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
                        b.params.max_growth_rate_heterotrophic
                        * (
                            b.conc_mass_comp_ref["substrate"]
                            / (
                                b.params.half_saturation_coeff_heterotrophic
                                + b.conc_mass_comp_ref["substrate"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["oxygen"]
                            / (
                                b.params.half_saturation_coeff_oxygen_heterotrophic
                                + b.conc_mass_comp_ref["oxygen"]
                            )
                        )
                        * b.conc_mass_comp_ref["heterotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R2":
                    # R2: Anoxic growth of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.max_growth_rate_heterotrophic
                        * (
                            b.conc_mass_comp_ref["substrate"]
                            / (
                                b.params.half_saturation_coeff_heterotrophic
                                + b.conc_mass_comp_ref["substrate"]
                            )
                        )
                        * (
                            b.params.half_saturation_coeff_oxygen_heterotrophic
                            / (
                                b.params.half_saturation_coeff_oxygen_heterotrophic
                                + b.conc_mass_comp_ref["oxygen"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["nitrates"]
                            / (
                                b.params.half_saturation_coeff_nitrate_heterotrophic
                                + b.conc_mass_comp_ref["nitrates"]
                            )
                        )
                        * b.params.eta_g
                        * b.conc_mass_comp_ref["heterotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R3":
                    # R3: Aerobic growth of autotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.max_growth_rate_autotrophic
                        * (
                            b.conc_mass_comp_ref["ammonium"]
                            / (
                                b.params.half_saturation_coeff_ammonium
                                + b.conc_mass_comp_ref["ammonium"]
                            )
                        )
                        * (
                            b.conc_mass_comp_ref["oxygen"]
                            / (
                                b.params.half_saturation_coeff_oxygen_autotrophic
                                + b.conc_mass_comp_ref["oxygen"]
                            )
                        )
                        * b.conc_mass_comp_ref["autotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R4":
                    # R4: Decay of heterotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.decay_coeff_heterotrophic
                        * b.conc_mass_comp_ref["heterotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R5":
                    # R5: Decay of autotrophs
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.decay_coeff_autotrophic
                        * b.conc_mass_comp_ref["autotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R6":
                    # R6: Ammonification of soluble organic nitrogen
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.ammonification_rate
                        * b.conc_mass_comp_ref["nitrogen_soluble"]
                        * b.conc_mass_comp_ref["heterotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R7":
                    # R7: Hydrolysis of entrapped organics
                    return b.reaction_rate[r] == pyo.units.convert(
                        b.params.max_hydrolysis_rate
                        * (
                            b.conc_mass_comp_ref["biodegradable"]
                            / b.conc_mass_comp_ref["heterotrophic"]
                        )
                        / (
                            b.params.half_saturation_coeff_hydrolysis
                            + (
                                b.conc_mass_comp_ref["biodegradable"]
                                / b.conc_mass_comp_ref["heterotrophic"]
                            )
                        )
                        * (
                            (
                                b.conc_mass_comp_ref["oxygen"]
                                / (
                                    b.params.half_saturation_coeff_oxygen_heterotrophic
                                    + b.conc_mass_comp_ref["oxygen"]
                                )
                            )
                            + b.params.eta_h
                            * b.params.half_saturation_coeff_oxygen_heterotrophic
                            / (
                                b.params.half_saturation_coeff_oxygen_heterotrophic
                                + b.conc_mass_comp_ref["oxygen"]
                            )
                            * (
                                b.conc_mass_comp_ref["nitrates"]
                                / (
                                    b.params.half_saturation_coeff_nitrate_heterotrophic
                                    + b.conc_mass_comp_ref["nitrates"]
                                )
                            )
                        )
                        * b.conc_mass_comp_ref["heterotrophic"],
                        to_units=pyo.units.kg / pyo.units.m**3 / pyo.units.s,
                    )
                elif r == "R8":
                    # R8: Hydrolysis of entrapped organic nitrogen
                    return b.reaction_rate[r] == (
                        b.reaction_rate["R7"]
                        * (
                            b.conc_mass_comp_ref["nitrogen_particulate"]
                            / b.conc_mass_comp_ref["biodegradable"]
                        )
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
