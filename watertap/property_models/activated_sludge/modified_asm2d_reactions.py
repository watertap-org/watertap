#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterData,
    _ASM2dReactionBlock,
    ASM2dReactionBlockData,
)

# Some more information about this module
__author__ = "Chenyu Wang"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ModifiedASM2dReactionParameterBlock")
class ModifiedASM2dReactionParameterData(ASM2dReactionParameterData):
    def build(self):
        super().build()

        self._reaction_block_class = ModifiedASM2dReactionBlock

        rate_reaction_stoichiometry_K_Mg = {
            ("R1", "Liq", "S_K"): 0,
            ("R1", "Liq", "S_Mg"): 0,
            ("R2", "Liq", "S_K"): 0,
            ("R2", "Liq", "S_Mg"): 0,
            ("R3", "Liq", "S_K"): 0,
            ("R3", "Liq", "S_Mg"): 0,
            ("R4", "Liq", "S_K"): 0,
            ("R4", "Liq", "S_Mg"): 0,
            ("R5", "Liq", "S_K"): 0,
            ("R5", "Liq", "S_Mg"): 0,
            ("R6", "Liq", "S_K"): 0,
            ("R6", "Liq", "S_Mg"): 0,
            ("R7", "Liq", "S_K"): 0,
            ("R7", "Liq", "S_Mg"): 0,
            ("R8", "Liq", "S_K"): 0,
            ("R8", "Liq", "S_Mg"): 0,
            ("R9", "Liq", "S_K"): 0,
            ("R9", "Liq", "S_Mg"): 0,
            ("R10", "Liq", "S_K"): 0,
            ("R10", "Liq", "S_Mg"): 0,
            ("R11", "Liq", "S_K"): 0,
            ("R11", "Liq", "S_Mg"): 0,
            ("R12", "Liq", "S_K"): 0,
            ("R12", "Liq", "S_Mg"): 0,
            ("R13", "Liq", "S_K"): 0,
            ("R13", "Liq", "S_Mg"): 0,
            ("R14", "Liq", "S_K"): 0,
            ("R14", "Liq", "S_Mg"): 0,
            ("R15", "Liq", "S_K"): 0,
            ("R15", "Liq", "S_Mg"): 0,
            ("R16", "Liq", "S_K"): 0,
            ("R16", "Liq", "S_Mg"): 0,
            ("R17", "Liq", "S_K"): 0,
            ("R17", "Liq", "S_Mg"): 0,
            ("R18", "Liq", "S_K"): 0,
            ("R18", "Liq", "S_Mg"): 0,
            ("R19", "Liq", "S_K"): 0,
            ("R19", "Liq", "S_Mg"): 0,
            ("R20", "Liq", "S_K"): 0,
            ("R20", "Liq", "S_Mg"): 0,
            ("R21", "Liq", "S_K"): 0,
            ("R21", "Liq", "S_Mg"): 0,
        }
        new_rate_reaction_stoichiometry = {
            **self.rate_reaction_stoichiometry,
            **rate_reaction_stoichiometry_K_Mg,
        }
        self.rate_reaction_stoichiometry = new_rate_reaction_stoichiometry


class _ModifiedASM2dReactionBlock(_ASM2dReactionBlock):
    pass


@declare_process_block_class(
    "ModifiedASM2dReactionBlock", block_class=_ModifiedASM2dReactionBlock
)
class ModifiedASM2dReactionBlockData(ASM2dReactionBlockData):
    pass
