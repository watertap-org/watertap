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
Modified ASM2d reaction package for compatibility with modified ADM1.


Reference:

[1] Henze, M., Gujer, W., Mino, T., Matsuo, T., Wentzel, M.C., Marais, G.v.R.,
Van Loosdrecht, M.C.M., "Activated Sludge Model No.2D, ASM2D", 1999,
Wat. Sci. Tech. Vol. 39, No. 1, pp. 165-182
[2] Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J.,
"Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of
anaerobic digestion processes", 2016, Water Research, 95, pp.370-382.
"""

# Import IDAES cores
from idaes.core import declare_process_block_class
import idaes.logger as idaeslog
from watertap.property_models.unit_specific.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterData,
    _ASM2dReactionBlock,
    ASM2dReactionBlockData,
)

# Some more information about this module
__author__ = "Chenyu Wang"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SimpleModifiedASM2dReactionParameterBlock")
class SimpleModifiedASM2dReactionParameterData(ASM2dReactionParameterData):
    def build(self):
        super().build()

        self._reaction_block_class = SimpleModifiedASM2dReactionBlock

        if self.config.property_package.config.additional_solute_list is not None:
            for j in self.config.property_package.config.additional_solute_list:

                rate_reaction_stoichiometry_additional = {
                    ("R1", "Liq", j): 0,
                    ("R2", "Liq", j): 0,
                    ("R3", "Liq", j): 0,
                    ("R4", "Liq", j): 0,
                    ("R5", "Liq", j): 0,
                    ("R6", "Liq", j): 0,
                    ("R7", "Liq", j): 0,
                    ("R8", "Liq", j): 0,
                    ("R9", "Liq", j): 0,
                    ("R10", "Liq", j): 0,
                    ("R11", "Liq", j): 0,
                    ("R12", "Liq", j): 0,
                    ("R13", "Liq", j): 0,
                    ("R14", "Liq", j): 0,
                    ("R15", "Liq", j): 0,
                    ("R16", "Liq", j): 0,
                    ("R17", "Liq", j): 0,
                    ("R18", "Liq", j): 0,
                    ("R19", "Liq", j): 0,
                    ("R20", "Liq", j): 0,
                    ("R21", "Liq", j): 0,
                }
                new_rate_reaction_stoichiometry = {
                    **self.rate_reaction_stoichiometry,
                    **rate_reaction_stoichiometry_additional,
                }
                self.rate_reaction_stoichiometry = new_rate_reaction_stoichiometry


class _SimpleModifiedASM2dReactionBlock(_ASM2dReactionBlock):
    pass


@declare_process_block_class(
    "SimpleModifiedASM2dReactionBlock", block_class=_SimpleModifiedASM2dReactionBlock
)
class SimpleModifiedASM2dReactionBlockData(ASM2dReactionBlockData):
    pass
