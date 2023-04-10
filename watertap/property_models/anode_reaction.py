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
Chlor-Alkali Anode Reaction Package
"""

from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    NonNegativeReals,
    log,
    Var,
    Param,
    Set,
    Suffix,
    value,
    check_optimal_termination,
    units as pyunits,
)
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
)
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog


__author__ = "Hunter Barber"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("AnodeReactionParameterBlock")
class AnodeReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super(AnodeReactionParameterData, self).build()

        self._reaction_block_class = ReactionBlock

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {
            ("R1", "Liq", "Cl-"): -1,
            ("R1", "Liq", "CL2-v"): 0.5,
        }

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({})
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _AnodeReactionBlock(ReactionBlockBase):
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


@declare_process_block_class("ReactionBlock", block_class=_AnodeReactionBlock)
class AnodeReactionBlockData(ReactionBlockDataBase):
    """
    Anode reaction package
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(AnodeReactionBlockData, self).build()

        # Create references to state vars

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar
