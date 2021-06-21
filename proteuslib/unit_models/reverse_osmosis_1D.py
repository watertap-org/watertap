###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################


from enum import Enum
# Import Pyomo libraries
from pyomo.environ import (Var,
                           Set,
                           Param,
                           SolverFactory,
                           Suffix,
                           NonNegativeReals,
                           NegativeReals,
                           Reference,
                           Block,
                           units as pyunits,
                           exp,
                           value)
from pyomo.common.config import ConfigBlock, ConfigValue, In
# Import IDAES cores
from idaes.core import (ControlVolume1DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault,
                        FlowDirection)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import get_solver, scaling as iscale

import idaes.logger as idaeslog


__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis1D")
class ReverseOsmosis1DData(UnitModelBlockData):
    """Standard 1D Reverse Osmosis Unit Model Class."""

    CONFIG = UnitModelBlockData.CONFIG()

