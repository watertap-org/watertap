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

from pyomo.environ import Constraint
# Import IDAES cores
from idaes.generic_models.unit_models.mixer import MixerData
from idaes.core import declare_process_block_class
from pyomo.environ import Block

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

@declare_process_block_class("Mixer")
class MixerData(MixerData):
    """
    Standard Mixer Unit Model Class with get_costing method added in.
    """

    def build(self):
        super().build()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

    def get_costing(self, module=None, **kwargs):
        self.costing = Block()
        module.Mixer_costing(self.costing, **kwargs)

