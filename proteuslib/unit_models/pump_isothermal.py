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
from idaes.generic_models.unit_models.pressure_changer import PumpData
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

@declare_process_block_class("Pump")
class PumpIsothermalData(PumpData):
    """
    Standard Isothermal Pump Unit Model Class
    """

    def build(self):
        super().build()

        self.control_volume.del_component(self.control_volume.enthalpy_balances)

        @self.control_volume.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal constraint")
        def isothermal_balance(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for ind, c in self.control_volume.isothermal_balance.items():
            sf = iscale.get_scaling_factor(self.control_volume.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)
