##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

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

        del self.control_volume.enthalpy_balances

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
