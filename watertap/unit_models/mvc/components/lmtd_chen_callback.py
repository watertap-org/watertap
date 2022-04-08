###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from pyomo.environ import units as pyunits, ExternalFunction
from idaes.core.util.functions import functions_lib


def delta_temperature_chen_callback(b):
    r"""
    This is a callback for a temperature difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD) approximation given by Chen (1987).  It can be
    supplied to "delta_temperature_callback" HeatExchanger configuration option.
    This uses a cube root function that works with negative numbers returning
    the real negative root. This should always evaluate successfully.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    temp_units = pyunits.get_units(dT1)

    # external function that ruturns the real root, for the cuberoot of negitive
    # numbers, so it will return without error for positive and negitive dT.
    b.cbrt = ExternalFunction(
        library=functions_lib(), function="cbrt", arg_units=[temp_units**3]
    )

    @b.Expression(b.flowsheet().time)
    def delta_temperature(b, t):
        return b.cbrt(dT1[t] * dT2[t] * 0.5 * (dT1[t] + dT2[t])) * temp_units
