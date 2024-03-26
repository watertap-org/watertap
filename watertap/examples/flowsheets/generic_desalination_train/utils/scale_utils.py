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

import math
from pyomo.environ import (
    value,
)
import logging

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "scale utils %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)

__author__ = "Alexander V. Dudchenko"


def calc_scale(value, factor=1):
    if value == 0:
        return 10 ** (-1 * math.log(abs(1), 10))
    else:
        return 10 ** (-1 * math.log(abs(value), 10) / factor)


def set_default_scaling(
    block_with_values,
    properties_block,
    component="flow_mol_phase_comp",
    factor=1,
    min_value=None,
    max_value=None,
    log_floor_scale=True,
    inverse_scale=False,
):
    for index in block_with_values.find_component(component):
        val = value(block_with_values.find_component(component)[index])
        if min_value != None and val < min_value:
            val = min_value
        if max_value != None and val > max_value:
            val = max_value
        if log_floor_scale:
            scale = calc_scale(val, factor=factor)
        elif inverse_scale:
            scale = 1 / abs(val)
        else:
            raise TypeError(
                "select type of scaling (log_floor_scale or inverse_scale must be True)"
            )
        _logger.info("Default scale for {}: {} val is: {}".format(index, scale, val))
        properties_block.set_default_scaling(component, scale, index=index)
