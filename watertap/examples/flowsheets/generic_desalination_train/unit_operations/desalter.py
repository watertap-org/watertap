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

from watertap.unit_models.generic_units.generic_desalter import (
    GenericDesalter,
)
from pyomo.environ import (
    value,
)
from watertap.examples.flowsheets.generic_desalination_train.costing import (
    desalter_costing,
)
import logging

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "desalter %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


__author__ = "Alexander V. Dudchenko"


def build(m, block, base_cost=1, recovery_cost=0):
    block.desalter = GenericDesalter(property_package=m.fs.properties)
    desalter_costing.cost_desalter(
        m.fs.costing,
        block.desalter,
        base_cost,
        recovery_cost,
        opt_name=block.process_name,
    )


def initialize(m, blk, solver):
    blk.desalter.initialize(optarg=solver.options)


def unfix_opt_vars(m, blk):
    blk.desalter.water_recovery.unfix()


def display(m, blk):
    _logger.info(
        f"Feed flow {value(blk.desalter.brine_unit.properties_in[0].flow_vol_phase['Liq'])}"
    )
    _logger.info(
        f"Brine flow {value(blk.desalter.brine_unit.properties_out[0].flow_vol_phase['Liq'])}"
    )
    _logger.info(
        f"Product flow {value(blk.desalter.product_properties[0].flow_vol_phase['Liq'])}"
    )
    _logger.info(f"Recovery (%) {value(blk.desalter.water_recovery)}")
    _logger.info(f"Annual cost ($) {value(blk.desalter.annual_cost)}")
