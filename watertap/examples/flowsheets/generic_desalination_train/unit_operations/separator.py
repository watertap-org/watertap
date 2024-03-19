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

from watertap.unit_models.generic_units.generic_separation import (
    GenericSeparation,
)
from pyomo.environ import (
    value,
    units as pyunits,
)
from watertap.examples.flowsheets.generic_desalination_train.costing import (
    separator_costing,
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


def build(m, block, base_cost=1, additive_cost=0, additive_dose=0):
    block.separator = GenericSeparation(property_package=m.fs.properties)
    separator_costing.cost_unit(
        m.fs.costing,
        block.separator,
        base_cost=base_cost,
        additive_cost=additive_cost,
        opt_name=block.process_name,
    )
    block.separator.additive_dose.fix(additive_dose)


def initialize(m, blk, solver):
    blk.separator.initialize(optarg=solver.options)


def display(m, blk):
    for (phase, ion), var in blk.separator.product_properties[
        0
    ].flow_mass_phase_comp.items():
        _logger.info(f"{ion} mass flow= {value(var)} kg/s")
    for (phase, ion), var in blk.separator.product_properties[
        0
    ].flow_mass_phase_comp.items():
        annual_ion_cost = value(
            blk.separator.separation_cost[ion]
            * pyunits.convert(
                blk.separator.product_properties[0].flow_mass_phase_comp["Liq", ion],
                to_units=pyunits.kg / pyunits.year,
            )
        )
        _logger.info(f"{ion} annual removal cost= {annual_ion_cost} $/year")
    annual_chem_cost = value(
        blk.separator.additive_cost
        * pyunits.convert(
            blk.separator.additive_mass_flow,
            to_units=pyunits.kg / pyunits.year,
        )
    )
    _logger.info(f"Annual chemical cost= {annual_chem_cost} $/year")
    _logger.info(f"Total annual cost= {value(blk.separator.annual_cost)} $/year")
