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

from watertap.unit_models.generic_desalter import (
    GenericDesalter,
)
from pyomo.environ import (
    value,
    Var,
    Constraint,
    Reals,
    units as pyunits,
)
from watertap.examples.flowsheets.generic_desalination_train.costing import (
    desalter_costing,
)
import logging
import idaes.core.util.scaling as iscale

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "desalter %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


__author__ = "Alexander V. Dudchenko"


def build(
    m,
    block,
    base_cost=1,
    recovery_cost=0,
    tracked_solids=None,
    min_recovery=0,
    max_recovery=99,
    default_recovery=80,
):
    block.desalter = GenericDesalter(property_package=m.fs.properties)

    block.desalter.water_recovery.setub(max_recovery)
    block.desalter.water_recovery.setlb(min_recovery)

    block.desalter.water_recovery.fix(default_recovery)
    desalter_costing.cost_desalter(
        m.fs.costing,
        block.desalter,
        base_cost,
        recovery_cost,
        opt_name=block.process_name,
    )

    block.desalter.brine_solids_concentration = Var(
        initialize=80,
        bounds=(None, None),
        domain=Reals,
        units=pyunits.kg / pyunits.m**3,
        doc="water recovery",
    )
    iscale.set_scaling_factor(block.desalter.brine_solids_concentration, 1)
    block.desalter.brine_solids_concentration.unfix()
    block.desalter.brine_water_percent = Var(
        initialize=80,
        bounds=(None, None),
        domain=Reals,
        units=pyunits.dimensionless,
        doc="water recovery",
    )
    iscale.set_scaling_factor(block.desalter.brine_water_percent, 1 / 100)
    block.desalter.brine_water_percent.unfix()
    block.desalter.brine_unit.properties_out[0].flow_vol_phase[...]

    block.desalter.solids_concentration_eq = Constraint(
        expr=block.desalter.brine_solids_concentration
        * block.desalter.brine_unit.properties_out[0].flow_vol_phase["Liq"]
        == sum(
            [
                block.desalter.brine_unit.properties_out[0].flow_mass_phase_comp[
                    "Liq", tds
                ]
                for tds in tracked_solids
            ]
        )
    )
    iscale.constraint_scaling_transform(block.desalter.solids_concentration_eq, 1)
    block.desalter.water_percent_eq = Constraint(
        expr=block.desalter.brine_water_percent
        * (
            sum(
                [
                    block.desalter.brine_unit.properties_out[0].flow_mass_phase_comp[
                        "Liq", tds
                    ]
                    for tds in tracked_solids
                ]
            )
            + block.desalter.brine_unit.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )
        == block.desalter.brine_unit.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ]
        * 100
    )
    iscale.constraint_scaling_transform(block.desalter.water_percent_eq, 1 / 100)


def initialize(m, blk, solver):
    blk.desalter.initialize(optarg=solver.options)


def unfix_opt_vars(m, blk):
    blk.desalter.water_recovery.unfix()


def display(m, blk):
    _logger.info(
        f"Feed flow {value(blk.desalter.brine_unit.properties_in[0].flow_mass_phase_comp['Liq','H2O'])}"
    )
    _logger.info(
        f"Brine flow {value(blk.desalter.brine_unit.properties_out[0].flow_mass_phase_comp['Liq','H2O'])}"
    )
    _logger.info(
        f"Product flow {value(blk.desalter.product_properties[0].flow_mass_phase_comp['Liq','H2O'])}"
    )
    _logger.info(f"Brine water content (%) {value(blk.desalter.brine_water_percent)}")
    _logger.info(f"Recovery (%) {value(blk.desalter.water_recovery)}")
    _logger.info(f"Annual cost ($) {value(blk.desalter.annual_cost)}")
