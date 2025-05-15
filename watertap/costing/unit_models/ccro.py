#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from idaes.core.util.math import smooth_max
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_ccro_cost_param_block(blk):

    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.2,
        doc="Membrane replacement factor [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.membrane_cost = pyo.Var(
        initialize=30,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.high_pressure_membrane_cost = pyo.Var(
        initialize=75,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.pump_cost = pyo.Var(
        initialize=53 / 1e5 * 3600,
        bounds=(0, None),
        doc="High pressure pump cost",
        units=pyo.units.USD_2018 / pyo.units.watt,
    )


@register_costing_parameter_block(
    build_rule=build_ccro_cost_param_block,
    parameter_block_name="ccro",
)
def cost_ccro(blk, feed_pump=None, mp=None):
    """
    CCRO costing method
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    bs = mp.get_active_process_blocks()
    feed_pumps = [b.fs.P1 for b in bs]
    recirc_pumps = [b.fs.P2 for b in bs]
    feed_pump_work = feed_pumps[0].work_mechanical[0]
    recirc_pump_work = recirc_pumps[-1].work_mechanical[0]
    acc_time = bs[0].fs.dead_volume.accumulation_time[0]
    total_time = acc_time * len(bs)

    blk.feed_pumping_energy = [
        pyo.units.convert(
            b.fs.P1.work_mechanical[0] * acc_time,
            to_units=pyo.units.joule,
        )
        for b in bs
    ]
    blk.recirc_pumping_energy = [
        pyo.units.convert(
            b.fs.P2.work_mechanical[0] * acc_time,
            to_units=pyo.units.joule,
        )
        for b in bs
    ]

    blk.feed_pumping_power = pyo.Expression(
        expr=pyo.units.convert(
            sum(blk.feed_pumping_energy) / total_time,
            to_units=pyo.units.kilowatt,
        )
    )
    blk.recirc_pumping_power = pyo.Expression(
        expr=pyo.units.convert(
            sum(blk.recirc_pumping_energy) / total_time,
            to_units=pyo.units.kilowatt,
        )
    )
    blk.total_pumping_power = pyo.Expression(
        expr=pyo.units.convert(
            sum(blk.feed_pumping_energy + blk.recirc_pumping_energy) / total_time,
            to_units=pyo.units.kilowatt,
        )
    )

    blk.costing_package.add_cost_factor(blk, "TIC")

    blk.capital_cost_membrane = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of membrane",
    )

    blk.capital_cost_feed_pump = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of feed pump",
    )

    blk.capital_cost_recirculation_pump = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of recirculation pump",
    )

    blk.capital_cost_side_conduit = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of side conduit pressure vessel",
    )

    # blk.max_feed_pump_work = pyo.Var(
    #     initialize=1e5,
    #     domain=pyo.NonNegativeReals,
    #     units=pyo.units.kilowatt,
    #     doc="Maximum recirculation pump work",
    # )

    blk.max_recirculation_pump_work = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=pyo.units.watt,
        doc="Maximum recirculation pump work",
    )

    for t in mp.TIME:
        if t == 0:
            max_rpw = smooth_max(0, recirc_pumps[t].work_mechanical[0])
            # max_fpw = smooth_max(0, feed_pumps[t].work_mechanical[0])
        else:
            max_rpw = smooth_max(max_rpw, recirc_pumps[t].work_mechanical[0])
            # max_fpw = smooth_max(max_fpw, feed_pumps[t].work_mechanical[0])

    blk.max_recirculation_pump_work_constraint = pyo.Constraint(
        expr=blk.max_recirculation_pump_work == max_rpw
    )

    # blk.max_feed_pump_work_constraint = pyo.Constraint(
    #     expr=blk.max_feed_pump_work == max_fpw
    # )

    capital_cost_expr = 0

    blk.capital_cost_membrane_constraint = pyo.Constraint(
        expr=blk.capital_cost_membrane
        == pyo.units.convert(
            blk.costing_package.ccro.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_membrane

    blk.capital_cost_feed_pump_constraint = pyo.Constraint(
        expr=blk.capital_cost_feed_pump
        == pyo.units.convert(
            blk.costing_package.ccro.pump_cost
            * pyo.units.convert(feed_pump_work, to_units=pyo.units.watt),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_feed_pump

    blk.capital_cost_recirc_pump_constraint = pyo.Constraint(
        expr=blk.capital_cost_recirculation_pump
        == pyo.units.convert(
            blk.costing_package.ccro.pump_cost
            * pyo.units.convert(
                blk.max_recirculation_pump_work, to_units=pyo.units.watt
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_recirculation_pump

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.cost_factor * capital_cost_expr
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.costing_package.ccro.factor_membrane_replacement
            * blk.costing_package.ccro.membrane_cost
            * blk.unit_model.area,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.costing_package.cost_flow(blk.total_pumping_power, "electricity")
