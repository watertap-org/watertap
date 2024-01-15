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

# Import Pyomo libraries
from pyomo.environ import (
    Param,
    Constraint,
    units as pyunits,
)
from ..util import register_costing_parameter_block, make_capital_cost_var


def build_stoichiometric_reactor_cost_param_block(blk):
    blk.capital_cost_softening = Param(
        initialize=374.9,
        units=pyunits.USD_2021 / (pyunits.lb / pyunits.day),
        mutable=True,
        doc="Cost for typical mid sized softening reactor",
    )

    blk.capital_cost_acid_addition = Param(
        initialize=127.8,
        units=pyunits.USD_2021 / (pyunits.gallon / pyunits.day),
        doc="Cost for acid mixer",
    )


@register_costing_parameter_block(
    build_rule=build_stoichiometric_reactor_cost_param_block,
    parameter_block_name="stoichiometric_reactor",
)
def cost_stoichiometric_reactor(blk):
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    if (
        blk.unit_model.has_precipitation_reaction
        and blk.unit_model.has_dissolution_reaction
    ):
        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor
            * pyunits.convert(
                blk.costing_package.stoichiometric_reactor.capital_cost_softening,
                to_units=blk.costing_package.base_currency / (pyunits.lb / pyunits.day),
            )
            * sum(
                pyunits.convert(
                    obj,
                    to_units=pyunits.lb / pyunits.day,
                )
                for reagent, obj in blk.unit_model.flow_mass_reagent.items()
            ),
        )
    elif blk.unit_model.has_dissolution_reaction:
        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.cost_factor
            * pyunits.convert(
                blk.costing_package.stoichiometric_reactor.capital_cost_acid_addition,
                to_units=blk.costing_package.base_currency
                / (pyunits.gallon / pyunits.day),
            )
            * sum(
                pyunits.convert(
                    obj,
                    to_units=pyunits.gallon / pyunits.day,
                )
                for reagent, obj in blk.unit_model.flow_vol_reagent.items()
            )
        )
