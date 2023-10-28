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


def build_chemical_precipitator_cost_param_block(blk):
    blk.capital_cost_param = Param(
        initialize=2000,
        units=pyunits.USD_2021 / (pyunits.kg / pyunits.day),
        mutable=True,
    )


@register_costing_parameter_block(
    build_rule=build_chemical_precipitator_cost_param_block,
    parameter_block_name="chemical_precipitator",
)
def cost_chemical_precipitator(blk):
    t0 = blk.flowsheet().time.first()
    make_capital_cost_var(blk)

    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            blk.costing_package.chemical_precipitator.capital_cost_param,
            to_units=blk.costing_package.base_currency / (pyunits.kg / pyunits.day),
        )
        * sum(
            pyunits.convert(
                obj,
                to_units=pyunits.kg / pyunits.day,
            )
            for reagent, obj in blk.unit_model.flow_mass_reagent.items()
        )
    )
