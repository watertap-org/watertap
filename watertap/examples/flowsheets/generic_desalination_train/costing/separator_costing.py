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


from pyomo.environ import (
    Var,
    Constraint,
    units as pyunits,
)

__author__ = "Alexander V. Dudchenko"


def cost_unit(
    costing_block,
    separator_unit,
    base_cost=1,
    additive_cost=0,
    species_cost=0,
    opt_name=None,
):
    separator_unit.base_cost = Var(
        initialize=base_cost, units=costing_block.base_currency / pyunits.m**3
    )
    separator_unit.additive_cost = Var(
        initialize=additive_cost, units=costing_block.base_currency / pyunits.kg
    )

    separator_unit.separation_cost = Var(
        list(separator_unit.component_removal_percent.keys()),
        initialize=species_cost,
        units=costing_block.base_currency / pyunits.kg,
    )
    separator_unit.product_properties[0].flow_mass_phase_comp[...]
    separator_unit.base_cost.fix()
    separator_unit.separation_cost.fix()
    separator_unit.additive_cost.fix()
    separator_unit.annual_cost = Var(
        initialize=1, units=costing_block.base_currency / pyunits.year
    )
    separator_unit.absolute_cost_eq = Constraint(
        expr=separator_unit.annual_cost
        == separator_unit.base_cost
        * pyunits.convert(
            separator_unit.separator_unit.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.year,
        )
        + separator_unit.additive_cost
        * pyunits.convert(
            separator_unit.additive_mass_flow,
            to_units=pyunits.kg / pyunits.year,
        )
        + sum(
            [
                separator_unit.separation_cost[ion]
                * pyunits.convert(
                    separator_unit.product_properties[0].flow_mass_phase_comp[
                        "Liq", ion
                    ],
                    to_units=pyunits.kg / pyunits.year,
                )
                for ion in separator_unit.separation_cost.keys()
            ]
        )
    )
    costing_block.annual_costs.append(
        {"name": opt_name, "annual_cost": separator_unit.annual_cost}
    )
