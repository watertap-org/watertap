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

import functools

import pyomo.environ as pyo


def register_costing_parameter_block(build_rule, parameter_block_name):
    def register_costing_parameter_block_decorator(func):
        @functools.wraps(func)
        def add_costing_parameter_block(blk, *args, **kwargs):
            parameter_block = blk.costing_package.component(parameter_block_name)
            if parameter_block is None:
                parameter_block = pyo.Block(rule=build_rule)
                blk.costing_package.add_component(parameter_block_name, parameter_block)
                # fix the parameters in case the build_rule did not
                parameter_block.fix_all_vars()
            elif parameter_block._rule is None or not hasattr(
                parameter_block._rule, "_fcn"
            ):
                raise RuntimeError(
                    "Use the register_costing_parameter_block decorator for specifying"
                    "costing-package-level parameters"
                )
            elif parameter_block._rule._fcn is not build_rule:
                other_rule = parameter_block._rule._fcn
                raise RuntimeError(
                    "Attempting to add identically named costing parameter blocks with "
                    "different build rules to the costing package "
                    f"{blk.costing_package}. Parameter block named "
                    f"{parameter_block_name} was previously built by function "
                    f"{other_rule.__name__} from module {other_rule.__module__}"
                )
            # else parameter_block was constructed by build_rule previously
            return func(blk, *args, **kwargs)

        return add_costing_parameter_block

    return register_costing_parameter_block_decorator


def make_capital_cost_var(blk):
    blk.capital_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit capital cost",
    )


def make_fixed_operating_cost_var(blk):
    blk.fixed_operating_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit fixed operating cost",
    )


def cost_membrane(blk, membrane_cost, factor_membrane_replacement):
    """
    Generic function for costing a membrane. Assumes the unit_model
    has an `area` variable or parameter.

    Args:
        membrane_cost - The cost of the membrane in currency per area
        factor_membrane_replacement - Membrane replacement factor
                                      [fraction of membrane replaced/year]
    """

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.membrane_cost = pyo.Expression(expr=membrane_cost)
    blk.factor_membrane_replacement = pyo.Expression(expr=factor_membrane_replacement)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.factor_membrane_replacement * blk.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )


def cost_by_flow_volume(blk, flow_cost, flow_to_cost):
    """
    Generic function for costing by flow volume.

    Args:
        flow_cost - The cost of the device in [currency]/([volume]/[time])
        flow_to_cost - The flow costed in [volume]/[time]
    """
    make_capital_cost_var(blk)
    blk.flow_cost = pyo.Expression(expr=flow_cost)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.flow_cost * flow_to_cost, to_units=blk.costing_package.base_currency
        )
    )
