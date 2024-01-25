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

import pyomo.environ as pyo
from ..util import register_costing_parameter_block, make_capital_cost_var


def build_compressor_cost_param_block(blk):
    blk.unit_cost = pyo.Var(
        initialize=7364,
        doc="Compressor unit cost (El-Sayed et al., 2001)",
        units=pyo.units.USD_2001,
    )

    blk.exponent = pyo.Var(
        initialize=0.7,
        doc="Compressor exponent (El-Sayed et al., 2001)",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_compressor_cost_param_block,
    parameter_block_name="compressor",
)
def cost_compressor(blk, cost_electricity_flow=True):
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.compressor.unit_cost
            * blk.unit_model.control_volume.properties_in[0].flow_mass_phase_comp[
                "Vap", "H2O"
            ]  # Add a convert to kg/s
            / (pyo.units.kg / pyo.units.s)
            * blk.unit_model.pressure_ratio
            * (blk.unit_model.efficiency / (1 - blk.unit_model.efficiency))
            ** blk.costing_package.compressor.exponent,
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.control_volume.work[0], to_units=pyo.units.kW
            ),
            "electricity",
        )
