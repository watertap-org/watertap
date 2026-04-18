#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from ..util import register_costing_parameter_block, make_capital_cost_var


def build_heat_exchanger_cost_param_block(blk):

    blk.unit_cost = pyo.Var(
        initialize=300,
        doc="Estimated from multiple sources",
        units=pyo.units.USD_2020,
    )

    blk.material_factor_cost = pyo.Var(
        initialize=1,
        doc="Material factor",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_heat_exchanger_cost_param_block,
    parameter_block_name="heat_exchanger",
)
def cost_heat_exchanger(blk, cost_steam_flow=False):
    """
    Heat Exchanger Costing Method

    TODO: describe equations
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.heat_exchanger.unit_cost
            * blk.costing_package.heat_exchanger.material_factor_cost
            * (
                pyo.units.convert(blk.unit_model.area, to_units=(pyo.units.m**2))
                / pyo.units.m**2
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    if cost_steam_flow:

        pressure_sat = blk.unit_model.hot_side_inlet.pressure[0]
        # 1. Compute saturation temperature of steam: computed from El-Dessouky expression
        tsat_constants = [
            42.6776 * pyo.units.K,
            -3892.7 * pyo.units.K,
            1000 * pyo.units.kPa,
            -9.48654 * pyo.units.dimensionless,
        ]
        psat = (
            pyo.units.convert(pressure_sat, to_units=pyo.units.kPa)
            + 101.325 * pyo.units.kPa
        )
        temperature_sat = tsat_constants[0] + tsat_constants[1] / (
            pyo.log(psat / tsat_constants[2]) + tsat_constants[3]
        )
        # 3. Compute specific volume: computed from Affandi expression (Eq 5)
        t_critical = 647.096 * pyo.units.K
        t_red = temperature_sat / t_critical  # Reduced temperature
        sp_vol_constants = [
            -7.75883 * pyo.units.dimensionless,
            3.23753 * pyo.units.dimensionless,
            2.05755 * pyo.units.dimensionless,
            -0.06052 * pyo.units.dimensionless,
            0.00529 * pyo.units.dimensionless,
        ]

        log_sp_vol = (
            sp_vol_constants[0]
            + sp_vol_constants[1] * (pyo.log(1 / t_red)) ** 0.4
            + sp_vol_constants[2] / (t_red**2)
            + sp_vol_constants[3] / (t_red**4)
            + sp_vol_constants[4] / (t_red**5)
        )
        sp_vol = pyo.exp(log_sp_vol) * pyo.units.m**3 / pyo.units.kg

        blk.costing_package.cost_flow(
            pyo.units.convert(
                (blk.unit_model.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"])
                * sp_vol,
                to_units=pyo.units.m**3 / pyo.units.s,
            ),
            "steam",
        )
