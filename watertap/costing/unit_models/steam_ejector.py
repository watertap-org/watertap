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


def build_steam_ejector_cost_param_block(blk):
    """
    Build the costing parameter block for the steam ejector.
    """
    blk.base_cost = pyo.Var(
        initialize=1949,
        doc="Base cost coefficient for steam ejector",
        units=pyo.units.USD_2020,
    )
    blk.cost_exponent = pyo.Var(
        initialize=0.3,
        doc="Cost scaling exponent",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_steam_ejector_cost_param_block,
    parameter_block_name="steam_ejector",
)
def cost_steam_ejector(blk, cost_steam_flow=False):
    """
    Thermo Compressor (Steam Ejector) Costing Method.

    Capital cost is calculated using the equation:
        Capital Cost (USD) = 1949 × (S + EV)^0.3  (Gabriel 2015, Desalination)
    where:
        S  = Motive steam flow rate (kg/h)
        EV = Entrained vapor flow rate (kg/h)

    """

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    S = pyo.units.convert(
        blk.unit_model.properties_motive_steam[0.0].flow_mass_phase_comp["Vap", "H2O"],
        to_units=pyo.units.kg / pyo.units.hour,
    )

    EV = pyo.units.convert(
        blk.unit_model.properties_entrained_vapor[0.0].flow_mass_phase_comp[
            "Vap", "H2O"
        ],
        to_units=pyo.units.kg / pyo.units.hour,
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.steam_ejector.base_cost
            * (S + EV) ** blk.costing_package.steam_ejector.cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    if cost_steam_flow:

        pressure_sat = blk.unit_model.properties_motive_steam[0.0].pressure[0]
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
                (
                    blk.unit_model.properties_motive_steam[0.0].flow_mass_phase_comp[
                        "Vap", "H2O"
                    ]
                )
                * sp_vol,
                to_units=pyo.units.m**3 / pyo.units.s,
            ),
            "steam",
        )
