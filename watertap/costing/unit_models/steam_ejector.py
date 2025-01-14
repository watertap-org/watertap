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

    blk.steam_cost = pyo.Var(
        initialize=0.008,
        units=pyo.units.USD_2018 / pyo.units.kg,
        doc="Steam cost per kg",
    )

    blk.parent_block().register_flow_type("steam", blk.steam_cost)


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
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.properties_motive_steam[0.0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ],
                to_units=pyo.units.kg / pyo.units.s,
            ),
            "steam",
        )
