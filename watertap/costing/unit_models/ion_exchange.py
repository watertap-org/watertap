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
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_hcl_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0.17,
        doc="HCl cost",  # for 37% sol'n - CatCost v 1.0.4
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Param(
        mutable=True,
        initialize=0.37,
        doc="HCl purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("HCl", blk.cost / blk.purity)


def build_naoh_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0.59,
        doc="NaOH cost",  # for 30% sol'n - iDST
        units=pyo.units.USD_2020 / pyo.units.kg,
    )

    blk.purity = pyo.Param(
        mutable=True,
        initialize=0.30,
        doc="NaOH purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("NaOH", blk.cost / blk.purity)


def build_meoh_cost_param_block(blk):
    # MeOH = Methanol
    blk.cost = pyo.Param(
        mutable=True,
        initialize=3.395,
        doc="MeOH cost",  # for 100% purity - ICIS
        units=pyo.units.USD_2008 / pyo.units.kg,
    )

    blk.purity = pyo.Param(
        mutable=True,
        initialize=1,
        doc="MeOH purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("MeOH", blk.cost / blk.purity)


def build_nacl_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0.09,
        doc="NaCl cost",  # for solid, 100% purity - CatCost
        units=pyo.units.USD_2020 / pyo.units.kg,
    )

    blk.purity = pyo.Param(
        mutable=True,
        initialize=1,
        doc="NaCl purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("NaCl", blk.cost / blk.purity)


def build_ion_exhange_cost_param_block(blk):
    blk.anion_exchange_resin_cost = pyo.Var(
        initialize=205,
        units=pyo.units.USD_2020 / pyo.units.ft**3,
        doc="Anion exchange resin cost per cubic ft. Assumes strong base polystyrenic gel-type Type II. From EPA-WBS cost model.",
    )
    blk.cation_exchange_resin_cost = pyo.Var(
        initialize=153,
        units=pyo.units.USD_2020 / pyo.units.ft**3,
        doc="Cation exchange resin cost per cubic ft. Assumes strong acid polystyrenic gel-type. From EPA-WBS cost model.",
    )
    # Ion exchange pressure vessels costed with power equation, col_vol in gallons:
    #   pressure_vessel_cost = A * col_vol ** b

    blk.vessel_A_coeff = pyo.Var(
        initialize=1596.499333,
        units=pyo.units.USD_2020,
        doc="Ion exchange pressure vessel cost equation - A coeff., Carbon steel w/ stainless steel internals",
    )
    blk.vessel_b_coeff = pyo.Var(
        initialize=0.459496809,
        units=pyo.units.dimensionless,
        doc="Ion exchange pressure vessel cost equation - b coeff., Carbon steel w/ stainless steel internals",
    )

    # Ion exchange backwash/rinse tank costed with power equation, tank_vol in gallons:
    #   bw_tank_cost = A * tank_vol ** b

    blk.backwash_tank_A_coeff = pyo.Var(
        initialize=308.9371309,
        units=pyo.units.USD_2020,
        doc="Ion exchange backwash tank cost equation - A coeff., Steel tank",
    )
    blk.backwash_tank_b_coeff = pyo.Var(
        initialize=0.501467571,
        units=pyo.units.dimensionless,
        doc="Ion exchange backwash tank cost equation - b coeff., Steel tank",
    )
    # Ion exchange regeneration solution tank costed with power equation, tank_vol in gallons:
    #   regen_tank_cost = A * tank_vol ** b

    blk.regen_tank_A_coeff = pyo.Var(
        initialize=57.02158923,
        units=pyo.units.USD_2020,
        doc="Ion exchange regen tank cost equation - A coeff. Stainless steel",
    )
    blk.regen_tank_b_coeff = pyo.Var(
        initialize=0.729325391,
        units=pyo.units.dimensionless,
        doc="Ion exchange regen tank cost equation - b coeff. Stainless steel",
    )
    blk.annual_resin_replacement_factor = pyo.Var(
        initialize=0.05,
        units=pyo.units.year**-1,
        doc="Fraction of ion excange resin replaced per year, 4-5% of bed volume - EPA",
    )
    blk.hazardous_min_cost = pyo.Var(
        initialize=3240,
        units=pyo.units.USD_2020 / pyo.units.year,
        doc="Min cost per hazardous waste shipment - EPA",
    )
    blk.hazardous_resin_disposal = pyo.Var(
        initialize=347.10,
        units=pyo.units.USD_2020 * pyo.units.ton**-1,
        doc="Hazardous resin disposal cost - EPA",
    )
    blk.hazardous_regen_disposal = pyo.Var(
        initialize=3.64,
        units=pyo.units.USD_2020 * pyo.units.gal**-1,
        doc="Hazardous liquid disposal cost - EPA",
    )
    blk.regen_recycle = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        doc="Number of cycles the regenerant can be reused before disposal",
    )


@register_costing_parameter_block(
    build_rule=build_hcl_cost_param_block,
    parameter_block_name="hcl",
)
@register_costing_parameter_block(
    build_rule=build_naoh_cost_param_block,
    parameter_block_name="naoh",
)
@register_costing_parameter_block(
    build_rule=build_meoh_cost_param_block,
    parameter_block_name="meoh",
)
@register_costing_parameter_block(
    build_rule=build_nacl_cost_param_block,
    parameter_block_name="nacl",
)
@register_costing_parameter_block(
    build_rule=build_ion_exhange_cost_param_block,
    parameter_block_name="ion_exchange",
)
def cost_ion_exchange(blk):
    """
    Volume-based capital cost for Ion Exchange
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    ion_exchange_params = blk.costing_package.ion_exchange
    # Conversions to use units from cost equations in reference
    tot_num_col = blk.unit_model.number_columns + blk.unit_model.number_columns_redund
    col_vol_gal = pyo.units.convert(blk.unit_model.col_vol_per, to_units=pyo.units.gal)
    bed_vol_ft3 = pyo.units.convert(blk.unit_model.bed_vol, to_units=pyo.units.ft**3)

    ix_type = blk.unit_model.ion_exchange_type
    blk.regen_soln_dens = pyo.Param(
        initialize=1000,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Density of regeneration solution",
    )
    blk.regen_dose = pyo.Param(
        initialize=300,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]",
    )
    blk.capital_cost_vessel = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost for one vessel",
    )
    blk.capital_cost_resin = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost for resin for one vessel",
    )
    blk.capital_cost_regen_tank = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost for regeneration solution tank",
    )
    blk.capital_cost_backwash_tank = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost for backwash + rinse solution tank",
    )
    blk.operating_cost_hazardous = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Operating cost for hazardous waste disposal",
    )
    blk.flow_mass_regen_soln = pyo.Var(
        initialize=1,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kg / pyo.units.year,
        doc="Regeneration solution flow",
    )
    blk.total_pumping_power = pyo.Var(
        initialize=1,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kilowatt,
        doc="Total pumping power required",
    )

    if ix_type == "cation":
        resin_cost = ion_exchange_params.cation_exchange_resin_cost

    elif ix_type == "anion":
        resin_cost = ion_exchange_params.anion_exchange_resin_cost

    blk.capital_cost_vessel_constraint = pyo.Constraint(
        expr=blk.capital_cost_vessel
        == pyo.units.convert(
            (
                ion_exchange_params.vessel_A_coeff
                * (col_vol_gal / pyo.units.gallon) ** ion_exchange_params.vessel_b_coeff
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.capital_cost_resin_constraint = pyo.Constraint(
        expr=blk.capital_cost_resin
        == pyo.units.convert(
            resin_cost * bed_vol_ft3, to_units=blk.costing_package.base_currency
        )
    )
    if blk.unit_model.config.regenerant == "single_use":
        blk.capital_cost_regen_tank.fix(0)
        blk.flow_mass_regen_soln.fix(0)
        blk.flow_vol_resin = pyo.Var(
            initialize=1e5,
            bounds=(0, None),
            units=pyo.units.m**3 / blk.costing_package.base_period,
            doc="Volumetric flow of resin per cycle",  # assumes you are only replacing the operational columns, t_cycle = t_breakthru
        )
        blk.single_use_resin_replacement_cost = pyo.Var(
            initialize=1e5,
            bounds=(0, None),
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Operating cost for using single-use resin (i.e., no regeneration)",
        )

        blk.flow_vol_resin_constraint = pyo.Constraint(
            expr=blk.flow_vol_resin
            == pyo.units.convert(
                blk.unit_model.bed_vol_tot / blk.unit_model.t_breakthru,
                to_units=pyo.units.m**3 / blk.costing_package.base_period,
            )
        )
        blk.mass_flow_resin = pyo.units.convert(
            blk.flow_vol_resin * blk.unit_model.resin_bulk_dens,
            to_units=pyo.units.ton / blk.costing_package.base_period,
        )
    else:

        blk.regeneration_tank_vol = pyo.Expression(
            expr=pyo.units.convert(
                blk.unit_model.regen_tank_vol,
                to_units=pyo.units.gal,
            )
        )

        blk.capital_cost_regen_tank_constraint = pyo.Constraint(
            expr=blk.capital_cost_regen_tank
            == pyo.units.convert(
                ion_exchange_params.regen_tank_A_coeff
                * (blk.regeneration_tank_vol / pyo.units.gallon)
                ** ion_exchange_params.regen_tank_b_coeff,
                to_units=blk.costing_package.base_currency,
            )
        )

    blk.backwash_tank_vol = pyo.Expression(
        expr=pyo.units.convert(
            (
                blk.unit_model.bw_flow * blk.unit_model.t_bw
                + blk.unit_model.rinse_flow * blk.unit_model.t_rinse
            ),
            to_units=pyo.units.gal,
        )
    )

    blk.capital_cost_backwash_tank_constraint = pyo.Constraint(
        expr=blk.capital_cost_backwash_tank
        == pyo.units.convert(
            ion_exchange_params.backwash_tank_A_coeff
            * (blk.backwash_tank_vol / pyo.units.gallon)
            ** ion_exchange_params.backwash_tank_b_coeff,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            (
                ((blk.capital_cost_vessel + blk.capital_cost_resin) * tot_num_col)
                + blk.capital_cost_backwash_tank
                + blk.capital_cost_regen_tank
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    if blk.unit_model.config.hazardous_waste:

        if blk.unit_model.config.regenerant == "single_use":
            blk.operating_cost_hazardous_constraint = pyo.Constraint(
                expr=blk.operating_cost_hazardous
                == pyo.units.convert(
                    (blk.mass_flow_resin * ion_exchange_params.hazardous_resin_disposal)
                    + ion_exchange_params.hazardous_min_cost,
                    to_units=blk.costing_package.base_currency
                    / blk.costing_package.base_period,
                )
            )
        else:
            bed_mass_ton = pyo.units.convert(
                blk.unit_model.bed_vol * blk.unit_model.resin_bulk_dens,
                to_units=pyo.units.ton,
            )
            blk.operating_cost_hazardous_constraint = pyo.Constraint(
                expr=blk.operating_cost_hazardous
                == pyo.units.convert(
                    (
                        bed_mass_ton
                        * tot_num_col
                        * ion_exchange_params.hazardous_resin_disposal
                    )
                    * ion_exchange_params.annual_resin_replacement_factor
                    + pyo.units.convert(
                        blk.flow_mass_regen_soln / blk.regen_soln_dens,
                        to_units=pyo.units.gal / pyo.units.year,
                    )
                    * ion_exchange_params.hazardous_regen_disposal
                    + ion_exchange_params.hazardous_min_cost,
                    to_units=blk.costing_package.base_currency
                    / blk.costing_package.base_period,
                )
            )
    else:
        blk.operating_cost_hazardous.fix(0)
    if blk.unit_model.config.regenerant == "single_use":

        blk.single_use_resin_replacement_cost_constraint = pyo.Constraint(
            expr=blk.single_use_resin_replacement_cost
            == pyo.units.convert(
                blk.flow_vol_resin * resin_cost,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == blk.single_use_resin_replacement_cost + blk.operating_cost_hazardous
        )

    else:
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                (
                    (
                        bed_vol_ft3
                        * tot_num_col
                        * ion_exchange_params.annual_resin_replacement_factor
                        * resin_cost
                    )
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
            + blk.operating_cost_hazardous
        )

        blk.flow_mass_regen_soln_constraint = pyo.Constraint(
            expr=blk.flow_mass_regen_soln
            == pyo.units.convert(
                (
                    (blk.regen_dose * blk.unit_model.bed_vol * tot_num_col)
                    / (blk.unit_model.t_cycle)
                )
                / ion_exchange_params.regen_recycle,
                to_units=pyo.units.kg / pyo.units.year,
            )
        )

        blk.costing_package.cost_flow(
            blk.flow_mass_regen_soln, blk.unit_model.config.regenerant
        )

    power_expr = (
        blk.unit_model.main_pump_power
        + blk.unit_model.bw_pump_power
        + blk.unit_model.rinse_pump_power
    )
    if blk.unit_model.config.regenerant != "single_use":
        power_expr += blk.unit_model.regen_pump_power

    blk.total_pumping_power_constr = pyo.Constraint(
        expr=blk.total_pumping_power == power_expr
    )

    blk.costing_package.cost_flow(blk.total_pumping_power, "electricity")
