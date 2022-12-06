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

import pyomo.environ as pyo
from idaes.core.util.math import smooth_min
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_gac_cost_param_block(blk):

    blk.num_contactors_op = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        doc="Number of GAC contactors in operation in parallel",
    )

    blk.num_contactors_redundant = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        doc="Number of off-line redundant GAC contactors in parallel",
    )

    blk.contactor_cost_coeff_0 = pyo.Var(
        initialize=10010.9,
        units=pyo.units.USD_2020,
        doc="GAC contactor polynomial cost coefficient 0",
    )

    blk.contactor_cost_coeff_1 = pyo.Var(
        initialize=2204.95,
        units=pyo.units.USD_2020 * (pyo.units.m**3) ** -1,
        doc="GAC contactor polynomial cost coefficient 1",
    )

    blk.contactor_cost_coeff_2 = pyo.Var(
        initialize=-15.9378,
        units=pyo.units.USD_2020 * (pyo.units.m**3) ** -2,
        doc="GAC contactor polynomial cost coefficient 2",
    )

    blk.contactor_cost_coeff_3 = pyo.Var(
        initialize=0.110592,
        units=pyo.units.USD_2020 * (pyo.units.m**3) ** -3,
        doc="GAC contactor polynomial cost coefficient 3",
    )

    blk.bed_mass_max_ref = pyo.Var(
        initialize=18143.7,
        units=pyo.units.kg,
        doc="Reference maximum value of GAC mass needed for initial charge where "
        "economy of scale no longer discounts the unit price",
    )

    blk.adsorbent_unit_cost_coeff = pyo.Var(
        initialize=4.58342,
        units=pyo.units.USD_2020 * pyo.units.kg**-1,
        doc="GAC adsorbent exponential cost pre-exponential coefficient",
    )

    blk.adsorbent_unit_cost_exp_coeff = pyo.Var(
        initialize=-1.25311e-5,
        units=pyo.units.kg**-1,
        doc="GAC adsorbent exponential cost parameter coefficient",
    )

    blk.other_cost_coeff = pyo.Var(
        initialize=16660.7,
        units=pyo.units.USD_2020,
        doc="GAC other cost power law coefficient",
    )

    blk.other_cost_exp = pyo.Var(
        initialize=0.552207,
        units=pyo.units.dimensionless,
        doc="GAC other cost power law exponent",
    )

    blk.regen_frac = pyo.Var(
        initialize=0.70,
        units=pyo.units.dimensionless,
        doc="Fraction of spent GAC adsorbent that can be regenerated for reuse",
    )

    blk.regen_unit_cost = pyo.Var(
        initialize=4.28352,
        units=pyo.units.USD_2020 * pyo.units.kg**-1,
        doc="Unit cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
    )

    blk.makeup_unit_cost = pyo.Var(
        initialize=4.58223,
        units=pyo.units.USD_2020 * pyo.units.kg**-1,
        doc="Unit cost to makeup spent GAC adsorbent with fresh adsorbent",
    )


@register_costing_parameter_block(
    build_rule=build_gac_cost_param_block,
    parameter_block_name="gac",
)
def cost_gac(blk):
    """
    3 equation capital cost estimation for GAC systems with: (i), contactor/pressure vessel cost by polynomial
    as a function of individual contactor volume; (ii), initial charge of GAC adsorbent cost by exponential as a
    function of required mass of GAC adsorbent; and (iii), other process costs (vessels, pipes, instrumentation, and
    controls) calculated by power law as a function of total contactor(s) volume. Operating costs calculated as the
    required makeup and regeneration of GAC adsorbent. Energy for backwash and booster pumps considered negligible
    compared to regeneration costs
    """
    make_capital_cost_var(blk)
    blk.contactor_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit contactor(s) capital cost",
    )
    blk.bed_mass_gac_ref = pyo.Var(
        initialize=4,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kg,
        doc="Reference value of GAC mass needed for initial charge where "
        "economy of scale no longer discounts the unit price",
    )
    blk.adsorbent_unit_cost = pyo.Var(
        initialize=2,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency * pyo.units.kg**-1,
        doc="GAC adsorbent cost per unit mass",
    )
    blk.adsorbent_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit adsorbent capital cost",
    )
    blk.other_process_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit other process capital cost",
    )

    blk.contactor_cost_constraint = pyo.Constraint(
        expr=blk.contactor_cost
        == (
            blk.costing_package.gac.num_contactors_op
            + blk.costing_package.gac.num_contactors_redundant
        )
        * pyo.units.convert(
            (
                blk.costing_package.gac.contactor_cost_coeff_3
                * (
                    blk.unit_model.bed_volume
                    / blk.costing_package.gac.num_contactors_op
                )
                ** 3
                + blk.costing_package.gac.contactor_cost_coeff_2
                * (
                    blk.unit_model.bed_volume
                    / blk.costing_package.gac.num_contactors_op
                )
                ** 2
                + blk.costing_package.gac.contactor_cost_coeff_1
                * (
                    blk.unit_model.bed_volume
                    / blk.costing_package.gac.num_contactors_op
                )
                ** 1
                + blk.costing_package.gac.contactor_cost_coeff_0
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.bed_mass_gac_ref_constraint = pyo.Constraint(
        expr=blk.bed_mass_gac_ref
        == smooth_min(
            blk.costing_package.gac.bed_mass_max_ref / pyo.units.kg,
            pyo.units.convert(blk.unit_model.bed_mass_gac, to_units=pyo.units.kg)
            / pyo.units.kg,
        )
        * pyo.units.kg
    )
    blk.adsorbent_unit_cost_constraint = pyo.Constraint(
        expr=blk.adsorbent_unit_cost
        == pyo.units.convert(
            blk.costing_package.gac.adsorbent_unit_cost_coeff
            * pyo.exp(
                blk.bed_mass_gac_ref
                * blk.costing_package.gac.adsorbent_unit_cost_exp_coeff
            ),
            to_units=blk.costing_package.base_currency * pyo.units.kg**-1,
        )
    )
    blk.adsorbent_cost_constraint = pyo.Constraint(
        expr=blk.adsorbent_cost == blk.adsorbent_unit_cost * blk.unit_model.bed_mass_gac
    )
    blk.other_process_cost_constraint = pyo.Constraint(
        expr=blk.other_process_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.other_cost_coeff
                * ((pyo.units.m**3) ** -blk.costing_package.gac.other_cost_exp)
                * (
                    (
                        blk.costing_package.gac.num_contactors_op
                        + blk.costing_package.gac.num_contactors_redundant
                    )
                    * (
                        blk.unit_model.bed_volume
                        / blk.costing_package.gac.num_contactors_op
                    )
                )
                ** blk.costing_package.gac.other_cost_exp
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.contactor_cost + blk.adsorbent_cost + blk.other_process_cost
    )

    make_fixed_operating_cost_var(blk)
    blk.gac_regen_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
    )
    blk.gac_makeup_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Cost to makeup spent GAC adsorbent with fresh adsorbent",
    )

    blk.gac_regen_cost_constraint = pyo.Constraint(
        expr=blk.gac_regen_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.regen_unit_cost
                * (
                    blk.costing_package.gac.regen_frac
                    * blk.unit_model.gac_mass_replace_rate
                )
            ),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
    blk.gac_makeup_cost_constraint = pyo.Constraint(
        expr=blk.gac_makeup_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.makeup_unit_cost
                * (
                    (1 - blk.costing_package.gac.regen_frac)
                    * blk.unit_model.gac_mass_replace_rate
                )
            ),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == blk.gac_regen_cost + blk.gac_makeup_cost
    )
