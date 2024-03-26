#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of a AOP addition
"""

import pyomo.environ as pyo


class AOPAdditionMixin:
    @staticmethod
    def _get_aop_capital_cost(blk, A, B):
        """
        Generate expression for capital cost due to AOP addition.
        """
        t0 = blk.flowsheet().time.first()

        chemical_flow_mass = pyo.units.convert(
            blk.unit_model.chemical_flow_mass[t0], to_units=pyo.units.lb / pyo.units.day
        )
        expr = pyo.units.convert(
            A
            * pyo.units.convert(
                chemical_flow_mass / (pyo.units.lb / pyo.units.day),
                to_units=pyo.units.dimensionless,
            )
            ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        return expr
