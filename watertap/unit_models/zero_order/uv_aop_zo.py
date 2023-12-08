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
"""
This module contains a zero-order representation of a UV-AOP unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.unit_models.zero_order.uv_zo import UVZOData
from watertap.unit_models.zero_order.aop_addition_zo import AOPAdditionMixin

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVAOPZO")
class UVAOPZOData(UVZOData, AOPAdditionMixin):
    """
    Zero-Order model for a UV-AOP unit operation.
    """

    CONFIG = UVZOData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "uv_aop"

        self.oxidant_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Oxidant dosage"
        )

        self.chemical_flow_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Mass flow rate of oxidant solution",
        )

        self._fixed_perf_vars.append(self.oxidant_dose)

        @self.Constraint(self.flowsheet().time, doc="Chemical mass flow constraint")
        def chemical_flow_mass_constraint(b, t):
            return b.chemical_flow_mass[t] == pyunits.convert(
                b.oxidant_dose[t] * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        self._perf_var_dict["Oxidant Dosage (mg/L)"] = self.oxidant_dose
        self._perf_var_dict["Oxidant Flow (kg/s)"] = self.chemical_flow_mass

    @property
    def default_costing_method(self):
        return self.cost_uv_aop

    @staticmethod
    def cost_uv_aop(blk):
        t0 = blk.flowsheet().time.first()

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )
        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "reactor_cost",
                "lamp_cost",
                "aop_capital_a_parameter",
                "aop_capital_b_parameter",
            ],
        )

        expr = blk.unit_model._get_uv_capital_cost(blk, A, B)
        expr += blk.unit_model._get_aop_capital_cost(blk, C, D)

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        # TODO: Check whether chemical flow cost was accounted for originally
        # and if should be in case study verification
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass[t0], "hydrogen_peroxide"
        )
