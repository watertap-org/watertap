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
This module contains a zero-order representation of a general unit that recovers
volatile fatty acids (VFAs).
"""

import pyomo.environ as pyo
from pyomo.environ import Reference, units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("VFARecoveryZO")
class VFARecoveryZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a VFA recovery unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "vfa_recovery"

        if "nonbiodegradable_cod" not in self.config.property_package.solute_set:
            raise ValueError(
                "nonbiodegradable_cod must be included in the solute list since"
                " this unit model computes heat requirement based on it."
            )

        build_sido(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)

        self.heat_required_per_vfa_mass = Var(
            self.flowsheet().time,
            units=pyunits.kJ / pyunits.kg,
            doc="Thermal energy required per mass VFA",
        )
        self._fixed_perf_vars.append(self.heat_required_per_vfa_mass)

        self.heat_consumption = Var(
            self.flowsheet().time,
            units=pyunits.kJ / pyunits.s,
            doc="Thermal energy required",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for heat consumption",
        )
        def eq_heat_consumption(b, t):
            return b.heat_consumption[t] == pyunits.convert(
                b.properties_in[t].flow_mass_comp["nonbiodegradable_cod"]
                * b.heat_required_per_vfa_mass[t],
                to_units=pyunits.kJ / pyunits.s,
            )

        self._perf_var_dict["Heat consumption"] = self.heat_consumption

    @property
    def default_costing_method(self):
        return self.cost_vfa_recovery

    @staticmethod
    def cost_vfa_recovery(blk):
        """
        Method for costing VFA recovery unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        unit_capex = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["unit_capex"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        capex_expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * unit_capex,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.heat_consumption[t0], "heat"
        )
