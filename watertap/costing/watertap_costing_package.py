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

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import register_idaes_currency_units

from watertap.costing.costing_base import WaterTAPCostingBlockData


@declare_process_block_class("WaterTAPCosting")
class WaterTAPCostingData(WaterTAPCostingBlockData):
    def build_global_params(self):
        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        self._build_common_global_params()

        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Build flowsheet level costing components
        # These are the global parameters
        self.factor_total_investment = pyo.Var(
            initialize=1.0,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.factor_maintenance_labor_chemical = pyo.Var(
            initialize=0.03,
            doc="Maintenance-labor-chemical factor [fraction of equipment cost/year]",
            units=pyo.units.year**-1,
        )
        self.factor_capital_annualization = pyo.Var(
            initialize=0.1,
            doc="Capital annualization factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )
        self.capital_recovery_factor.expr = self.factor_capital_annualization

        # fix the parameters
        self.fix_all_vars()

    def build_process_costs(self):
        # add total_captial_cost and total_operating_cost
        self._build_common_process_costs()

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.factor_total_investment * self.aggregate_capital_cost
        )

        self.maintenance_labor_chemical_operating_cost = pyo.Expression(
            expr=self.factor_maintenance_labor_chemical * self.aggregate_capital_cost,
            doc="Maintenance-labor-chemical operating cost",
        )

        self.total_fixed_operating_cost = pyo.Expression(
            expr=self.aggregate_fixed_operating_cost
            + self.maintenance_labor_chemical_operating_cost,
            doc="Total fixed operating costs",
        )

        self.total_variable_operating_cost = pyo.Expression(
            expr=(
                self.aggregate_variable_operating_cost
                + sum(self.aggregate_flow_costs[f] for f in self.used_flows)
                * self.utilization_factor
            )
            if self.used_flows
            else self.aggregate_variable_operating_cost,
            doc="Total variable operating cost of process per operating period",
        )

        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == (self.total_fixed_operating_cost + self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period",
        )

        self.total_annualized_cost = pyo.Expression(
            expr=(
                self.total_capital_cost * self.capital_recovery_factor
                + self.total_operating_cost
            ),
            doc="Total annualized cost of operation",
        )

    def initialize_build(self):
        calculate_variable_from_constraint(
            self.total_capital_cost, self.total_capital_cost_constraint
        )
        calculate_variable_from_constraint(
            self.total_operating_cost, self.total_operating_cost_constraint
        )
