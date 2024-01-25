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
General costing package for zero-order processes.
"""
import os
import yaml

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import register_idaes_currency_units
from watertap.costing.costing_base import WaterTAPCostingBlockData

# NOTE: some of these are defined in WaterTAPCostingBlockData
global_params = [
    "plant_lifetime",
    "utilization_factor",
    "land_cost_percent_FCI",
    "working_capital_percent_FCI",
    "salaries_percent_FCI",
    "benefit_percent_of_salary",
    "maintenance_costs_percent_FCI",
    "laboratory_fees_percent_FCI",
    "insurance_and_taxes_percent_FCI",
    "wacc",
    "TPEC",
    "TIC",
]


@declare_process_block_class("ZeroOrderCosting")
class ZeroOrderCostingData(WaterTAPCostingBlockData):
    """
    General costing package for zero-order processes.
    """

    CONFIG = WaterTAPCostingBlockData.CONFIG()
    CONFIG.declare(
        "case_study_definition",
        ConfigValue(
            default=None,
            doc="Path to YAML file defining global parameters for case study. If "
            "not provided, default values from the WaterTap database are used.",
        ),
    )

    def build_global_params(self):
        """
        To minimize overhead, only create global parameters for now.
        Unit-specific parameters will be added as sub-Blocks on a case-by-case
        basis as a unit of that type is costed.
        """
        # Load case study definition from file
        cs_def = _load_case_study_definition(self)

        # Register currency and conversion rates
        if "currency_definitions" in cs_def:
            pyo.units.load_definitions_from_strings(cs_def["currency_definitions"])
        else:
            register_idaes_currency_units()

        self._build_common_global_params()

        # Set the base year for all costs
        self.base_currency = getattr(pyo.units, cs_def["base_currency"])
        # Set a base period for all operating costs
        self.base_period = getattr(pyo.units, cs_def["base_period"])

        # Define expected flows
        for f, v in cs_def["defined_flows"].items():
            value = v["value"]
            units = getattr(pyo.units, v["units"])
            if self.component(f + "_cost") is not None:
                self.component(f + "_cost").fix(value * units)
            else:
                self.defined_flows[f] = value * units

        # Costing factors
        self.plant_lifetime = pyo.Var(units=self.base_period, doc="Plant lifetime")

        self.land_cost_percent_FCI = pyo.Var(
            units=pyo.units.dimensionless, doc="Land cost as % FCI"
        )
        self.working_capital_percent_FCI = pyo.Var(
            units=pyo.units.dimensionless, doc="Working capital as % FCI"
        )
        self.salaries_percent_FCI = pyo.Var(
            units=1 / self.base_period, doc="Salaries as % FCI"
        )
        self.benefit_percent_of_salary = pyo.Var(
            units=pyo.units.dimensionless, doc="Benefits as % salaries"
        )
        self.maintenance_costs_percent_FCI = pyo.Var(
            units=1 / self.base_period, doc="Maintenance and contingency costs as % FCI"
        )
        self.laboratory_fees_percent_FCI = pyo.Var(
            units=1 / self.base_period, doc="Laboratory fees as % FCI"
        )
        self.insurance_and_taxes_percent_FCI = pyo.Var(
            units=1 / self.base_period, doc="Insurance and taxes as % FCI"
        )

        self.wacc = pyo.Var(
            units=pyo.units.dimensionless, doc="Weighted Average Cost of Capital (WACC)"
        )
        self.capital_recovery_factor.expr = (
            (self.wacc * (1 + self.wacc) ** (self.plant_lifetime / self.base_period))
            / (((1 + self.wacc) ** (self.plant_lifetime / self.base_period)) - 1)
            / self.base_period
        )

        self.factor_total_investment = pyo.Expression(
            expr=1.0 + self.working_capital_percent_FCI + self.land_cost_percent_FCI,
            doc="Total investment factor [investment cost/equipment cost]",
        )

        # Fix all Vars from database
        for v in global_params:
            try:
                value = cs_def["global_parameters"][v]["value"]
                units = cs_def["global_parameters"][v]["units"]
                getattr(self, v).fix(value * getattr(pyo.units, units))
            except KeyError:
                raise KeyError(
                    f"Invalid case study definition file - no entry found "
                    f"for {v}, or entry lacks value and units."
                )

    def build_process_costs(self):
        """
        Calculating process wide costs.
        """
        # add total_captial_cost and total_operating_cost
        self._build_common_process_costs()

        # Other capital costs
        self.land_cost = pyo.Expression(
            expr=self.land_cost_percent_FCI * self.aggregate_capital_cost,
            doc="Land costs - based on aggregate capital costs",
        )
        self.working_capital = pyo.Expression(
            expr=self.working_capital_percent_FCI * self.aggregate_capital_cost,
            doc="Working capital - based on aggregate capital costs",
        )
        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.factor_total_investment * self.aggregate_capital_cost
        )

        # Other fixed costs
        self.salary_cost = pyo.Expression(
            expr=self.salaries_percent_FCI * self.aggregate_capital_cost,
            doc="Salary costs - based on aggregate capital costs",
        )
        self.benefits_cost = pyo.Expression(
            expr=self.benefit_percent_of_salary
            * self.salaries_percent_FCI
            * self.aggregate_capital_cost,
            doc="Benefits costs - based on percentage of salary costs",
        )
        self.maintenance_cost = pyo.Expression(
            expr=self.maintenance_costs_percent_FCI * self.aggregate_capital_cost,
            doc="Maintenance costs - based on aggregate capital costs",
        )
        self.laboratory_cost = pyo.Expression(
            expr=self.laboratory_fees_percent_FCI * self.aggregate_capital_cost,
            doc="Laboratory costs - based on aggregate capital costs",
        )
        self.insurance_and_taxes_cost = pyo.Expression(
            expr=self.insurance_and_taxes_percent_FCI * self.aggregate_capital_cost,
            doc="Insurance and taxes costs - based on aggregate capital costs",
        )
        self.factor_maintenance_labor_chemical = pyo.Expression(
            expr=self.salaries_percent_FCI
            + self.benefit_percent_of_salary * self.salaries_percent_FCI
            + self.maintenance_costs_percent_FCI
            + self.laboratory_fees_percent_FCI
            + self.insurance_and_taxes_percent_FCI,
            doc="Maintenance-labor-chemical factor [fraction of equipment cost/year]",
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

        # Other variable costs
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
        """
        Basic initialization for flowsheet level quantities
        """
        calculate_variable_from_constraint(
            self.total_capital_cost, self.total_capital_cost_constraint
        )
        calculate_variable_from_constraint(
            self.total_operating_cost, self.total_operating_cost_constraint
        )


def _load_case_study_definition(self):
    """
    Load data from case study definition file into a Python dict.
    If users did not provide a definition file as a config argument, the
    default definition from the WaterTap techno-economic database is used.
    """
    source_file = self.config.case_study_definition
    if source_file is None:
        source_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "techno_economic",
            "default_case_study.yaml",
        )

    try:
        with open(source_file, "r") as f:
            lines = f.read()
            f.close()
    except OSError:
        raise OSError(
            "Could not find specified case study definition file. "
            "Please check the path provided."
        )

    return yaml.load(lines, yaml.Loader)
