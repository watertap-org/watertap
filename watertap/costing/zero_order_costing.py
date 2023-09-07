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

        # Set the base year for all costs
        self.base_currency = getattr(pyo.units, cs_def["base_currency"])
        # Set a base period for all operating costs
        self.base_period = getattr(pyo.units, cs_def["base_period"])

        # Define expected flows
        for f, v in cs_def["defined_flows"].items():
            self.defined_flows[f] = v["value"] * getattr(pyo.units, v["units"])

        # Costing factors
        self.plant_lifetime = pyo.Var(units=self.base_period, doc="Plant lifetime")
        self.utilization_factor = pyo.Var(
            units=pyo.units.dimensionless, doc="Plant capacity utilization [%]"
        )

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
        self.capital_recovery_factor = pyo.Expression(
            expr=(
                (
                    self.wacc
                    * (1 + self.wacc) ** (self.plant_lifetime / self.base_period)
                )
                / (((1 + self.wacc) ** (self.plant_lifetime / self.base_period)) - 1)
                / self.base_period
            )
        )

        self.TPEC = pyo.Var(
            units=pyo.units.dimensionless, doc="Total Purchased Equipment Cost (TPEC)"
        )
        self.TIC = pyo.Var(
            units=pyo.units.dimensionless, doc="Total Installed Cost (TIC)"
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
        # Other capital costs
        self.land_cost = pyo.Var(
            initialize=0,
            units=self.base_currency,
            doc="Land costs - based on aggregate capital costs",
        )
        self.working_capital = pyo.Var(
            initialize=0,
            units=self.base_currency,
            doc="Working capital - based on aggregate capital costs",
        )
        self.total_capital_cost = pyo.Var(
            initialize=0, units=self.base_currency, doc="Total capital cost of process"
        )

        self.land_cost_constraint = pyo.Constraint(
            expr=self.land_cost
            == self.aggregate_capital_cost * self.land_cost_percent_FCI
        )
        self.working_capital_constraint = pyo.Constraint(
            expr=self.working_capital
            == self.aggregate_capital_cost * self.working_capital_percent_FCI
        )
        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.aggregate_capital_cost + self.land_cost + self.working_capital
        )

        # Other fixed costs
        self.salary_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Salary costs - based on aggregate capital costs",
        )
        self.benefits_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Benefits costs - based on percentage of salary costs",
        )
        self.maintenance_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Maintenance costs - based on aggregate capital costs",
        )
        self.laboratory_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Laboratory costs - based on aggregate capital costs",
        )
        self.insurance_and_taxes_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Insurance and taxes costs - based on aggregate capital costs",
        )
        self.total_fixed_operating_cost = pyo.Var(
            initialize=0,
            units=self.base_currency / self.base_period,
            doc="Total fixed operating costs",
        )

        self.salary_cost_constraint = pyo.Constraint(
            expr=self.salary_cost
            == self.aggregate_capital_cost * self.salaries_percent_FCI
        )
        self.benefits_cost_constraint = pyo.Constraint(
            expr=self.benefits_cost == self.salary_cost * self.benefit_percent_of_salary
        )
        self.maintenance_cost_constraint = pyo.Constraint(
            expr=self.maintenance_cost
            == self.aggregate_capital_cost * self.maintenance_costs_percent_FCI
        )
        self.laboratory_cost_constraint = pyo.Constraint(
            expr=self.laboratory_cost
            == self.aggregate_capital_cost * self.laboratory_fees_percent_FCI
        )
        self.insurance_and_taxes_cost_constraint = pyo.Constraint(
            expr=self.insurance_and_taxes_cost
            == self.aggregate_capital_cost * self.insurance_and_taxes_percent_FCI
        )

        self.total_fixed_operating_cost_constraint = pyo.Constraint(
            expr=self.total_fixed_operating_cost
            == self.aggregate_fixed_operating_cost
            + self.salary_cost
            + self.benefits_cost
            + self.maintenance_cost
            + self.laboratory_cost
            + self.insurance_and_taxes_cost
        )

        # Other variable costs
        self.total_variable_operating_cost = pyo.Expression(
            expr=self.aggregate_variable_operating_cost
            + sum(self.aggregate_flow_costs[f] for f in self.used_flows)
            * self.utilization_factor,
            doc="Total variable operating cost of process per operating period",
        )

        self.total_operating_cost = pyo.Expression(
            expr=(self.total_fixed_operating_cost + self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period",
        )

    def initialize_build(self):
        """
        Basic initialization for flowsheet level quantities
        """
        calculate_variable_from_constraint(self.land_cost, self.land_cost_constraint)
        calculate_variable_from_constraint(
            self.working_capital, self.working_capital_constraint
        )
        calculate_variable_from_constraint(
            self.total_capital_cost, self.total_capital_cost_constraint
        )

        calculate_variable_from_constraint(
            self.salary_cost, self.salary_cost_constraint
        )
        calculate_variable_from_constraint(
            self.benefits_cost, self.benefits_cost_constraint
        )
        calculate_variable_from_constraint(
            self.maintenance_cost, self.maintenance_cost_constraint
        )
        calculate_variable_from_constraint(
            self.laboratory_cost, self.laboratory_cost_constraint
        )
        calculate_variable_from_constraint(
            self.insurance_and_taxes_cost, self.insurance_and_taxes_cost_constraint
        )

        calculate_variable_from_constraint(
            self.total_fixed_operating_cost, self.total_fixed_operating_cost_constraint
        )

    def add_LCOW(self, flow_rate):
        """
        Add Levelized Cost of Water (LCOW) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOW
        """
        self.LCOW = pyo.Expression(
            expr=(
                self.total_capital_cost * self.capital_recovery_factor
                + self.total_operating_cost
            )
            / (
                pyo.units.convert(
                    flow_rate, to_units=pyo.units.m**3 / self.base_period
                )
                * self.utilization_factor
            ),
            doc="Levelized Cost of Water",
        )

    def add_electricity_intensity(self, flow_rate):
        """
        Add calculation of overall electricity intensity to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating electricity intensity
        """
        self.electricity_intensity = pyo.Expression(
            expr=pyo.units.convert(
                self.aggregate_flow_electricity / flow_rate,
                to_units=pyo.units.kWh / pyo.units.m**3,
            ),
            doc="Overall electricity intensity",
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
