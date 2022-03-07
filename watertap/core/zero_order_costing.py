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
"""
General costing package for zero-order processes.
"""
import os
import yaml

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue

from idaes.core import declare_process_block_class
from idaes.generic_models.costing.costing_base import (
    FlowsheetCostingBlockData, register_idaes_currency_units)

from watertap.core.zero_order_base import ZeroOrderBase


global_params = ["plant_lifetime",
                 "utilization_factor",
                 "land_cost_percent_FCI",
                 "working_capital_percent_FCI",
                 "salaries_percent_FCI",
                 "benefit_percent_of_salary",
                 "maintenance_costs_percent_FCI",
                 "laboratory_fees_percent_FCI",
                 "insurance_and_taxes_percent_FCI",
                 "wacc"]


@declare_process_block_class("ZeroOrderCosting")
class ZeroOrderCostingData(FlowsheetCostingBlockData):

    CONFIG = FlowsheetCostingBlockData.CONFIG()
    CONFIG.declare("case_study_definition", ConfigValue(
        default=None,
        doc="Path to YAML file defining global parameters for case study. If "
        "not provided, default values from the WaterTap database are used."))

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
            pyo.units.load_definitions_from_strings(
                cs_def["currency_definitions"])
        else:
            register_idaes_currency_units()

        # Set the base year for all costs
        self.base_currency = getattr(pyo.units, cs_def["base_currency"])
        # Set a base period for all operating costs
        self.base_period = getattr(pyo.units, cs_def["base_period"])

        # Define expected flows
        self.defined_flows = {}
        for f, v in cs_def["defined_flows"].items():
            self.defined_flows[f] = v["value"]*getattr(pyo.units, v["units"])

        # Costing factors
        self.plant_lifetime = pyo.Var(initialize=30,
                                      units=self.base_period,
                                      doc="Plant lifetime")
        self.utilization_factor = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            doc='Plant capacity utilization [%]')

        self.land_cost_percent_FCI = pyo.Var(initialize=0.07,
                                             units=pyo.units.dimensionless,
                                             doc="Land cost as % FCI")
        self.working_capital_percent_FCI = pyo.Var(
            initialize=0.008,
            units=pyo.units.dimensionless,
            doc="Working capital as % FCI")
        self.salaries_percent_FCI = pyo.Var(
            initialize=0.07,
            units=1/self.base_period,
            doc="Salaries as % FCI")
        self.benefit_percent_of_salary = pyo.Var(
            initialize=0.07,
            units=pyo.units.dimensionless,
            doc="Benefits as % salaries")
        self.maintenance_costs_percent_FCI = pyo.Var(
            initialize=0.07,
            units=1/self.base_period,
            doc="Maintenance and contingency costs as % FCI")
        self.laboratory_fees_percent_FCI = pyo.Var(
            initialize=0.07,
            units=1/self.base_period,
            doc="Laboratory fees as % FCI")
        self.insurance_and_taxes_percent_FCI = pyo.Var(
            initialize=0.07,
            units=1/self.base_period,
            doc="Insurance and taxes as % FCI")

        self.wacc = pyo.Var(initialize=0.05,
                            units=pyo.units.dimensionless,
                            doc='Weighted Average Cost of Capital (WACC)')
        self.capital_recovery_factor = pyo.Expression(
            expr=((self.wacc *
                   (1 + self.wacc)**(self.plant_lifetime/self.base_period)) /
                  (((1 + self.wacc)**(self.plant_lifetime/self.base_period)) -
                   1) / self.base_period))

        # Fix all Vars from database
        for v in global_params:
            try:
                value = cs_def["global_parameters"][v]["value"]
                units = cs_def["global_parameters"][v]["units"]
                getattr(self, v).fix(value*getattr(pyo.units, units))
            except KeyError:
                raise KeyError(
                    f"Invalid case study definition file - no entry found "
                    f"for {v}, or entyr lacks value and units.")

    def build_process_costs(self):
        """
        Calculating process wide costs.
        """
        # Other capital costs
        self.land_cost = pyo.Var(
            initialize=0,
            units=self.base_currency,
            doc="Land costs - based on aggregate captial costs")
        self.working_capital = pyo.Var(
            initialize=0,
            units=self.base_currency,
            doc="Working capital - based on aggregate captial costs")
        self.total_capital_cost = pyo.Var(
            initialize=0,
            units=self.base_currency,
            doc="Total capital cost of process")

        self.land_cost_constraint = pyo.Constraint(
            expr=self.land_cost ==
            self.aggregate_capital_cost*self.land_cost_percent_FCI)
        self.working_capital_constraint = pyo.Constraint(
            expr=self.working_capital ==
            self.aggregate_capital_cost*self.working_capital_percent_FCI)
        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost ==
            self.aggregate_capital_cost+self.land_cost+self.working_capital)

        # Other fixed costs
        self.salary_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Salary costs - based on aggregate captial costs")
        self.benefits_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Benefits costs - based on percentage of salary costs")
        self.maintenance_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Maintenance costs - based on aggregate captial costs")
        self.laboratory_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Laboratory costs - based on aggregate captial costs")
        self.insurance_and_taxes_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Insurance and taxes costs - based on aggregate captial costs")
        self.total_fixed_operating_cost = pyo.Var(
            initialize=0,
            units=self.base_currency/self.base_period,
            doc="Total fixed operating costs")

        self.salary_cost_constraint = pyo.Constraint(
            expr=self.salary_cost ==
            self.aggregate_capital_cost*self.salaries_percent_FCI)
        self.benefits_cost_constraint = pyo.Constraint(
            expr=self.benefits_cost ==
            self.salary_cost*self.benefit_percent_of_salary)
        self.maintenance_cost_constraint = pyo.Constraint(
            expr=self.maintenance_cost ==
            self.aggregate_capital_cost*self.maintenance_costs_percent_FCI)
        self.laboratory_cost_constraint = pyo.Constraint(
            expr=self.laboratory_cost ==
            self.aggregate_capital_cost*self.laboratory_fees_percent_FCI)
        self.insurance_and_taxes_cost_constraint = pyo.Constraint(
            expr=self.insurance_and_taxes_cost ==
            self.aggregate_capital_cost*self.insurance_and_taxes_percent_FCI)

        self.total_fixed_operating_cost_constraint = pyo.Constraint(
            expr=self.total_fixed_operating_cost ==
            self.aggregate_fixed_operating_cost +
            self.salary_cost +
            self.benefits_cost +
            self.maintenance_cost +
            self.laboratory_cost +
            self.insurance_and_taxes_cost)

        # Other variable costs
        self.total_variable_operating_cost = pyo.Expression(
            expr=self.aggregate_variable_operating_cost +
            sum(self.aggregate_flow_costs[f] for f in self.flow_types),
            doc="Total variable operating cost of process per operating period")

        self.total_operating_cost = pyo.Expression(
            expr=(self.total_fixed_operating_cost +
                  self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period")

    def initialize(self):
        """
        Not needed for now, but can add custom initialization here.
        """
        pass

    def add_LCOW(self, flow_rate):
        """
        Add Levelized Cost of Water (LCOW) to costing block.

        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOW
        """
        self.LCOW = pyo.Expression(
            expr=(self.total_capital_cost*self.capital_recovery_factor +
                  self.total_operating_cost) /
                 (pyo.units.convert(
                     flow_rate,
                     to_units=pyo.units.m**3/self.base_period) *
                  self.utilization_factor),
            doc='Levelized Cost of Water')

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
                to_units=pyo.units.kWh/pyo.units.m**3),
            doc='Overall electricity intensity')

    # -------------------------------------------------------------------------
    # Unit operation costing methods
    def exponential_flow_form(blk):
        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        blk.capital_cost = pyo.Var(
            initialize=1,
            units=getattr(pyo.units, parameter_dict["capital_cost"]["units"]),
            bounds=(0, None),
            doc="Capital cost of unit operation")

        t0 = blk.flowsheet().time.first()

        # Get reference state for capital calculation
        basis = parameter_dict["capital_cost"]["basis"]

        # Get costing parameter sub-block for this technology
        try:
            # Try to get parameter Block from costing package
            pblock = getattr(blk.config.flowsheet_costing_block,
                             blk.unit_model._tech_type)
        except AttributeError:
            # If not present, call method to create parameter Block
            pblock = _add_tech_parameter_block(blk, parameter_dict)

        A = pblock.capital_a_parameter
        B = pblock.capital_b_parameter

        # Get state block for flow bases
        try:
            sblock = blk.unit_model.properties_in[t0]
        except AttributeError:
            # Pass-through case
            sblock = blk.unit_model.properties[t0]

        # TODO: More bases
        if basis == "flow_vol":
            state = sblock.flow_vol
            state_ref = pblock.reference_state
            sizing_term = state/state_ref
        elif basis == "flow_mass":
            state = sum(sblock.flow_mass_comp[j]
                        for j in sblock.component_list)
            state_ref = pblock.reference_state
            sizing_term = state/state_ref
        else:
            raise ValueError(
                f"{blk.name} - unrecognized basis in parameter declaration: "
                f"{basis}.")

        # TODO: Include TPEC/TIC
        # TODO: Chemical addition
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost ==
            A*pyo.units.convert(sizing_term,
                                to_units=pyo.units.dimensionless)**B)

        # Register flows if present
        if hasattr(blk.unit_model, "electricity"):
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.electricity[t0], "electricity")

    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    unit_mapping = {ZeroOrderBase: exponential_flow_form}


def _add_tech_parameter_block(blk, parameter_dict):
    # Parameters for this technology haven't been added yet
    pblock = pyo.Block()

    # Add block to FlowsheetCostingBlock
    blk.config.flowsheet_costing_block.add_component(
        blk.unit_model._tech_type, pblock)

    # Add required parameters
    pblock.capital_a_parameter = pyo.Var(
        initialize=float(
            parameter_dict[
                "capital_cost"]["capital_a_parameter"]["value"]),
        units=getattr(
            pyo.units,
            parameter_dict[
                "capital_cost"]["capital_a_parameter"]["units"]),
        bounds=(0, None),
        doc="Pre-exponential factor for capital cost expression")
    pblock.capital_b_parameter = pyo.Var(
        initialize=float(
            parameter_dict[
                "capital_cost"]["capital_b_parameter"]["value"]),
        units=getattr(
            pyo.units,
            parameter_dict[
                "capital_cost"]["capital_b_parameter"]["units"]),
        doc="Exponential factor for capital cost expression")

    pblock.capital_a_parameter.fix()
    pblock.capital_b_parameter.fix()

    if parameter_dict["capital_cost"]["basis"] in ["flow_vol", "flow_mass"]:
        # Flow based costing requires a reference flow
        pblock.reference_state = pyo.Var(
            initialize=float(
                parameter_dict[
                    "capital_cost"]["reference_state"]["value"]),
            units=getattr(
                pyo.units,
                parameter_dict[
                    "capital_cost"]["reference_state"]["units"]),
            doc="Reference state for capital cost expression")

        pblock.reference_state.fix()

    return pblock


def _load_case_study_definition(self):
    source_file = self.config.case_study_definition
    if source_file is None:
        source_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "data", "techno_economic", "default_case_study.yaml")

    try:
        with open(source_file, "r") as f:
            lines = f.read()
            f.close()
    except OSError:
        raise KeyError(
            "Could not find specified case study definition file. "
            "Please check the path provided.")

    return yaml.load(lines, yaml.Loader)
