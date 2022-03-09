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
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.generic_models.costing.costing_base import (
    FlowsheetCostingBlockData, register_idaes_currency_units)

from watertap.core.zero_order_base import ZeroOrderBase
from watertap.unit_models.zero_order import ChemicalAdditionZO


global_params = ["plant_lifetime",
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
                 "TIC"]


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
        self.plant_lifetime = pyo.Var(units=self.base_period,
                                      doc="Plant lifetime")
        self.utilization_factor = pyo.Var(
            units=pyo.units.dimensionless,
            doc='Plant capacity utilization [%]')

        self.land_cost_percent_FCI = pyo.Var(units=pyo.units.dimensionless,
                                             doc="Land cost as % FCI")
        self.working_capital_percent_FCI = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Working capital as % FCI")
        self.salaries_percent_FCI = pyo.Var(
            units=1/self.base_period,
            doc="Salaries as % FCI")
        self.benefit_percent_of_salary = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Benefits as % salaries")
        self.maintenance_costs_percent_FCI = pyo.Var(
            units=1/self.base_period,
            doc="Maintenance and contingency costs as % FCI")
        self.laboratory_fees_percent_FCI = pyo.Var(
            units=1/self.base_period,
            doc="Laboratory fees as % FCI")
        self.insurance_and_taxes_percent_FCI = pyo.Var(
            units=1/self.base_period,
            doc="Insurance and taxes as % FCI")

        self.wacc = pyo.Var(units=pyo.units.dimensionless,
                            doc='Weighted Average Cost of Capital (WACC)')
        self.capital_recovery_factor = pyo.Expression(
            expr=((self.wacc *
                   (1 + self.wacc)**(self.plant_lifetime/self.base_period)) /
                  (((1 + self.wacc)**(self.plant_lifetime/self.base_period)) -
                   1) / self.base_period))

        self.TPEC = pyo.Var(units=pyo.units.dimensionless,
                            doc='Total Purchased Equipment Cost (TPEC)')
        self.TIC = pyo.Var(units=pyo.units.dimensionless,
                           doc='Total Installed Cost (TIC)')

        # Fix all Vars from database
        for v in global_params:
            try:
                value = cs_def["global_parameters"][v]["value"]
                units = cs_def["global_parameters"][v]["units"]
                getattr(self, v).fix(value*getattr(pyo.units, units))
            except KeyError:
                raise KeyError(
                    f"Invalid case study definition file - no entry found "
                    f"for {v}, or entry lacks value and units.")

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

    def initialize_build(self):
        """
        Basic initialization for flowsheet level quantities
        """
        calculate_variable_from_constraint(
            self.land_cost, self.land_cost_constraint)
        calculate_variable_from_constraint(
            self.working_capital, self.working_capital_constraint)
        calculate_variable_from_constraint(
            self.total_capital_cost, self.total_capital_cost_constraint)

        calculate_variable_from_constraint(
            self.salary_cost, self.salary_cost_constraint)
        calculate_variable_from_constraint(
            self.benefits_cost, self.benefits_cost_constraint)
        calculate_variable_from_constraint(
            self.maintenance_cost, self.maintenance_cost_constraint)
        calculate_variable_from_constraint(
            self.laboratory_cost, self.laboratory_cost_constraint)
        calculate_variable_from_constraint(
            self.insurance_and_taxes_cost,
            self.insurance_and_taxes_cost_constraint)

        calculate_variable_from_constraint(
            self.total_fixed_operating_cost,
            self.total_fixed_operating_cost_constraint)

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
    def cost_power_law_flow(blk):
        """
        General method for costing equipment based on power law form. This is
        the most common costing form for zero-order models.

        CapCost = A*(F/Fref)**B

        This method also registers electricity demand as a costed flow (if
        present in the unit operation model).
        """
        t0 = blk.flowsheet().time.first()
        ZeroOrderCostingData._general_power_law_form(blk, time=t0)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_chemical_addition(blk):
        """
        General method for costing chemical addition processes. Capital cost is
        based on the mass flow rate of chemcial added.

        This method also registers the chemical flow and electricity demand as
        costed flows.
        """
        chem_name = blk.unit_model.config.process_subtype

        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (blk.unit_model.chemical_dosage[t0] *
                          blk.unit_model.properties[t0].flow_vol /
                          blk.unit_model.ratio_in_solution)
        sizing_term = (blk.unit_model.chemical_flow_vol[t0] /
                       (pyo.units.gal/pyo.units.day))

        ZeroOrderCostingData._general_power_law_form(
            blk, time=t0, sizing_term=sizing_term)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")
        blk.config.flowsheet_costing_block.cost_flow(
            chem_flow_mass, chem_name)

    def _general_power_law_form(blk, time=None, sizing_term=None):
        """
        General method for bulding power law costing expressions.
        """
        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        # Get costing parameter sub-block for this technology
        A, B, pblock = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype)

        if sizing_term is None:
            # Get reference state for capital calculation
            basis = parameter_dict["capital_cost"]["basis"]

            # Get state block for flow bases
            try:
                sblock = blk.unit_model.properties_in[time]
            except AttributeError:
                # Pass-through case
                sblock = blk.unit_model.properties[time]

            state_ref = pblock.reference_state[
                blk.unit_model.config.process_subtype]
            if basis == "flow_vol":
                state = sblock.flow_vol
                sizing_term = state/state_ref
            elif basis == "flow_mass":
                state = sum(sblock.flow_mass_comp[j]
                            for j in sblock.component_list)
                sizing_term = state/state_ref
            else:
                raise ValueError(
                    f"{blk.name} - unrecognized basis in parameter "
                    f"declaration: {basis}.")

        expr = pyo.units.convert(
            A*pyo.units.convert(sizing_term,
                                to_units=pyo.units.dimensionless)**B,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        if parameter_dict["capital_cost"]["cost_factor"] == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif parameter_dict["capital_cost"]["cost_factor"] == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    unit_mapping = {ZeroOrderBase: cost_power_law_flow,
                    ChemicalAdditionZO: cost_chemical_addition}


def _get_tech_parameters(blk, parameter_dict, subtype):
    """
    First, need to check to see if a Block with parameters for this technology
    exists.
    Second, to handle technology subtypes all parameters need to be indexed by
    subtype. We will dynamically add subtypes to the indexing set and Vars as
    required.
    """
    # Check to see in parameter Block already exists
    try:
        # Try to get parameter Block from costing package
        pblock = getattr(blk.config.flowsheet_costing_block,
                         blk.unit_model._tech_type)
    except AttributeError:
        # Parameter Block for this technology haven't been added yet, create
        pblock = pyo.Block()

        # Add block to FlowsheetCostingBlock
        blk.config.flowsheet_costing_block.add_component(
            blk.unit_model._tech_type, pblock)

        # Add subtype Set to Block
        pblock.subtype_set = pyo.Set()

        # Add required Vars
        pblock.capital_a_parameter = pyo.Var(
            pblock.subtype_set,
            units=getattr(
                pyo.units,
                parameter_dict[
                    "capital_cost"]["capital_a_parameter"]["units"]),
            bounds=(0, None),
            doc="Pre-exponential factor for capital cost expression")
        pblock.capital_b_parameter = pyo.Var(
            pblock.subtype_set,
            units=getattr(
                pyo.units,
                parameter_dict[
                    "capital_cost"]["capital_b_parameter"]["units"]),
            doc="Exponential factor for capital cost expression")

        if "reference_state" in parameter_dict["capital_cost"]:
            pblock.reference_state = pyo.Var(
                pblock.subtype_set,
                units=getattr(
                    pyo.units,
                    parameter_dict[
                        "capital_cost"]["reference_state"]["units"]),
                doc="Reference state for capital cost expression")

    # Check to see if requried subtype is in subtype_set
    if subtype not in pblock.subtype_set:
        # Need to add subtype and set Vars
        pblock.subtype_set.add(subtype)

        # Set vars
        pblock.capital_a_parameter[subtype].fix(
            float(parameter_dict[
                "capital_cost"]["capital_a_parameter"]["value"]) *
            getattr(pyo.units,
                    parameter_dict[
                        "capital_cost"]["capital_a_parameter"]["units"]))

        pblock.capital_b_parameter[subtype].fix(
            float(parameter_dict[
                "capital_cost"]["capital_b_parameter"]["value"]) *
            getattr(pyo.units,
                    parameter_dict[
                        "capital_cost"]["capital_b_parameter"]["units"]))

        if "reference_state" in parameter_dict["capital_cost"]:
            pblock.reference_state[subtype].fix(
                float(parameter_dict[
                    "capital_cost"]["reference_state"]["value"]) *
                getattr(pyo.units,
                        parameter_dict[
                            "capital_cost"]["reference_state"]["units"]))

    # Get requried element from Vars
    A = pblock.capital_a_parameter[subtype]
    B = pblock.capital_b_parameter[subtype]

    return A, B, pblock


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
            "..", "data", "techno_economic", "default_case_study.yaml")

    try:
        with open(source_file, "r") as f:
            lines = f.read()
            f.close()
    except OSError:
        raise OSError(
            "Could not find specified case study definition file. "
            "Please check the path provided.")

    return yaml.load(lines, yaml.Loader)
