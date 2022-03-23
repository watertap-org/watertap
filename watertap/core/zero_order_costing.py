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
from watertap.unit_models.zero_order import (
    BrineConcentratorZO,
    ChemicalAdditionZO,
    ChlorinationZO,
    CoagulationFlocculationZO,
    LandfillZO,
    IonExchangeZO,
    IronManganeseRemovalZO,
    OzoneZO,
    OzoneAOPZO,
    SedimentationZO,
    StorageTankZO,
    UVZO,
    UVAOPZO,
    WellFieldZO,
    )


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

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, state_ref = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "reference_state"])

        # Get state block for flow bases
        basis = parameter_dict["capital_cost"]["basis"]
        try:
            sblock = blk.unit_model.properties_in[t0]
        except AttributeError:
            # Pass-through case
            sblock = blk.unit_model.properties[t0]

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

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_brine_concentrator(blk):
        """
        General method for costing brine concentration processes. Capital cost
        is based on the volumetirc flowrate and TDS of the incoming stream and
        the water recovery.

        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()
        inlet_state = blk.unit_model.properties_in[t0]

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter",
             "capital_b_parameter",
             "capital_c_parameter",
             "capital_d_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        expr = (
            pyo.units.convert(
                A,
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                B*inlet_state.conc_mass_comp["tds"],
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                C*blk.unit_model.recovery_frac_mass_H2O[t0],
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                D*inlet_state.flow_vol,
                to_units=blk.config.flowsheet_costing_block.base_currency))

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_chemical_addition(blk):
        """
        General method for costing chemical addition processes. Capital cost is
        based on the mass flow rate of chemical added.

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

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")
        blk.config.flowsheet_costing_block.cost_flow(
            chem_flow_mass, chem_name)

    def cost_chlorination(blk):
        """
        General method for costing chlorination units. Capital cost is based on
        the both inlet flow and dosage of chlorine.

        This method also registers the chemical flow and electricity demand as
        costed flows.
        """
        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (blk.unit_model.chlorine_dose[t0] *
                          blk.unit_model.properties_in[t0].flow_vol)

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter",
             "capital_b_parameter",
             "capital_c_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        ln_Q = pyo.log(pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol /
            (pyo.units.m**3/pyo.units.hour),
            to_units=pyo.units.dimensionless))
        ln_D = pyo.log(pyo.units.convert(
            blk.unit_model.chlorine_dose[t0] / (pyo.units.mg/pyo.units.liter),
            to_units=pyo.units.dimensionless))

        expr = pyo.units.convert(
            A*ln_Q + B*ln_D + C*ln_Q*ln_D,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")
        blk.config.flowsheet_costing_block.cost_flow(
            chem_flow_mass, "chlorine")

    def cost_coag_and_floc(blk):
        """
        General method for costing coagulation/flocculation processes. Capital cost
        is based on the alum flowrate and the polymer flowrate of the incoming stream.

        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F, G, H = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_mix_a_parameter",
             "capital_mix_b_parameter",
             "capital_floc_a_parameter",
             "capital_floc_b_parameter",
             "capital_coag_inj_a_parameter",
             "capital_coag_inj_b_parameter",
             "capital_floc_inj_a_parameter",
             "capital_floc_inj_b_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        cost_rapid_mix = (A*pyo.units.convert(
            blk.unit_model.rapid_mix_basin_vol,
            to_units=pyo.units.gallons) + B)*blk.unit_model.num_rapid_mix_processes

        cost_floc = (C*pyo.units.convert(
            blk.unit_model.floc_basin_vol,
            to_units=pyo.units.Mgallons) + D)*blk.unit_model.num_floc_processes

        cost_coag_inj = (E*pyo.units.convert(
            blk.unit_model.chemical_flow_mass[t0, "alum"],
            to_units=(pyo.units.lb/pyo.units.hr)) + F)*blk.unit_model.num_coag_processes

        cost_floc_inj = (G*pyo.units.convert(
            blk.unit_model.chemical_flow_mass[t0, "polymer"],
            to_units=(pyo.units.lb/pyo.units.day)) + H)*blk.unit_model.num_floc_injection_processes

        expr = (
            pyo.units.convert(
                cost_rapid_mix,
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                cost_floc,
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                cost_coag_inj,
                to_units=blk.config.flowsheet_costing_block.base_currency) +
            pyo.units.convert(
                cost_floc_inj,
                to_units=blk.config.flowsheet_costing_block.base_currency))

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_ion_exchange(blk):
        """
        General method for costing ion exchange units. Capital cost is based on
        the both inlet flow and TDS.

        This method also registers the NaCl demand, resin replacement and
        electricity demand as costed flows.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter",
             "capital_b_parameter",
             "capital_c_parameter",
             "capital_d_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        ln_Q = pyo.log(pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol /
            (pyo.units.m**3/pyo.units.hour),
            to_units=pyo.units.dimensionless))
        T = pyo.units.convert(
            blk.unit_model.properties_in[t0].conc_mass_comp["tds"] /
            (pyo.units.mg/pyo.units.liter),
            to_units=pyo.units.dimensionless)

        expr = pyo.units.convert(
            pyo.exp(A + B*ln_Q + C*T + D*ln_Q*T)*pyo.units.USD_2017,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.NaCl_flowrate[t0], "sodium_chloride")
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.resin_demand[t0], "ion_exchange_resin")

    def cost_iron_and_manganese_removal(blk):
        """
        General method for costing iron and manganese removal processes. Capital cost
        is based on the cost of air blower, backwash and dual media filter.

        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_blower_a_parameter",
             "capital_backwash_a_parameter",
             "capital_backwash_b_parameter",
             "capital_filter_a_parameter",
             "capital_filter_b_parameter",
             "flow_exponent"])

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        cost_blower = A

        cost_backwash = B + C*pyo.units.convert(blk.unit_model.filter_surf_area,
                                                to_units= pyo.units.ft**2)

        cost_filter = D + E * pyo.units.convert(blk.unit_model.filter_surf_area,
                                                 to_units=pyo.units.ft ** 2)

        cost_total = pyo.units.convert(cost_blower + cost_backwash + cost_filter*blk.unit_model.num_filter_units,
                                       to_units=blk.config.flowsheet_costing_block.base_currency)

        Q = pyo.units.convert(blk.unit_model.properties_in[t0].flow_vol,
                              to_units=pyo.units.m**3/pyo.units.hour)

        sizing_term = (Q / blk.unit_model.flow_basis[t0])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, cost_total, F, sizing_term, factor)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_sedimentation(blk):
        """
        General method for costing sedimentaion processes. Capital cost is
        based on the surface area of the basin.
        """
        t0 = blk.flowsheet().time.first()
        sizing_term = (blk.unit_model.basin_surface_area[t0] /
                       pyo.units.foot**2)

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_storage_tank(blk):
        """
        General method for costing storage tanks. Capital cost is based on the
        volume of the tank.
        """
        t0 = blk.flowsheet().time.first()
        sizing_term = (blk.unit_model.tank_volume[t0] / pyo.units.m**3)

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor)

    def cost_ozonation(blk):
        """
        General method for costing ozone addition. Capital cost is
        based on the inlet flowrate and dosage of ozone.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["ozone_capital_a_parameter",
             "ozone_capital_b_parameter",
             "ozone_capital_c_parameter",
             "ozone_capital_d_parameter"])
        # Get costing term for ozone addition
        expr = ZeroOrderCostingData._get_ozone_capital_cost(blk, A, B, C, D)
        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]
        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_ozonation_aop(blk):
        """
        General method for costing ozonation with AOP. Capital cost is
        based on the inlet flowrate, dosage of ozone and flow rate of H2O2.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["ozone_capital_a_parameter",
             "ozone_capital_b_parameter",
             "ozone_capital_c_parameter",
             "ozone_capital_d_parameter",
             "aop_capital_a_parameter",
             "aop_capital_b_parameter"])

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        # Get costing term for ozone addition
        expr = ZeroOrderCostingData._get_ozone_capital_cost(blk, A, B, C, D)

        # Add costing term for AOP addition
        expr += ZeroOrderCostingData._get_aop_capital_cost(blk, E, F)

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass[t0], "hydrogen_peroxide")

    def cost_uv(blk):
        """
        General method for costing UV reactor units. Capital cost is based on
        the inlet flow, UV reduced equivalent dosage, and UV transmittance at
        the inlet.
        """
        t0 = blk.flowsheet().time.first()
        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["uv_capital_a_parameter",
             "uv_capital_b_parameter",
             "uv_capital_c_parameter",
             "uv_capital_d_parameter"])

        expr = ZeroOrderCostingData._get_uv_capital_cost(blk, A, B, C, D)

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_uv_aop(blk):
        t0 = blk.flowsheet().time.first()

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")
        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["uv_capital_a_parameter",
             "uv_capital_b_parameter",
             "uv_capital_c_parameter",
             "uv_capital_d_parameter",
             "aop_capital_a_parameter",
             "aop_capital_b_parameter"])

        expr = ZeroOrderCostingData._get_uv_capital_cost(blk, A, B, C, D)
        expr += ZeroOrderCostingData._get_aop_capital_cost(blk, E, F)

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]
        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

        #TODO: Check whether chemical flow cost was accounted for originally
        # and if should be in case study verification
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass[t0], "hydrogen_peroxide")


    def cost_landfill(blk):
        """
        General method for costing landfill. Capital cost is based on the total mass and
        capacity basis.
        """

        t0 = blk.flowsheet().time.first()
        sizing_term = (blk.unit_model.total_mass[t0] / blk.unit_model.capacity_basis[t0])

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")

    def cost_well_field(blk):
        """
        General method for costing well fields. Capital cost is based on well field
        cosntruction and pipe costs.
        """

        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = \
            blk.unit_model.config.database.get_unit_operation_parameters(
                blk.unit_model._tech_type,
                subtype=blk.unit_model.config.process_subtype)

        # Get costing parameter sub-block for this technology
        A, B, pipe_cost_basis = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "pipe_cost_basis"])

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        Q = pyo.units.convert(
            blk.unit_model.properties[t0].flow_vol /
            (pyo.units.m**3/pyo.units.hour),
            to_units=pyo.units.dimensionless)
        expr = pyo.units.convert(
            A*Q**B + (pipe_cost_basis *
            blk.unit_model.pipe_distance[t0] * blk.unit_model.pipe_diameter[t0]),
            to_units=blk.config.flowsheet_costing_block.base_currency)

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity")


    def _get_ozone_capital_cost(blk, A, B, C, D):
        """
        Generate expressions for capital cost of ozonation system.
        """
        t0 = blk.flowsheet().time.first()

        ln_Q = pyo.log(pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol /
                (pyo.units.m**3/pyo.units.hour),
                to_units=pyo.units.dimensionless))
        dosage = pyo.units.convert(
            blk.unit_model.ozone_consumption[t0] /
            (pyo.units.mg/pyo.units.liter),
            to_units=pyo.units.dimensionless)

        expr = pyo.units.convert(
            A + B*dosage + C*ln_Q + D*dosage*ln_Q,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        return expr

    def _get_uv_capital_cost(blk, A, B, C, D):
        """
        Generate expression for capital cost of UV reactor.
        """
        t0 = blk.flowsheet().time.first()

        Q = pyo.units.convert(
                pyo.units.convert(blk.unit_model.properties_in[t0].flow_vol,
                                  to_units=pyo.units.Mgallons/pyo.units.day)
                / (pyo.units.Mgallons/pyo.units.day),
                to_units=pyo.units.dimensionless)

        uv_dose = pyo.units.convert(
            blk.unit_model.uv_reduced_equivalent_dose[t0] /
            (pyo.units.mJ / pyo.units.cm ** 2),
            to_units=pyo.units.dimensionless)

        uvt_in = blk.unit_model.uv_transmittance_in[t0]

        expr = pyo.units.convert(
            A * Q + B * uv_dose * Q + C * (Q * uvt_in) ** 7 +
            D * uv_dose * Q * uvt_in,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        return expr

    def _get_aop_capital_cost(blk, A, B):
        """
        Generate expression for capital cost due to AOP addition.
        """
        t0 = blk.flowsheet().time.first()

        chemical_flow_mass = pyo.units.convert(
            blk.unit_model.chemical_flow_mass[t0],
            to_units=pyo.units.lb/pyo.units.day)
        expr = pyo.units.convert(
            A*pyo.units.convert(chemical_flow_mass /
                                (pyo.units.lb/pyo.units.day),
                                to_units=pyo.units.dimensionless) ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        return expr

    def _general_power_law_form(blk, A, B, sizing_term, factor=None):
        """
        General method for building power law costing expressions.
        """
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation")

        expr = pyo.units.convert(
            A * pyo.units.convert(sizing_term,
                                  to_units=pyo.units.dimensionless) ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency)

        if factor == "TPEC":
            expr *= blk.config.flowsheet_costing_block.TPEC
        elif factor == "TIC":
            expr *= blk.config.flowsheet_costing_block.TIC

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == expr)
    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    unit_mapping = {ZeroOrderBase: cost_power_law_flow,
                    BrineConcentratorZO: cost_brine_concentrator,
                    ChemicalAdditionZO: cost_chemical_addition,
                    ChlorinationZO: cost_chlorination,
                    CoagulationFlocculationZO: cost_coag_and_floc,
                    LandfillZO: cost_landfill,
                    IonExchangeZO: cost_ion_exchange,
                    IronManganeseRemovalZO: cost_iron_and_manganese_removal,
                    OzoneZO: cost_ozonation,
                    OzoneAOPZO: cost_ozonation_aop,
                    SedimentationZO: cost_sedimentation,
                    StorageTankZO: cost_storage_tank,
                    UVZO: cost_uv,
                    UVAOPZO: cost_uv_aop,
                    WellFieldZO: cost_well_field,
                   }


def _get_tech_parameters(blk, parameter_dict, subtype, param_list):
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
        for p in param_list:
            try:
                vobj = pyo.Var(
                    pblock.subtype_set,
                    units=getattr(
                        pyo.units,
                        parameter_dict[
                            "capital_cost"][p]["units"]))
                pblock.add_component(p, vobj)
            except KeyError:
                raise KeyError(
                    "Error when trying to retrieve costing parameter"
                    " from {p} database. Please check the YAML "
                    "file for this technology for errors.")

    # Check to see if required subtype is in subtype_set
    vlist = []
    if subtype not in pblock.subtype_set:
        # Need to add subtype and set Vars
        pblock.subtype_set.add(subtype)

        # Set vars
        for p in param_list:
            vobj = getattr(pblock, p)
            vobj[subtype].fix(
                float(parameter_dict[
                    "capital_cost"][p]["value"]) *
                getattr(pyo.units,
                        parameter_dict[
                            "capital_cost"][p]["units"]))
            vlist.append(vobj[subtype])
    else:
        for p in param_list:
            vobj = getattr(pblock, p)
            vlist.append(vobj[subtype])

    return tuple(x for x in vlist)


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
