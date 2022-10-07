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
from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

from watertap.core.zero_order_base import ZeroOrderBase
from watertap.unit_models.zero_order import (
    AnaerobicMBRMECZO,
    ATHTLZO,
    BrineConcentratorZO,
    CANDOPZO,
    CentrifugeZO,
    ChemicalAdditionZO,
    ChlorinationZO,
    ClarifierZO,
    ClothMediaFiltrationZO,
    CoagulationFlocculationZO,
    CofermentationZO,
    ConstructedWetlandsZO,
    DeepWellInjectionZO,
    DMBRZO,
    ElectroNPZO,
    EvaporationPondZO,
    FilterPressZO,
    FixedBedZO,
    GACZO,
    HRCSZO,
    HTGZO,
    IonExchangeZO,
    IronManganeseRemovalZO,
    LandfillZO,
    MABRZO,
    MagprexZO,
    MembraneEvaporatorZO,
    MetabZO,
    MicrobialBatteryZO,
    NanofiltrationZO,
    OzoneZO,
    OzoneAOPZO,
    PeraceticAcidDisinfectionZO,
    PhotothermalMembraneZO,
    PumpElectricityZO,
    SaltPrecipitationZO,
    SedimentationZO,
    StorageTankZO,
    StruviteClassifierZO,
    SuboxicASMZO,
    SurfaceDischargeZO,
    UVZO,
    UVAOPZO,
    VFARecoveryZO,
    WellFieldZO,
)

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
class ZeroOrderCostingData(FlowsheetCostingBlockData):
    """
    General costing package for zero-order processes.
    """

    CONFIG = FlowsheetCostingBlockData.CONFIG()
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
        self.defined_flows = {}
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
            + sum(self.aggregate_flow_costs[f] for f in self.flow_types)
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

    # -------------------------------------------------------------------------
    # Unit operation costing methods
    def cost_power_law_flow(blk, number_of_parallel_units=1):
        """
        General method for costing equipment based on power law form. This is
        the most common costing form for zero-order models.
        CapCost = A*(F/Fref)**B
        This method also registers electricity demand as a costed flow (if
        present in the unit operation model).
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, state_ref = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "reference_state"],
        )

        # Get state block for flow bases
        basis = parameter_dict["capital_cost"]["basis"]
        try:
            sblock = blk.unit_model.properties_in[t0]
        except AttributeError:
            # Pass-through case
            sblock = blk.unit_model.properties[t0]

        if basis == "flow_vol":
            state = sblock.flow_vol
            sizing_term = state / state_ref
        elif basis == "flow_mass":
            state = sum(sblock.flow_mass_comp[j] for j in sblock.component_list)
            sizing_term = state / state_ref
        else:
            raise ValueError(
                f"{blk.name} - unrecognized basis in parameter "
                f"declaration: {basis}."
            )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def _add_cost_factor(blk, factor):
        if factor == "TPEC":
            blk.cost_factor = pyo.Expression(
                expr=blk.config.flowsheet_costing_block.TPEC
            )
        elif factor == "TIC":
            blk.cost_factor = pyo.Expression(
                expr=blk.config.flowsheet_costing_block.TIC
            )
        else:
            blk.cost_factor = pyo.Expression(expr=1.0)
        blk.direct_capital_cost = pyo.Expression(
            expr=blk.capital_cost / blk.cost_factor
        )

    def cost_anaerobic_mbr_mec(blk):
        """
        Method for costing anaerobic membrane bioreactor integrated with
        microbial electrolysis cell.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        unit_capex, unit_opex = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["unit_capex", "unit_opex"],
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
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        # Add fixed operating cost variable and constraint
        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit",
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_autothermal_hydrothermal_liquefaction(blk):
        """
        General method for costing autothermal-hydrothermal liquefaction unit. Capital cost
        is based on the HTL reactor, booster pump, solid filter, other equipment, and
        heat oil system.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            A,
            B,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
            J,
            K,
            L,
            M,
            N,
            O,
            P,
            Q,
            R,
            S,
            T,
        ) = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "installation_factor_reactor",
                "equipment_cost_reactor",
                "base_flowrate_reactor",
                "scaling_exponent_reactor",
                "installation_factor_pump",
                "equipment_cost_pump",
                "base_flowrate_pump",
                "scaling_exponent_pump",
                "installation_factor_other",
                "equipment_cost_other",
                "base_flowrate_other",
                "scaling_exponent_other",
                "installation_factor_solid_filter",
                "equipment_cost_solid_filter",
                "base_flowrate_solid_filter",
                "scaling_exponent_solid_filter",
                "installation_factor_heat",
                "equipment_cost_heat",
                "base_flowrate_heat",
                "scaling_exponent_heat",
            ],
        )

        sizing_term_reactor = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / C),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_pump = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / G),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_other = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / K),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_solid_filter = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / O),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_heat = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / S),
            to_units=pyo.units.dimensionless,
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        reactor_cost = pyo.units.convert(
            A * B * sizing_term_reactor**D,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        pump_cost = pyo.units.convert(
            E * F * sizing_term_pump**H,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        other_cost = pyo.units.convert(
            I * J * sizing_term_other**L,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        solid_filter_cost = pyo.units.convert(
            M * N * sizing_term_solid_filter**P,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        heat_cost = pyo.units.convert(
            Q * R * sizing_term_heat**T,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = reactor_cost + pump_cost + other_cost + solid_filter_cost + heat_cost

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.catalyst_flow[t0], "catalyst_ATHTL"
        )

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
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_a_parameter",
                "capital_b_parameter",
                "capital_c_parameter",
                "capital_d_parameter",
            ],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = (
            pyo.units.convert(
                A, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                B * inlet_state.conc_mass_comp["tds"],
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                C * blk.unit_model.recovery_frac_mass_H2O[t0],
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                D * inlet_state.flow_vol,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_chemical_addition(blk, number_of_parallel_units=1):
        """
        General method for costing chemical addition processes. Capital cost is
        based on the mass flow rate of chemical added.
        This method also registers the chemical flow and electricity demand as
        costed flows.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        chem_name = blk.unit_model.config.process_subtype

        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (
            blk.unit_model.chemical_dosage[t0]
            * blk.unit_model.properties[t0].flow_vol
            / blk.unit_model.ratio_in_solution
        )
        sizing_term = blk.unit_model.chemical_flow_vol[t0] / (
            pyo.units.gal / pyo.units.day
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(chem_flow_mass, chem_name)

    def cost_chlorination(blk):
        """
        General method for costing chlorination units. Capital cost is based on
        the both inlet flow and dosage of chlorine.
        This method also registers the chemical flow and electricity demand as
        costed flows.
        """
        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (
            blk.unit_model.chlorine_dose[t0] * blk.unit_model.properties_in[t0].flow_vol
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "capital_c_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        ln_Q = pyo.log(
            pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol
                / (pyo.units.m**3 / pyo.units.hour),
                to_units=pyo.units.dimensionless,
            )
        )
        ln_D = pyo.log(
            pyo.units.convert(
                blk.unit_model.chlorine_dose[t0] / (pyo.units.mg / pyo.units.liter),
                to_units=pyo.units.dimensionless,
            )
        )

        expr = pyo.units.convert(
            A * ln_Q + B * ln_D + C * ln_Q * ln_D,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(chem_flow_mass, "chlorine")

    def cost_coag_and_floc(blk):
        """
        General method for costing coagulation/flocculation processes. Capital cost
        is based on the alum flowrate and the polymer flowrate of the incoming stream.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F, G, H = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_mix_a_parameter",
                "capital_mix_b_parameter",
                "capital_floc_a_parameter",
                "capital_floc_b_parameter",
                "capital_coag_inj_a_parameter",
                "capital_coag_inj_b_parameter",
                "capital_floc_inj_a_parameter",
                "capital_floc_inj_b_parameter",
            ],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        cost_rapid_mix = (
            A
            * pyo.units.convert(
                blk.unit_model.rapid_mix_basin_vol, to_units=pyo.units.gallons
            )
            + B
        ) * blk.unit_model.num_rapid_mix_processes

        cost_floc = (
            C
            * pyo.units.convert(
                blk.unit_model.floc_basin_vol, to_units=pyo.units.Mgallons
            )
            + D
        ) * blk.unit_model.num_floc_processes

        cost_coag_inj = (
            E
            * pyo.units.convert(
                blk.unit_model.chemical_flow_mass[t0, "alum"],
                to_units=(pyo.units.lb / pyo.units.hr),
            )
            + F
        ) * blk.unit_model.num_coag_processes

        cost_floc_inj = (
            G
            * pyo.units.convert(
                blk.unit_model.chemical_flow_mass[t0, "polymer"],
                to_units=(pyo.units.lb / pyo.units.day),
            )
            + H
        ) * blk.unit_model.num_floc_injection_processes

        expr = (
            pyo.units.convert(
                cost_rapid_mix,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                cost_floc, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                cost_coag_inj, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                cost_floc_inj, to_units=blk.config.flowsheet_costing_block.base_currency
            )
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_cofermentation(blk):
        """
        Method for costing cofermentation unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        unit_capex, unit_opex = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["unit_capex", "unit_opex"],
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
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        # Add fixed operating cost variable and constraint
        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit",
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_constructed_wetlands(blk):
        """
        Method for costing constructed wetlands.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )
        # Get costing parameter sub-block for this technology
        unit_capex = _get_tech_parameters(
            blk, parameter_dict, blk.unit_model.config.process_subtype, ["unit_capex"]
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
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

    def cost_deep_well_injection(blk, number_of_parallel_units=1):
        """
        General method for costing deep well injection processes. Capital cost
        is based on the cost of pump and pipe.
        This method also registers the electricity demand as a costed flow.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["well_pump_cost", "pipe_cost_basis", "flow_exponent"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        cost_well_pump = A

        cost_pipe = (
            B * blk.unit_model.pipe_distance[t0] * blk.unit_model.pipe_diameter[t0]
        )

        cost_total = pyo.units.convert(
            cost_well_pump + cost_pipe,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        Q = pyo.units.convert(
            blk.unit_model.properties[t0].flow_vol,
            to_units=pyo.units.m**3 / pyo.units.hour,
        )

        sizing_term = Q / blk.unit_model.flow_basis[t0]

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk,
            cost_total,
            C,
            sizing_term,
            factor,
            number_of_parallel_units,
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_dmbr(blk):
        """
        General method for costing dynamic membrane bioreactor. Capital cost
        is based on the cost of membrane area.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["water_flux", "reactor_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_treated[t0].flow_vol / A * B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_electrochemical_nutrient_removal(blk):
        """
        General method for costing electrochemical nutrient recovery. Capital cost
        is based on the volumetirc flowrate and HRT of the incoming stream. Chemical
        dosing cost is based on MgCl2 cost.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()
        inlet_state = blk.unit_model.properties_in[t0]

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "HRT",
                "sizing_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A * inlet_state.flow_vol * B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.MgCl2_flowrate[t0], "magnesium_chloride"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp["struvite"],
            "struvite_product",
        )

    def cost_fixed_bed(blk, number_of_parallel_units=1):
        """
        General method for costing fixed bed units. This primarily calls the
        cost_power_law_flow method.
        This method also registers demand for a number of additional material
        flows.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        ZeroOrderCostingData.cost_power_law_flow(blk, number_of_parallel_units)

        # Register flows - electricity already done by cost_power_law_flow
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.acetic_acid_demand[t0], "acetic_acid"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.phosphoric_acid_demand[t0], "phosphoric_acid"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.ferric_chloride_demand[t0], "ferric_chloride"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.activated_carbon_demand[t0], "activated_carbon"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.sand_demand[t0], "sand"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.anthracite_demand[t0], "anthracite"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.cationic_polymer_demand[t0], "cationic_polymer"
        )

    def cost_gac(blk):
        """
        General method for costing granular activated carbon processes. Capital
        cost is based on the inlet flow rate of liquid and the empty bed
        contacting time.
        This method also registers electricity and activated carbon consumption
        as costed flows.
        """
        t0 = blk.flowsheet().time.first()

        Q = blk.unit_model.properties_in[t0].flow_vol
        T = blk.unit_model.empty_bed_contact_time

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        A, B, C = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "capital_c_parameter"],
        )

        # Call general power law costing method
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = (
            pyo.units.convert(
                A * Q, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                B * T, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                C * Q * T, to_units=blk.config.flowsheet_costing_block.base_currency
            )
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.activated_carbon_demand[t0], "activated_carbon"
        )

    def cost_hydrothermal_gasification(blk):
        """
        General method for costing hydrothermal gasification unit. Capital cost
        is based on the CHG reactor and other wastewater treatment equipment including
        a feed pump, a booster pump, a feed/product exchanger, a fired heater,
        a hydrocyclone, and a product air fin cooler.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            IF_reactor,
            EP_reactor,
            F0_reactor,
            SE_reactor,
            IF_pump,
            EP_pump,
            F0_pump,
            SE_pump,
            IF_booster,
            EP_booster,
            F0_booster,
            SE_booster,
            IF_hydrocyclone,
            EP_hydrocyclone,
            F0_hydrocyclone,
            SE_hydrocyclone,
            IF_cooler,
            EP_cooler,
            F0_cooler,
            SE_cooler,
            IF_exchanger,
            EP_exchanger,
            F0_exchanger,
            Fnew_exchanger,
            SE_exchanger,
            IF_heater,
            EP_heater,
            F0_heater,
            Fnew_heater,
            SE_heater,
        ) = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "installation_factor_reactor",
                "equipment_cost_reactor",
                "base_flowrate_reactor",
                "scaling_exponent_reactor",
                "installation_factor_pump",
                "equipment_cost_pump",
                "base_flowrate_pump",
                "scaling_exponent_pump",
                "installation_factor_booster",
                "equipment_cost_booster",
                "base_flowrate_booster",
                "scaling_exponent_booster",
                "installation_factor_hydrocyclone",
                "equipment_cost_hydrocyclone",
                "base_flowrate_hydrocyclone",
                "scaling_exponent_hydrocyclone",
                "installation_factor_cooler",
                "equipment_cost_cooler",
                "base_flowrate_cooler",
                "scaling_exponent_cooler",
                "installation_factor_exchanger",
                "equipment_cost_exchanger",
                "base_area_exchanger",
                "new_area_exchanger",
                "scaling_exponent_exchanger",
                "installation_factor_heater",
                "equipment_cost_heater",
                "base_heat_duty_heater",
                "new_heat_duty_heater",
                "scaling_exponent_heater",
            ],
        )

        sizing_term_reactor = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_reactor),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_pump = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_pump),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_booster = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_booster),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_hydrocyclone = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_hydrocyclone),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_cooler = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_cooler),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_exchanger = pyo.units.convert(
            (Fnew_exchanger / F0_exchanger),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_heater = pyo.units.convert(
            (Fnew_heater / F0_heater),
            to_units=pyo.units.dimensionless,
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        reactor_cost = pyo.units.convert(
            IF_reactor * EP_reactor * sizing_term_reactor**SE_reactor,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        pump_cost = pyo.units.convert(
            IF_pump * EP_pump * sizing_term_pump**SE_pump,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        booster_cost = pyo.units.convert(
            IF_booster * EP_booster * sizing_term_booster**SE_booster,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        hydrocyclone_cost = pyo.units.convert(
            IF_hydrocyclone
            * EP_hydrocyclone
            * sizing_term_hydrocyclone**SE_hydrocyclone,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        cooler_cost = pyo.units.convert(
            IF_cooler * EP_cooler * sizing_term_cooler**SE_cooler,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        exchanger_cost = pyo.units.convert(
            IF_exchanger * EP_exchanger * sizing_term_exchanger**SE_exchanger,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        heater_cost = pyo.units.convert(
            IF_heater * EP_heater * sizing_term_heater**SE_heater,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = (
            reactor_cost
            + pump_cost
            + booster_cost
            + hydrocyclone_cost
            + cooler_cost
            + exchanger_cost
            + heater_cost
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.catalyst_flow[t0], "catalyst_HTG"
        )

    def cost_mabr(blk):
        """
        General method for costing membrane aerated biofilm reactor. Capital cost
        is based on the cost of reactor and blower.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["reactor_cost", "blower_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        DCC_reactor = pyo.units.convert(
            blk.unit_model.properties_treated[t0].flow_mass_comp["ammonium_as_nitrogen"]
            / blk.unit_model.nitrogen_removal_rate
            * A,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        DCC_blower = pyo.units.convert(
            blk.unit_model.reactor_area * blk.unit_model.air_flow_rate[t0] * B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = DCC_reactor + DCC_blower

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_metab(blk):
        """
        General method for costing the metab reactor. Capital cost
        is based on the cost of reactor, mixer, METAB beads, membrane,
        and vacuum pump.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            reactor_cost,
            mixer_cost,
            bead_bulk_density,
            bead_cost,
            bead_replacement_factor,
            membrane_sidestream_fraction,
            membrane_specific_size,
            membrane_cost,
            vacuum_cost,
        ) = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "reactor_cost",
                "mixer_cost",
                "bead_bulk_density",
                "bead_cost",
                "bead_replacement_factor",
                "membrane_sidestream_fraction",
                "membrane_specific_size",
                "membrane_cost",
                "vacuum_cost",
            ],
        )

        # Add capital cost variables and constraints
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )
        blk.DCC_reactor = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of reactor",
        )
        blk.DCC_mixer = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of mixer",
        )
        blk.DCC_bead = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of beads",
        )
        blk.DCC_membrane = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of membrane",
        )
        blk.DCC_vacuum = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of vacuum pump",
        )
        blk.eq_DCC_reactor = pyo.Constraint(
            expr=blk.DCC_reactor
            == pyo.units.convert(
                blk.unit_model.volume * reactor_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_mixer = pyo.Constraint(
            expr=blk.DCC_mixer
            == pyo.units.convert(
                blk.unit_model.energy_electric_mixer_vol
                * blk.unit_model.volume
                * mixer_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_bead = pyo.Constraint(
            expr=blk.DCC_bead
            == pyo.units.convert(
                blk.unit_model.volume * bead_bulk_density * bead_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_membrane = pyo.Constraint(
            expr=blk.DCC_membrane
            == pyo.units.convert(
                blk.unit_model.get_inlet_flow(t0)
                * membrane_sidestream_fraction
                * membrane_specific_size
                * membrane_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_vacuum = pyo.Constraint(
            expr=blk.DCC_vacuum
            == pyo.units.convert(
                blk.unit_model.properties_byproduct[t0].flow_mass_comp[
                    blk.unit_model._gas_comp
                ]
                * vacuum_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )

        expr = (
            blk.DCC_reactor
            + blk.DCC_mixer
            + blk.DCC_bead
            + blk.DCC_membrane
            + blk.DCC_vacuum
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Add fixed operating cost variable and constraint
        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit",
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.DCC_bead * bead_replacement_factor,
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

        # Register operating cost flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(blk.unit_model.heat[t0], "heat")
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp[
                blk.unit_model._gas_comp
            ],
            blk.unit_model._gas_comp + "_product",
        )

    def cost_ion_exchange(blk):
        """
        Two methods for costing ion exchange:
        (1) General method for costing ion exchange units. Capital cost is based on
        the both inlet flow and TDS.
        This method also registers the NaCl demand, resin replacement and
        electricity demand as costed flows.
        (2) General method using unit capex and unit opex cost parameters, tailored for AMO
        wastewater resource recovery (process subtype: clinoptilolite)
        """
        t0 = blk.flowsheet().time.first()

        if blk.unit_model.config.process_subtype != "clinoptilolite":
            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )
            # Get costing parameter sub-block for this technology
            A, B, C, D = _get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                [
                    "capital_a_parameter",
                    "capital_b_parameter",
                    "capital_c_parameter",
                    "capital_d_parameter",
                ],
            )

            # Add cost variable and constraint
            blk.capital_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency,
                bounds=(0, None),
                doc="Capital cost of unit operation",
            )

            ln_Q = pyo.log(
                pyo.units.convert(
                    blk.unit_model.properties_in[t0].flow_vol
                    / (pyo.units.m**3 / pyo.units.hour),
                    to_units=pyo.units.dimensionless,
                )
            )
            T = pyo.units.convert(
                blk.unit_model.properties_in[t0].conc_mass_comp["tds"]
                / (pyo.units.mg / pyo.units.liter),
                to_units=pyo.units.dimensionless,
            )

            expr = pyo.units.convert(
                pyo.exp(A + B * ln_Q + C * T + D * ln_Q * T) * pyo.units.USD_2017,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )

            ZeroOrderCostingData._add_cost_factor(
                blk, parameter_dict["capital_cost"]["cost_factor"]
            )

            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost == blk.cost_factor * expr
            )

            # Register flows
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.electricity[t0], "electricity"
            )
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.NaCl_flowrate[t0], "sodium_chloride"
            )
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.resin_demand[t0], "ion_exchange_resin"
            )
        else:
            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )
            # Get costing parameter sub-block for this technology
            unit_capex, unit_opex = _get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                ["unit_capex", "unit_opex"],
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
            ZeroOrderCostingData._add_cost_factor(
                blk, parameter_dict["capital_cost"]["cost_factor"]
            )

            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost == blk.cost_factor * capex_expr
            )

            # Add fixed operating cost variable and constraint
            blk.fixed_operating_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
                bounds=(0, None),
                doc="Fixed operating cost of unit",
            )
            blk.fixed_operating_cost_constraint = pyo.Constraint(
                expr=blk.fixed_operating_cost
                == pyo.units.convert(
                    blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                    to_units=blk.config.flowsheet_costing_block.base_currency
                    / blk.config.flowsheet_costing_block.base_period,
                )
            )

            # Register flows
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.electricity[t0], "electricity"
            )

    def cost_iron_and_manganese_removal(blk, number_of_parallel_units=1):
        """
        General method for costing iron and manganese removal processes. Capital cost
        is based on the cost of air blower, backwash and dual media filter.
        This method also registers the electricity demand as a costed flow.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_blower_a_parameter",
                "capital_backwash_a_parameter",
                "capital_backwash_b_parameter",
                "capital_filter_a_parameter",
                "capital_filter_b_parameter",
                "flow_exponent",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        cost_blower = A

        cost_backwash = B + C * pyo.units.convert(
            blk.unit_model.filter_surf_area, to_units=pyo.units.ft**2
        )

        cost_filter = D + E * pyo.units.convert(
            blk.unit_model.filter_surf_area, to_units=pyo.units.ft**2
        )

        cost_total = pyo.units.convert(
            cost_blower + cost_backwash + cost_filter * blk.unit_model.num_filter_units,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.m**3 / pyo.units.hour,
        )

        sizing_term = Q / blk.unit_model.flow_basis[t0]

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk,
            cost_total,
            F,
            sizing_term,
            factor,
            number_of_parallel_units,
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_pump_electricity(blk):
        """
        General method for costing low pressure pump. Capital cost
        is based on the cost of inlet flow rate.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A = _get_tech_parameters(
            blk, parameter_dict, blk.unit_model.config.process_subtype, ["pump_cost"]
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties[t0].flow_vol * A,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_sedimentation(blk, number_of_parallel_units=1):
        """
        General method for costing sedimentaion processes. Capital cost is
        based on the surface area of the basin.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        if blk.unit_model.config.process_subtype != "phosphorus_capture":
            sizing_term = blk.unit_model.basin_surface_area[t0] / pyo.units.foot**2

            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )

            A, B = _get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                ["capital_a_parameter", "capital_b_parameter"],
            )

            # Determine if a costing factor is required
            factor = parameter_dict["capital_cost"]["cost_factor"]

            # Call general power law costing method
            ZeroOrderCostingData._general_power_law_form(
                blk, A, B, sizing_term, factor, number_of_parallel_units
            )
        else:
            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )

            # Get costing parameter sub-block for this technology
            unit_capex, unit_opex = _get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                ["unit_capex", "unit_opex"],
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
            ZeroOrderCostingData._add_cost_factor(
                blk, parameter_dict["capital_cost"]["cost_factor"]
            )

            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost == blk.cost_factor * capex_expr
            )

            # Add fixed operating cost variable and constraint
            blk.fixed_operating_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
                bounds=(0, None),
                doc="Fixed operating cost of unit",
            )
            blk.fixed_operating_cost_constraint = pyo.Constraint(
                expr=blk.fixed_operating_cost
                == pyo.units.convert(
                    blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                    to_units=blk.config.flowsheet_costing_block.base_currency
                    / blk.config.flowsheet_costing_block.base_period,
                )
            )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_storage_tank(blk, number_of_parallel_units=1):
        """
        General method for costing storage tanks. Capital cost is based on the
        volume of the tank.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()
        sizing_term = blk.unit_model.tank_volume[t0] / pyo.units.m**3

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_suboxic_asm(blk):
        """
        General method for costing suboxic activated sludge process unit. Capital cost
        is based on the aeration basin, other equipments including mixers, blowers, MLR pumps,
        RAS pumps and automated valves, and instrumentation and control system including
        probes (dissolved oxygen, nitrate and ammonia), phosphorus analyzer and air flowmeter.
        """
        t0 = blk.flowsheet().time.first()
        flow_in = blk.unit_model.properties_in[t0].flow_vol

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "aeration_basin_cost",
                "other_equipment_cost",
                "control_system_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        aeration_basin_cost = pyo.units.convert(
            A * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        other_equipment_cost = pyo.units.convert(
            B * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        control_system_cost = pyo.units.convert(
            C * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = aeration_basin_cost + other_equipment_cost + control_system_cost

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_surface_discharge(blk):
        """
        General method for costing surface discharge. Capital cost is based on
        construction and pipe costs.
        """

        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, pipe_cost_basis, ref_state = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_a_parameter",
                "capital_b_parameter",
                "pipe_cost_basis",
                "reference_state",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A
            * pyo.units.convert(
                blk.unit_model.properties[t0].flow_vol / ref_state,
                to_units=pyo.units.dimensionless,
            )
            ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        ) + pyo.units.convert(
            pipe_cost_basis
            * blk.unit_model.pipe_distance[t0]
            * blk.unit_model.pipe_diameter[t0],
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_ozonation(blk):
        """
        General method for costing ozone addition. Capital cost is
        based on the inlet flowrate and dosage of ozone.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "ozone_capital_a_parameter",
                "ozone_capital_b_parameter",
                "ozone_capital_c_parameter",
                "ozone_capital_d_parameter",
            ],
        )
        # Get costing term for ozone addition
        expr = ZeroOrderCostingData._get_ozone_capital_cost(blk, A, B, C, D)

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_ozonation_aop(blk):
        """
        General method for costing ozonation with AOP. Capital cost is
        based on the inlet flowrate, dosage of ozone and flow rate of H2O2.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "ozone_capital_a_parameter",
                "ozone_capital_b_parameter",
                "ozone_capital_c_parameter",
                "ozone_capital_d_parameter",
                "aop_capital_a_parameter",
                "aop_capital_b_parameter",
            ],
        )

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        # Get costing term for ozone addition
        expr = ZeroOrderCostingData._get_ozone_capital_cost(blk, A, B, C, D)

        # Add costing term for AOP addition
        expr += ZeroOrderCostingData._get_aop_capital_cost(blk, E, F)

        blk.capital_cost_constraint = pyo.Constraint(expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass[t0], "hydrogen_peroxide"
        )

    def cost_supercritical_salt_precipitation(blk):
        """
        General method for costing supercritical salt precipitation unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "installation_factor",
                "equipment_cost",
                "base_flowrate",
                "scaling_exponent",
            ],
        )

        sizing_term = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / C),
            to_units=pyo.units.dimensionless,
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A * B * sizing_term**D,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

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
            doc="Capital cost of unit operation",
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "reactor_cost",
                "lamp_cost",
            ],
        )

        expr = ZeroOrderCostingData._get_uv_capital_cost(blk, A, B)

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

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
        A, B, C, D = _get_tech_parameters(
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

        expr = ZeroOrderCostingData._get_uv_capital_cost(blk, A, B)
        expr += ZeroOrderCostingData._get_aop_capital_cost(blk, C, D)

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
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

    def cost_evaporation_pond(blk):
        """
        General method for costing evaporation pond. Capital cost is based on the pond area and
        other pond construction parameters.
        """

        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            A,
            B,
            C,
            D,
            E,
            liner_thickness,
            land_cost,
            land_clearing_cost,
        ) = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "cost_per_acre_a_parameter",
                "cost_per_acre_b_parameter",
                "cost_per_acre_c_parameter",
                "cost_per_acre_d_parameter",
                "cost_per_acre_e_parameter",
                "liner_thickness",
                "land_cost",
                "land_clearing_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.adj_area[t0]
            * (
                A
                + B * liner_thickness
                + C * land_cost
                + D * land_clearing_cost
                + E * blk.unit_model.dike_height[t0]
            ),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        blk.capital_cost_constraint = pyo.Constraint(expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_filter_press(blk):
        """
        General method for costing belt filter press. Capital cost is a function
        of flow in gal/hr.
        """
        t0 = blk.flowsheet().time.first()
        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.gal / pyo.units.hr,
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        expr = pyo.units.convert(
            A * Q + B, to_units=blk.config.flowsheet_costing_block.base_currency
        )

        blk.capital_cost_constraint = pyo.Constraint(expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_landfill(blk, number_of_parallel_units=1):
        """
        General method for costing landfill. Capital cost is based on the total mass and
        capacity basis.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """

        t0 = blk.flowsheet().time.first()
        sizing_term = blk.unit_model.total_mass[t0] / blk.unit_model.capacity_basis[t0]

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        ZeroOrderCostingData._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_well_field(blk):
        """
        General method for costing well fields. Capital cost is based on well field
        construction and pipe costs.
        """

        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, pipe_cost_basis = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "pipe_cost_basis"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        Q = pyo.units.convert(
            blk.unit_model.properties[t0].flow_vol
            / (pyo.units.m**3 / pyo.units.hour),
            to_units=pyo.units.dimensionless,
        )
        expr = pyo.units.convert(
            A * Q**B
            + (
                pipe_cost_basis
                * blk.unit_model.pipe_distance[t0]
                * blk.unit_model.pipe_diameter[t0]
            ),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_nanofiltration(blk, number_of_parallel_units=1):
        """
        General method for costing nanofiltration. Costing is carried out
        using either the general_power_law form or the standard form which
        computes membrane cost and replacement rate.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        # Get cost method for this technology
        cost_method = _get_unit_cost_method(blk)
        valid_methods = ["cost_power_law_flow", "cost_membrane"]
        if cost_method == "cost_power_law_flow":
            ZeroOrderCostingData.cost_power_law_flow(blk, number_of_parallel_units)
        elif cost_method == "cost_membrane":
            # NOTE: number of units does not matter for cost_membrane
            #       as its a linear function of membrane area
            ZeroOrderCostingData.cost_membrane(blk)
        else:
            raise KeyError(
                f"{cost_method} is not a relevant cost method for "
                f"{blk.unit_model._tech_type}. Specify one of the following "
                f"cost methods in the unit's YAML file: {valid_methods}"
            )

    def cost_membrane(blk):
        """
        Get membrane cost based on membrane area and unit membrane costs
        as well as fixed operating cost for membrane replacement.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )
        # Get costing parameter sub-block for this technology
        mem_cost, rep_rate = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["membrane_cost", "membrane_replacement_rate"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        blk.variable_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit operation",
        )

        capex_expr = pyo.units.convert(
            mem_cost
            * pyo.units.convert(blk.unit_model.area, to_units=pyo.units.m**2),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        blk.variable_operating_cost_constraint = pyo.Constraint(
            expr=blk.variable_operating_cost
            == pyo.units.convert(
                rep_rate
                * mem_cost
                * pyo.units.convert(blk.unit_model.area, to_units=pyo.units.m**2),
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

    def cost_photothermal_membrane(blk):
        """
        General method for costing photothermal membrane.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        memb_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["membrane_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp["H2O"]
            / blk.unit_model.water_flux
            * memb_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

    def cost_CANDOP(blk):
        """
        General method for costing CANDO+P reactor.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        size_param, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["sizing_parameter", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * size_param * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_microbial_battery(blk):
        """
        General method for costing microbial battery treatment unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        sizing_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * sizing_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_in[t0].flow_mass_comp["filtration_media"],
            "filtration_media",
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp["filtration_media"],
            "filtration_media_disposal",
        )

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
        unit_capex = _get_tech_parameters(
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
        ZeroOrderCostingData._add_cost_factor(
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

    def cost_magprex(blk):
        """
        Method for costing Magprex reactor unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        HRT, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["HRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * HRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.MgCl2_flowrate[t0], "magnesium_chloride"
        )

    def cost_hrcs(blk):
        """
        Method for costing high-rate contact stabilization unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        SRT, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["SRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * SRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_centrifuge(blk):
        """
        Method for costing centrifuge unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        HRT, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["HRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * HRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.polymer_demand[t0], "polymer"
        )

    def cost_clarifier(blk, number_of_parallel_units=1):
        """
        General method for costing clarifiers. Costing is carried out
        using either the general_power_law form or the standard form which
        computes HRT, sizing costs, and chemical input costs.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        # Get cost method for this technology
        cost_method = _get_unit_cost_method(blk)
        valid_methods = ["cost_power_law_flow", "cost_HRCS_clarifier"]
        if cost_method == "cost_power_law_flow":
            ZeroOrderCostingData.cost_power_law_flow(blk, number_of_parallel_units)
        elif cost_method == "cost_HRCS_clarifier":
            # NOTE: number of units does not matter for cost_HRCS_clarifier
            #       as its a linear function of membrane area
            ZeroOrderCostingData.cost_HRCS_clarifier(blk)
        else:
            raise KeyError(
                f"{cost_method} is not a relevant cost method for "
                f"{blk.unit_model._tech_type}. Specify one of the following "
                f"cost methods in the unit's YAML file: {valid_methods}"
            )

    def cost_HRCS_clarifier(blk):
        """
        Method for costing a clarifier unit in a high-rate contact stabilization (HRCS) process.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        HRT, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["HRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * HRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.ferric_chloride_demand[t0], "ferric_chloride"
        )

    def cost_peracetic_acid(blk):
        """
        General method for costing peracetic acid water disinfection.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type,
            subtype=blk.unit_model.config.process_subtype,
        )

        # Get costing parameter sub-block for this technology
        sizing_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * sizing_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.disinfection_solution_flow_vol[t0], "disinfection_solution"
        )

    def cost_struvite_classifier(blk):
        """
        Method for costing struvite classifier unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        HRT, size_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["HRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties[t0].flow_vol * HRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties[t0].flow_mass_comp["struvite"],
            "struvite_product",
        )

    def cost_cloth_media_filtration(blk):
        """
        General method for costing cloth media filtration unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type,
            subtype=blk.unit_model.config.process_subtype,
        )

        # Get costing parameter sub-block for this technology
        sizing_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * sizing_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def cost_membrane_evaporator(blk):
        """
        General method for costing membrane evaporator unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        memb_cost = _get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["membrane_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.membrane_area * memb_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        ZeroOrderCostingData._add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        # TODO: Consider adding heat as a registered flow, since the inlet
        #       stream to the membrane evaporator must be heated to ~37 C.
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

    def _get_ozone_capital_cost(blk, A, B, C, D):
        """
        Generate expressions for capital cost of ozonation system.
        """
        t0 = blk.flowsheet().time.first()

        ln_Q = pyo.log(
            pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol
                / (pyo.units.m**3 / pyo.units.hour),
                to_units=pyo.units.dimensionless,
            )
        )
        dosage = pyo.units.convert(
            blk.unit_model.ozone_consumption[t0] / (pyo.units.mg / pyo.units.liter),
            to_units=pyo.units.dimensionless,
        )

        expr = pyo.units.convert(
            A + B * dosage + C * ln_Q + D * dosage * ln_Q,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        return expr

    def _get_uv_capital_cost(blk, A, B):
        """
        Generate expression for capital cost of UV reactor.
        """
        t0 = blk.flowsheet().time.first()

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.m**3 / pyo.units.hr,
        )

        E = pyo.units.convert(blk.unit_model.electricity[t0], to_units=pyo.units.kW)

        expr = pyo.units.convert(
            A * Q + B * E,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        return expr

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

    def _general_power_law_form(
        blk, A, B, sizing_term, factor=None, number_of_parallel_units=1
    ):
        """
        General method for building power law costing expressions.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A
            * (
                pyo.units.convert(sizing_term, to_units=pyo.units.dimensionless)
                / number_of_parallel_units
            )
            ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr *= number_of_parallel_units

        ZeroOrderCostingData._add_cost_factor(blk, factor)

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    unit_mapping = {
        ZeroOrderBase: cost_power_law_flow,
        AnaerobicMBRMECZO: cost_anaerobic_mbr_mec,
        ATHTLZO: cost_autothermal_hydrothermal_liquefaction,
        BrineConcentratorZO: cost_brine_concentrator,
        ChemicalAdditionZO: cost_chemical_addition,
        ChlorinationZO: cost_chlorination,
        ClarifierZO: cost_clarifier,
        ClothMediaFiltrationZO: cost_cloth_media_filtration,
        CoagulationFlocculationZO: cost_coag_and_floc,
        CofermentationZO: cost_cofermentation,
        ConstructedWetlandsZO: cost_constructed_wetlands,
        DeepWellInjectionZO: cost_deep_well_injection,
        DMBRZO: cost_dmbr,
        ElectroNPZO: cost_electrochemical_nutrient_removal,
        FixedBedZO: cost_fixed_bed,
        GACZO: cost_gac,
        LandfillZO: cost_landfill,
        MABRZO: cost_mabr,
        HTGZO: cost_hydrothermal_gasification,
        IonExchangeZO: cost_ion_exchange,
        IronManganeseRemovalZO: cost_iron_and_manganese_removal,
        MembraneEvaporatorZO: cost_membrane_evaporator,
        MetabZO: cost_metab,
        NanofiltrationZO: cost_nanofiltration,
        OzoneZO: cost_ozonation,
        OzoneAOPZO: cost_ozonation_aop,
        PeraceticAcidDisinfectionZO: cost_peracetic_acid,
        PumpElectricityZO: cost_pump_electricity,
        SaltPrecipitationZO: cost_supercritical_salt_precipitation,
        SedimentationZO: cost_sedimentation,
        StorageTankZO: cost_storage_tank,
        SuboxicASMZO: cost_suboxic_asm,
        SurfaceDischargeZO: cost_surface_discharge,
        UVZO: cost_uv,
        UVAOPZO: cost_uv_aop,
        EvaporationPondZO: cost_evaporation_pond,
        FilterPressZO: cost_filter_press,
        WellFieldZO: cost_well_field,
        PhotothermalMembraneZO: cost_photothermal_membrane,
        CANDOPZO: cost_CANDOP,
        MicrobialBatteryZO: cost_microbial_battery,
        VFARecoveryZO: cost_vfa_recovery,
        HRCSZO: cost_hrcs,
        MagprexZO: cost_magprex,
        CentrifugeZO: cost_centrifuge,
        StruviteClassifierZO: cost_struvite_classifier,
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
        pblock = getattr(blk.config.flowsheet_costing_block, blk.unit_model._tech_type)
    except AttributeError:
        # Parameter Block for this technology hasn't been added yet. Create it.
        pblock = pyo.Block()

        # Add block to FlowsheetCostingBlock
        blk.config.flowsheet_costing_block.add_component(
            blk.unit_model._tech_type, pblock
        )

        # Add subtype Set to Block
        pblock.subtype_set = pyo.Set()

        # Add required Vars
        for p in param_list:
            try:
                vobj = pyo.Var(
                    pblock.subtype_set,
                    units=getattr(
                        pyo.units, parameter_dict["capital_cost"][p]["units"]
                    ),
                )
                pblock.add_component(p, vobj)
            except KeyError:
                raise KeyError(
                    "Error when trying to retrieve costing parameter "
                    "for {p}. Please check the YAML "
                    "file for this technology for errors.".format(p=p)
                )

    # Check to see if required subtype is in subtype_set
    vlist = []
    if subtype not in pblock.subtype_set:
        # Need to add subtype and set Vars
        pblock.subtype_set.add(subtype)

        # Set vars
        for p in param_list:
            vobj = getattr(pblock, p)
            vobj[subtype].fix(
                float(parameter_dict["capital_cost"][p]["value"])
                * getattr(pyo.units, parameter_dict["capital_cost"][p]["units"])
            )
            vlist.append(vobj[subtype])
    else:
        for p in param_list:
            vobj = getattr(pblock, p)
            vlist.append(vobj[subtype])

    # add conditional for cases where there is only one parameter returned
    if len(vlist) == 1:
        return vlist[0]
    else:
        return tuple(x for x in vlist)


def _get_unit_cost_method(blk):
    """
    Get a specified cost_method if one is defined in the YAML file.
    This is meant for units with different cost methods between subtypes.
    """
    # Get parameter dict from database
    parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
        blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
    )

    if "cost_method" not in parameter_dict["capital_cost"]:
        raise KeyError(
            f"Costing for {blk.unit_model._tech_type} requires a cost_method argument, however "
            f"this was not defined for process sub-type {blk.unit_model.config.process_subtype}."
        )

    return parameter_dict["capital_cost"]["cost_method"]


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
