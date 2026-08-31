#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from math import isclose
import pyomo.environ as pyo

from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.core.expr.visitor import identify_variables

from idaes.core.base.costing_base import register_idaes_currency_units
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.base.costing_base import FlowsheetCostingBlockData
from idaes.models.unit_models import Mixer, HeatExchanger, Heater, CSTR

import idaes.logger as idaeslog

from watertap.core.util.misc import is_constant_up_to_units
from watertap.costing.unit_models.mixer import cost_mixer
from watertap.costing.unit_models.heat_exchanger import cost_heat_exchanger
from watertap.costing.unit_models.cstr import cost_cstr
from watertap.costing.unit_models.heater_chiller import cost_heater_chiller

_log = idaeslog.getLogger(__name__)


class WaterTAPCostingBlockData(FlowsheetCostingBlockData):
    """
    Base class for creating WaterTAP costing packages. Allows
    unit models to "self-register" their default costing methods,
    and for anonymous expressions in flow costs.
    """

    # Define default mapping of costing methods to unit models
    unit_mapping = {
        Mixer: cost_mixer,
        HeatExchanger: cost_heat_exchanger,
        CSTR: cost_cstr,
        Heater: cost_heater_chiller,
    }

    def register_currency_definitions(self):
        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

    def set_base_currency_base_period(self):
        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

    def add_LCOW(self, flow_rate, name="LCOW"):
        """
        Add Levelized Cost of Water (LCOW) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOW
            name (optional) - name for the LCOW variable (default: LCOW)
        """

        denominator = (
            pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / self.base_period)
            * self.utilization_factor
        )

        self.add_component(
            name,
            pyo.Expression(
                expr=(
                    self.total_capital_cost * self.capital_recovery_factor
                    + self.total_operating_cost
                )
                / denominator,
                doc=f"Levelized Cost of Water based on flow {flow_rate.name}",
            ),
        )

        c_units = self.base_currency
        t_units = self.base_period
        direct_capex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} direct capital expenditure by component",
            initialize=0 * c_units / pyo.units.m**3,
        )
        indirect_capex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} indirect capital expenditure by component",
            initialize=0 * c_units / pyo.units.m**3,
        )
        fixed_opex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} fixed operating expenditure by component",
            initialize=0 * c_units / pyo.units.m**3,
        )
        variable_opex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} variable operating expenditure by component",
            initialize=0 * c_units / pyo.units.m**3,
        )
        self.add_component(name + "_component_direct_capex", direct_capex_lcows)
        self.add_component(name + "_component_indirect_capex", indirect_capex_lcows)
        self.add_component(name + "_component_fixed_opex", fixed_opex_lcows)
        self.add_component(name + "_component_variable_opex", variable_opex_lcows)

        agg_direct_capex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} direct capital expenditure by unit type",
            initialize=0 * c_units / pyo.units.m**3,
        )
        agg_indirect_capex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} indirect capital expenditure by unit type",
            initialize=0 * c_units / pyo.units.m**3,
        )
        agg_fixed_opex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} fixed operating expenditure by unit type",
            initialize=0 * c_units / pyo.units.m**3,
        )
        agg_variable_opex_lcows = pyo.Expression(
            pyo.Any,
            doc=f"Levelized Cost of Water based on flow {flow_rate.name} variable operating expenditure by unit type",
            initialize=0 * c_units / pyo.units.m**3,
        )
        self.add_component(name + "_aggregate_direct_capex", agg_direct_capex_lcows)
        self.add_component(name + "_aggregate_indirect_capex", agg_indirect_capex_lcows)
        self.add_component(name + "_aggregate_fixed_opex", agg_fixed_opex_lcows)
        self.add_component(name + "_aggregate_variable_opex", agg_variable_opex_lcows)
        for u in self._registered_unit_costing:
            direct_capex_numerator = 0 * c_units / t_units
            indirect_capex_numerator = 0 * c_units / t_units
            fixed_opex_numerator = 0 * c_units / t_units
            variable_opex_numerator = 0 * c_units / t_units
            if hasattr(u, "capital_cost"):
                direct_capital_cost = pyo.units.convert(
                    u.direct_capital_cost, to_units=c_units
                )
                capital_cost = pyo.units.convert(u.capital_cost, to_units=c_units)
                # capital costs w/ recovery factor
                direct_capex_numerator += (
                    self.capital_recovery_factor * direct_capital_cost
                )
                capex_numerator = self.capital_recovery_factor * (
                    self.total_investment_factor * capital_cost
                )
                indirect_capex_numerator += capex_numerator - direct_capex_numerator

                # maintenance_labor_chemical_operating_cost,
                # part of total_fixed_operating_cost
                fixed_opex_numerator += (
                    self.maintenance_labor_chemical_factor * capital_cost
                )
            if hasattr(u, "fixed_operating_cost"):
                fixed_opex_numerator += pyo.units.convert(
                    u.fixed_operating_cost, to_units=c_units / t_units
                )
            if hasattr(u, "variable_operating_cost"):
                variable_opex_numerator += pyo.units.convert(
                    u.variable_operating_cost, to_units=c_units / t_units
                )
            direct_capex_lcows[u.unit_model.name] = direct_capex_numerator / denominator
            indirect_capex_lcows[u.unit_model.name] = (
                indirect_capex_numerator / denominator
            )
            fixed_opex_lcows[u.unit_model.name] = fixed_opex_numerator / denominator
            variable_opex_lcows[u.unit_model.name] = (
                variable_opex_numerator / denominator
            )

            unit_model_class = (
                u.unit_model.parent_component().process_block_class().__name__
            )
            agg_direct_capex_lcows[unit_model_class] += direct_capex_lcows[
                u.unit_model.name
            ]
            agg_indirect_capex_lcows[unit_model_class] += indirect_capex_lcows[
                u.unit_model.name
            ]
            agg_fixed_opex_lcows[unit_model_class] += fixed_opex_lcows[
                u.unit_model.name
            ]
            agg_variable_opex_lcows[unit_model_class] += variable_opex_lcows[
                u.unit_model.name
            ]

        for ftype, flows in self._registered_flows.items():
            # part of total_variable_operating_cost
            if ftype in agg_variable_opex_lcows:
                raise RuntimeError(
                    f"Found unit model named {ftype} but want to name flow {ftype} for {name+'_component_variable_opex'}"
                )
            if ftype not in self.used_flows:
                assert len(flows) == 0
                continue
            agg_variable_opex_lcows[ftype] = (
                self.aggregate_flow_costs[ftype] * self.utilization_factor / denominator
            )
            cost_var = getattr(self, f"{ftype}_cost")
            for flow_expr in flows:
                # part of total_variable_operating_cost
                flow_cost = pyo.units.convert(
                    flow_expr * cost_var, to_units=c_units / t_units
                )
                unit = self._find_flow_unit(flow_expr)
                if unit is not None:
                    if unit.name not in variable_opex_lcows:
                        unit_model_class = (
                            unit.parent_component().process_block_class().__name__
                        )
                        agg_variable_opex_lcows[
                            unit_model_class
                        ] += variable_opex_lcows[unit.name]
                    variable_opex_lcows[unit.name] += (
                        flow_cost * self.utilization_factor
                    ) / denominator
                    continue
                _log.warning(f"Could not find unique unit for flow {flow_expr}")
                flow_name = self._get_flow_name(flow_expr)
                if flow_name in variable_opex_lcows:
                    raise RuntimeError(
                        f"Found unit model named {flow_name} but need to name flow {flow_expr} for {name+'_aggregate_variable_opex'}"
                    )
                variable_opex_lcows[flow_name] = (
                    flow_cost * self.utilization_factor
                ) / denominator

    @staticmethod
    def _find_flow_unit(flow_expr):
        """
        Attempts to identify a unit model associated with the flow expression
        `flow_expr`. Returns the unit model if found, otherwise returns `None`.
        """
        variables = identify_variables(flow_expr, include_fixed=True)
        unit = WaterTAPCostingBlockData._find_flow_unit_variables(variables)
        return unit

    @staticmethod
    def _find_flow_unit_variables(variables):
        """
        Try to find a unit model uniquely associated with the iterable
        of variables `variables`.
        """
        unit = None
        for v in variables:
            blk = v.parent_block()
            while blk is not None:
                if UnitModelBlockData in blk.__class__.__mro__:
                    if unit is None:
                        unit = blk
                        break
                    elif blk is unit:
                        break
                    else:
                        # blk is unit model but not *this* unit,
                        # so we have an inconsistency
                        return None
                blk = blk.parent_block()
        return unit

    @staticmethod
    def _get_flow_name(flow_expr):
        """
        Get a standardized name for a specific flow expression.

        This will attempt to get the specific variable determining the flow,
        but failing that will return a string representation of the entire
        expression.
        """
        variables = list(identify_variables(flow_expr, include_fixed=False))
        if not variables:
            # No unfixed variables
            variables = list(identify_variables(flow_expr, include_fixed=True))
        # If we found a single variable for the flow_expr, just return the
        # variable name. Otherwise return a string representing the whole
        # expression.
        if len(variables) == 1:
            return str(variables[0])
        return str(flow_expr)

    def add_specific_energy_consumption(
        self, flow_rate, name="specific_energy_consumption"
    ):
        """
        Add specific energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific energy consumption
            name (optional) - the name of the Expression for the specific
                              energy consumption (default: specific_energy_consumption)
        """

        self.add_component(
            name,
            pyo.Expression(
                expr=self.aggregate_flow_electricity
                / pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / pyo.units.hr),
                doc=f"Specific energy consumption based on flow {flow_rate.name}",
            ),
        )
        self.add_flow_component_breakdown(
            "electricity", flow_rate, name=name, period=pyo.units.hr
        )

    def add_annual_water_production(self, flow_rate, name="annual_water_production"):
        """
        Add annual water production to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating annual water production
            name (optional) - name for the annual water production
                              Expression (default: annual_water_production)
        """
        self.add_component(
            name,
            pyo.Expression(
                expr=(
                    pyo.units.convert(
                        flow_rate, to_units=pyo.units.m**3 / self.base_period
                    )
                    * self.utilization_factor
                ),
                doc=f"Annual water production based on flow {flow_rate.name}",
            ),
        )

    def add_SEC(self, flow_rate):
        self.add_specific_energy_consumption(flow_rate, name="SEC")

    def add_electricity_intensity(self, flow_rate, name="electricity_intensity"):
        """
        Add calculation of overall electricity intensity to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating electricity intensity
            name (optional) - the name of the Expression for the specific
                              electrical intensity (default: specific_electrical_carbon_intensity)
        """
        self.add_specific_energy_consumption(flow_rate, name=name)

    def add_specific_electrical_carbon_intensity(
        self, flow_rate, name="specific_electrical_carbon_intensity"
    ):
        """
        Add specific electrical carbon intensity (kg_CO2eq/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific electrical carbon intensity
            name (optional) - the name of the Expression for the specific
                              carbon intensity (default: specific_electrical_carbon_intensity)
        """

        self.add_component(
            name,
            pyo.Expression(
                expr=self.aggregate_flow_electricity
                * self.electrical_carbon_intensity
                / pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / pyo.units.hr),
                doc=f"Specific electrical carbon intensity based on flow {flow_rate.name}",
            ),
        )
        self.add_flow_component_breakdown(
            "electricity",
            flow_rate,
            name=name,
            period=pyo.units.hr,
            multiplier=self.electrical_carbon_intensity,
        )

    def add_flow_component_breakdown(
        self,
        flow_name,
        flow_rate,
        name=None,
        period=None,
        multiplier=1.0,
    ):
        """
        Add per-component breakdowns for a registered flow type normalized by a flow rate over a period
        Args:
            flow_name - name of registered flow
            flow_rate - flow rate of water (volumetric) to be used for normalization
            name (optional)- base name for the component Expression (default is flow_name)
            period (optional) - time period for normalization (default is base_period)
            multiplier (optional) - multiplier for the flow (default is 1.0)
        """
        # NOTE: utilization_factor cancels in final expression but is included for convention

        if period is None:
            period = self.base_period

        if name is None:
            name = flow_name

        denominator = (
            pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / period)
            * self.utilization_factor
        )
        f_units = pyo.units.get_units(getattr(self, f"aggregate_flow_{flow_name}"))
        c_units = f_units * pyo.units.get_units(multiplier)
        expr_units = c_units / pyo.units.get_units(denominator)

        try:
            flows = self._registered_flows[flow_name]
        except KeyError:
            raise RuntimeError(f"Unrecognized flow_name {flow_name}.")

        specific_flow_consumption = pyo.Expression(
            pyo.Any,
            doc=f"Specific {flow_name} consumption by component",
            initialize=0.0 * period * c_units / pyo.units.m**3,
        )
        self.add_component(name + "_component", specific_flow_consumption)

        for flow_expr in flows:
            flow_std = pyo.units.convert(flow_expr, to_units=f_units)
            unit = self._find_flow_unit(flow_expr)
            if unit is not None:
                specific_flow_consumption[unit.name] += pyo.units.convert(
                    (flow_std * self.utilization_factor * multiplier) / denominator,
                    to_units=expr_units,
                )
                continue
            _log.warning(f"Could not find unique unit for flow {flow_expr}")
            flow_name = self._get_flow_name(flow_expr)
            specific_flow_consumption[flow_name] += pyo.units.convert(
                (flow_std * self.utilization_factor * multiplier) / denominator,
                to_units=expr_units,
            )

    def build_process_costs(self):
        """
        Build the common process costs to WaterTAP Costing Packages.
        The currency units should already be registered.

        The derived class should add constraints for total_capital_cost
        and total_operating_cost
        """
        self.total_capital_cost = pyo.Var(
            initialize=0,
            doc="Total capital cost of the process",
            units=self.base_currency,
        )
        self.total_operating_cost = pyo.Var(
            initialize=0,
            doc="Total operating cost of process per operating period",
            units=self.base_currency / self.base_period,
        )

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.total_investment_factor * self.aggregate_capital_cost
        )

        self.maintenance_labor_chemical_operating_cost = pyo.Expression(
            expr=self.maintenance_labor_chemical_factor * self.aggregate_capital_cost,
            doc="Maintenance-labor-chemical operating cost",
        )

        self.total_fixed_operating_cost = pyo.Expression(
            expr=self.aggregate_fixed_operating_cost
            + self.maintenance_labor_chemical_operating_cost,
            doc="Total fixed operating costs",
        )

        self.total_variable_operating_cost = pyo.Expression(
            expr=(
                (
                    self.aggregate_variable_operating_cost
                    + sum(self.aggregate_flow_costs[f] for f in self.used_flows)
                    * self.utilization_factor
                )
                if self.used_flows
                else self.aggregate_variable_operating_cost
            ),
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

    def build_global_params(self):
        """
        Build the global parameters common to WaterTAP Costing Packages.
        The currency units should already be registered.
        """

        self.register_currency_definitions()
        self.set_base_currency_base_period()

        self.utilization_factor = pyo.Var(
            initialize=0.9,
            doc="Plant capacity utilization [fraction of uptime]",
            units=pyo.units.dimensionless,
        )

        self.electricity_cost = pyo.Var(
            initialize=0.07,
            doc="Electricity cost",
            units=self.base_currency / pyo.units.kWh,
        )
        self.defined_flows["electricity"] = self.electricity_cost

        self.electrical_carbon_intensity = pyo.Var(
            initialize=0.475,
            doc="Grid carbon intensity [kgCO2_eq/kWh]",
            units=pyo.units.kg / pyo.units.kWh,
        )

        self.plant_lifetime = pyo.Var(
            initialize=30, units=self.base_period, doc="Plant lifetime"
        )

        self.wacc = pyo.Var(
            # consistent with a 30 year plant_lifetime
            # and a capital_recovery_factor of 0.1
            initialize=0.09307339771758532,
            units=pyo.units.dimensionless,
            doc="Weighted Average Cost of Capital (WACC)",
        )

        self.capital_recovery_factor = pyo.Var(
            initialize=0.1,
            units=pyo.units.year**-1,
            doc="Capital annualization factor [fraction of investment cost/year]",
        )

        # used in initialize_build to check for fixed variable consistency
        self._annualization_vars = (
            self.plant_lifetime,
            self.wacc,
            self.capital_recovery_factor,
        )

        self.capital_recovery_factor_constraint = pyo.Constraint(
            expr=self.capital_recovery_factor
            == (
                (self.wacc / pyo.units.year)
                / (1 - 1 / ((1 + self.wacc) ** (self.plant_lifetime / pyo.units.year)))
            )
        )

        self.TPEC = pyo.Var(
            initialize=3.4 * (2.0 / 1.65),
            doc="Total Purchased Equipment Cost (TPEC)",
            units=pyo.units.dimensionless,
        )

        self.TIC = pyo.Var(
            initialize=2.0,
            doc="Total Installed Cost (TIC)",
            units=pyo.units.dimensionless,
        )

        self.fix_all_vars()
        self.capital_recovery_factor.unfix()

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
        # handle wacc / plant_lifetime / capital_recovery_factor
        unfixed_vars = []
        for v in self._annualization_vars:
            if not v.fixed:
                unfixed_vars.append(v)
        if len(unfixed_vars) != 1:
            msg = "Exactly two of the variables "
            msg += ", ".join(v.name for v in self._annualization_vars)
            msg += " should be fixed and the other unfixed."
            raise RuntimeError(msg)
        calculate_variable_from_constraint(
            unfixed_vars[0], self.capital_recovery_factor_constraint
        )

    @staticmethod
    def add_cost_factor(blk, factor):
        """
        For a unit model costing block `blk`, adds `blk.cost_factor`,
        an expression pointing to the appropriate indirect capital
        cost multiplier, and `blk.direct_capital_cost`, which is a expression
        defined to be `blk.capital_cost / blk.cost_factor`. Valid strings for
        `factor` are "TIC" and "TPEC"; all others will result in an indirect
        capital cost factor of 1.

        Args:
            blk: an ideas.core.UnitModelCosting block
            factor: a string representing the cost factor to use
        """
        if factor == "TPEC":
            blk.cost_factor = pyo.Expression(expr=blk.costing_package.TPEC)
        elif factor == "TIC":
            blk.cost_factor = pyo.Expression(expr=blk.costing_package.TIC)
        else:
            blk.cost_factor = pyo.Expression(expr=1.0)
        blk.direct_capital_cost = pyo.Expression(
            expr=blk.capital_cost / blk.cost_factor
        )

    def _get_costing_method_for(self, unit_model):
        """
        Allow the unit model to register its default costing method,
        either through an attribute named "default_costing_method"
        or by naming the default costing method "default_costing_method"
        """
        if hasattr(unit_model, "default_costing_method"):
            return unit_model.default_costing_method
        return super()._get_costing_method_for(unit_model)

    def aggregate_costs(self):
        """
        This method aggregates costs from all the unit models and flows
        registered with this FlowsheetCostingBlock and creates aggregate
        variables for these on the FlowsheetCostingBlock that can be used for
        further process-wide costing calculations.

        The following costing variables are aggregated from all the registered
        UnitModelCostingBlocks (if they exist):

        * capital_cost,
        * direct_capital_cost,
        * fixed_operating_cost, and
        * variable_operating_cost

        Additionally, aggregate flow variables are created for all registered
        flow types along with aggregate costs associated with each of these.

        Args:
            None
        """
        super().aggregate_costs()
        c_units = self.base_currency

        @self.Expression(doc="Aggregation Expression for direct capital cost")
        def aggregate_direct_capital_cost(blk):
            e = 0
            for u in self._registered_unit_costing:
                # Allow for units that might only have a subset of cost Vars
                if hasattr(u, "direct_capital_cost"):
                    e += pyo.units.convert(u.direct_capital_cost, to_units=c_units)
                elif hasattr(u, "capital_cost"):
                    raise RuntimeError(
                        f"WaterTAP models with a capital_cost must also supply a direct_capital_cost. Found unit {u.unit_model} with `capital_cost` but no `direct_capital_cost`."
                    )

            return e

    def register_flow_type(self, flow_type, cost):
        """
        This method allows users to register new material and utility flows
        with the FlowsheetCostingBlock for use when costing flows.
        If `cost` is a constant (up to units), then this method creates a new
        `Var` on the FlowsheetCostingBlock named f`{flow_type}_cost`.
        Otherwise `cost` is a non-constant expression and this method will
        create a new `Expression` on the FlowsheetCostingBlock named
        f`{flow_type}_cost` whose value is fixed to `cost`.

        If a component named f`{flow_type}_cost` already exists on the
        FlowsheetCostingBlock, then an error is raised unless f`{flow_type}_cost`
        is `cost`. If f`{flow_type}_cost` is `cost`, no error is raised and
        the existing component f`{flow_type}_cost` is used to cost the flow.

        Args:
            flow_type: string name to represent flow type
            cost: a Pyomo expression with units representing the flow cost
        """

        flow_cost_name = flow_type + "_cost"
        current_flow_cost = self.component(flow_cost_name)
        if (current_flow_cost is None) and (not is_constant_up_to_units(cost)):
            cost_expr = pyo.Expression(expr=cost)
            self.add_component(flow_cost_name, cost_expr)
            super().register_flow_type(flow_type, cost_expr)
        else:
            # all other cases are handled in the base class
            super().register_flow_type(flow_type, cost)

    def plot_LCOW_breakdown(
        self,
        lcow_name="LCOW",
        fs_name="fs",
        by="component",
        separate_flows=False,
        relative=False,
        colormap="tab20",
        save_as=None,
        close_fig=False,
        component_labels={},
        tol=1e-5,
        dx=0,
        **kwargs,
    ):
        """
        This method creates a bar chart showing the breakdown of the LCOW calculation
        into its contributing components. The breakdown can be done by individual unit
        models or by unit model type, and flows can be shown separately or included as
        part of variable OPEX.

        Args:
            lcow_name: name of the LCOW Expression to plot (default: "LCOW")
            fs_name: name of the flowsheet (default: "fs")
            by: breakdown by individual "component" unit or "aggregate" unit class (default: "component")
            separate_flows: whether to separate material/energy flows from OPEX portion (default: False)
            relative: whether to plot as fraction of total LCOW (default: False)
            colormap: colormap to use for the plot (default: "tab20")
            save_as: filename to save the plot (default: None)
            close_fig: whether to close the figure after creation (default: False)
            component_labels: dictionary mapping flowsheet unit names to custom unit names (default: {})
            tol: tolerance for LCOW calculation check (default: 1e-5)

        Returns:
            fig, ax: matplotlib figure and axis objects for the plot
        """

        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        hatch_dict = {
            "CAPEX": "",
            "OPEX": "\\\\",
        }

        if self.find_component(lcow_name) is None:
            msg = f"{self.name} does not have a component named {lcow_name}."
            msg += " Use the add_LCOW method to add an LCOW Expression to the costing"
            msg += " block before calling plot_LCOW_breakdown."
            raise RuntimeError(msg)

        if by not in ("component", "aggregate"):
            raise ValueError("'by' argument must be either 'component' or 'aggregate'.")

        if relative:
            d = pyo.value(self.find_component(lcow_name)) * 1e-2
        else:
            d = 1

        capex_comp_exprs = [
            f"{lcow_name}_{by}_direct_capex",
            f"{lcow_name}_{by}_indirect_capex",
        ]

        # CAPEX portion of LCOW
        capex = dict(
            zip(
                [
                    k.replace(f"{fs_name}.", "")
                    for k in self.find_component(
                        f"{lcow_name}_{by}_direct_capex"
                    ).keys()
                ],
                [0] * len(self.find_component(f"{lcow_name}_{by}_direct_capex")),
            )
        )

        # CAPEX = direct + indirect
        for ce in capex_comp_exprs:
            expr = self.find_component(ce)
            for k, v in expr.items():
                capex[k.replace(f"{fs_name}.", "")] += pyo.value(v)

        # OPEX portion of LCOW
        opex = dict(
            zip(
                [
                    k.replace(f"{fs_name}.", "")
                    for k in self.find_component(f"{lcow_name}_{by}_fixed_opex").keys()
                ],
                [0] * len(self.find_component(f"{lcow_name}_{by}_fixed_opex")),
            )
        )
        # OPEX will always include fixed OPEX,
        # but may or may not include variable OPEX depending on `separate_flows`
        for k, v in self.find_component(f"{lcow_name}_{by}_fixed_opex").items():
            opex[k.replace(f"{fs_name}.", "")] += pyo.value(v)

        # Sort OPEX by value with lowest values first.
        # This is necessary visualize negative contributions (e.g., ERDs)
        opex = dict(sorted(opex.items(), key=lambda item: item[1], reverse=False))

        # Unique units contributing to LCOW as they appear on flowsheet
        fs_units = set(
            list(k.replace(f"{fs_name}.", "") for k in capex.keys())
            + list(k.replace(f"{fs_name}.", "") for k in opex.keys())
        )

        # We want this list to start with units that have negative OPEX values
        sorted_units = list(opex.keys()) + [
            u for u in capex.keys() if u not in opex.keys()
        ]
        if separate_flows:
            # Aggregate flows are plotted separately from variable OPEX
            flows = dict(zip(self.used_flows, [0] * len(self.used_flows)))

            for k, v in self.find_component(
                f"{lcow_name}_aggregate_variable_opex"
            ).items():
                if k not in self.used_flows:
                    continue
                flows[k] += pyo.value(v)

            hatch_dict["Flows"] = ".."
        else:
            # Aggregate flows are included as part of variable OPEX
            flows = {}
            for k, v in self.find_component(f"{lcow_name}_{by}_variable_opex").items():
                if k.replace(f"{fs_name}.", "") in fs_units:
                    opex[k.replace(f"{fs_name}.", "")] += pyo.value(v)

        # Unique components contributing to LCOW
        unique_components = list(flows.keys()) + list(sorted_units)

        cmap = plt.get_cmap(colormap)
        colors = [
            cmap(i / len(unique_components)) for i in range(len(unique_components))
        ]
        color_dict = dict(zip(unique_components, colors))

        # Start plotting
        fig, ax = plt.subplots()
        x = -1  # location of bar

        bottom = (
            0
            if all(v / d >= 0 for v in opex.values())
            else sum(v / d for v in opex.values() if v < 0) * 2
        )
        # Tracking minimum and maximum y-values
        bottoms = [bottom]

        # Create legend
        handles = list()
        labels = list()

        # For checking LCOW calculation
        lcow_check = 0

        # Starting from the bottom: add flows, then opex, then capex
        for f in flows.keys():
            v = flows[f] / d
            if f in component_labels.keys():
                label = component_labels[f]
            else:
                label = f
            ax.bar(
                [x],
                [v],
                bottom=bottom,
                hatch=hatch_dict["Flows"],
                color=color_dict[f],
                width=0.5,
                edgecolor="black",
            )
            bottom += v
            bottoms.append(bottom)

            ax.hlines(bottom, x - dx, x + dx, color="k", lw=2)
            x += dx

            lcow_check += v

            handles.append(Patch(facecolor=color_dict[f], edgecolor="k"))
            labels.append(label)

        for u in sorted_units:

            o = opex.get(u, 0) / d
            c = capex.get(u, 0) / d

            if u in component_labels.keys():
                label = component_labels[u]
            else:
                label = u

            ax.bar(
                [x],
                [abs(o)],
                bottom=bottom,
                hatch=hatch_dict["OPEX"],
                color=color_dict[u],
                edgecolor="black",
                width=0.5,
            )
            bottom += abs(o)

            ax.hlines(bottom, x - dx, x + dx, color="k", lw=2)
            x += dx

            ax.bar(
                [x],
                [c],
                bottom=bottom,
                hatch=hatch_dict["CAPEX"],
                color=color_dict[u],
                edgecolor="black",
                width=0.5,
            )
            bottom += abs(c)
            bottoms.append(bottom)

            ax.hlines(bottom, x - dx, x + dx, color="k", lw=2)
            x += dx

            lcow_check += c + o

            handles.append(Patch(facecolor=color_dict[u], edgecolor="k"))
            labels.append(label)

        # Reverse order to align with plot order (bottom to top)
        handles = handles[::-1]
        labels = labels[::-1]

        # Add hatch legends
        handles[0:0] = [
            Patch(facecolor="white", edgecolor="k", label=k, hatch=v)
            for k, v in hatch_dict.items()
        ]
        labels[0:0] = list(hatch_dict.keys())

        ax.set_axisbelow(True)
        ax.grid(visible=True)
        ax.legend(handles=handles, labels=labels)

        if dx == 0:
            ax.set_xlim(-1.5, 0.5)

        ax.set_xticks([])

        if bottoms[0] < 0:
            ax.hlines(0, -1.5, 0.5, color="k", lw=2, zorder=-100)
            ax.set_ylim(bottoms[0] * 1.1, ax.get_ylim()[1])

        ax.set_ylabel(
            f"{lcow_name} (\\$/m$^3$)" if not relative else f"Relative {lcow_name} (%)"
        )
        fig.tight_layout()

        if save_as is not None:
            fig.savefig(f"{save_as}.png", dpi=300, bbox_inches="tight")

        if not isclose(lcow_check, pyo.value(self.LCOW / d), rel_tol=tol):
            print(f"Calculated LCOW: {lcow_check}")
            print(f"Actual LCOW: {pyo.value(self.LCOW / d)}")
            raise ValueError(
                f"Calculated LCOW and actual LCOW differ by more than tolerance ({tol})"
            )

        if close_fig:
            plt.close(fig)

        return fig, ax, capex, opex

    def plot_SEC_breakdown(
        self,
        sec_name="specific_energy_consumption",
        fs_name="fs",
        relative=False,
        colormap="tab20",
        save_as=None,
        close_fig=False,
        tol=1e-5,
        component_labels={},
        **kwargs,
    ):
        """
        This method creates a bar chart showing the breakdown of the specific energy consumption calculation
        into its contributing components. The breakdown can be done by individual unit
        models or by unit model type.

        Args:
            sec_name: name of the specific energy consumption Expression to plot (default: "specific_energy_consumption")
            fs_name: name of the flowsheet (default: "fs")
            relative: whether to plot as fraction of total SEC (default: False)
            colormap: colormap to use for the plot (default: "tab20")
            save_as: filename to save the plot (default: None)
            close_fig: whether to close the figure after creation (default: False)
            component_labels: dictionary mapping flowsheet unit names to custom unit names (default: {})
            tol: tolerance for LCOW calculation check (default: 1e-5)
        """

        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        if self.find_component(sec_name) is None:
            msg = f"{self.name} does not have a component named {sec_name}."
            msg += " Use the add_specific_energy_consumption method to add a specific energy consumption Expression to the costing"
            msg += " block before calling plot_SEC_breakdown."
            raise RuntimeError(msg)

        if relative:
            d = pyo.value(self.find_component(sec_name)) * 1e-2
        else:
            d = 1

        sec_expr = self.find_component(f"{sec_name}_component")
        sec_comp = {
            k.replace(f"{fs_name}.", ""): pyo.value(v) for k, v in sec_expr.items()
        }
        # Sort so plotting starts with negative contributions
        sec_comp = dict(sorted(sec_comp.items(), key=lambda i: i[1], reverse=False))

        sorted_units = list(sec_comp.keys())

        cmap = plt.get_cmap(colormap)
        colors = [cmap(i / len(sorted_units)) for i in range(len(sorted_units))]
        color_dict = dict(zip(sorted_units, colors))

        # Start plotting
        fig, ax = plt.subplots()
        x = -1  # location of bar

        # Create legend
        handles = []
        labels = []

        sec_check = 0

        # Need to start at lowest point
        bottom = sum(ec for ec in sec_comp.values() if ec < 0) / d * 2
        bottoms = [bottom]

        for u in sorted_units:
            ec = sec_comp.get(u, 0) / d

            if u in component_labels.keys():
                label = component_labels[u]
            else:
                label = u

            ax.bar(
                [x],
                [abs(ec)],
                bottom=bottom,
                color=color_dict[u],
                edgecolor="black",
                width=0.5,
            )
            bottom += abs(ec)
            bottoms.append(bottom)
            sec_check += ec

            handles.append(Patch(facecolor=color_dict[u], edgecolor="k"))
            labels.append(label)

        # Reverse order to align with plot order (bottom to top)
        handles = handles[::-1]
        labels = labels[::-1]

        if bottoms[0] < 0:
            ax.hlines(0, -1.5, 0.5, color="k", lw=2, zorder=-1)
            ax.set_ylim(bottoms[0] * 1.1, ax.get_ylim()[1])

        ax.set_axisbelow(True)
        ax.grid(visible=True)
        ax.legend(handles=handles, labels=labels)
        ax.set_xlim(-1.5, 0.5)
        ax.set_xticks([])
        ax.set_ylabel(f"SEC (kWh/m$^3$)" if not relative else f"Relative SEC (%)")

        fig.tight_layout()

        if save_as is not None:
            fig.savefig(f"{save_as}.png", dpi=300, bbox_inches="tight")

        if not isclose(
            sec_check, pyo.value(self.find_component(sec_name) / d), rel_tol=tol
        ):
            print(f"Calculated SEC: {sec_check}")
            print(f"Actual SEC: {pyo.value(self.find_component(sec_name) / d)}")
            raise ValueError(
                f"Calculated SEC and actual SEC differ by more than tolerance ({tol})"
            )

        if close_fig:
            plt.close(fig)

        return fig, ax


@declare_process_block_class("WaterTAPCosting")
class WaterTAPCostingData(WaterTAPCostingBlockData):
    def build_global_params(self):

        # Build flowsheet level costing components
        # These are the global parameters
        self.total_investment_factor = pyo.Var(
            initialize=1.0,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.maintenance_labor_chemical_factor = pyo.Var(
            initialize=0.03,
            doc="Maintenance-labor-chemical factor [fraction of equipment cost/year]",
            units=pyo.units.year**-1,
        )

        super().build_global_params()


@declare_process_block_class("WaterTAPCostingDetailed")
class WaterTAPCostingDetailedData(WaterTAPCostingBlockData):
    def build_global_params(self):
        """
        To minimize overhead, only create global parameters for now.
        Unit-specific parameters will be added as sub-Blocks on a case-by-case
        basis as a unit of that type is costed.
        """

        # Costing factors
        self.land_cost_percent_FCI = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Land cost as % FCI",
            initialize=0.0,
        )
        self.working_capital_percent_FCI = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Working capital as % FCI",
            initialize=0.0,
        )
        self.salaries_percent_FCI = pyo.Var(
            units=1 / self.base_period,
            doc="Salaries as % FCI",
            initialize=0.001 * (0.03 / 0.0149),
        )
        self.benefit_percent_of_salary = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Benefits as % salaries",
            initialize=0.9,
        )
        self.maintenance_costs_percent_FCI = pyo.Var(
            units=1 / self.base_period,
            doc="Maintenance and contingency costs as % FCI",
            initialize=0.008 * (0.03 / 0.0149),
        )
        self.laboratory_fees_percent_FCI = pyo.Var(
            units=1 / self.base_period,
            doc="Laboratory fees as % FCI",
            initialize=0.003 * (0.03 / 0.0149),
        )
        self.insurance_and_taxes_percent_FCI = pyo.Var(
            units=1 / self.base_period,
            doc="Insurance and taxes as % FCI",
            initialize=0.002 * (0.03 / 0.0149),
        )

        self.total_investment_factor = pyo.Expression(
            expr=1.0 + self.working_capital_percent_FCI + self.land_cost_percent_FCI,
            doc="Total investment factor [investment cost/equipment cost]",
        )

        self.maintenance_labor_chemical_factor = pyo.Expression(
            expr=self.salaries_percent_FCI
            + self.benefit_percent_of_salary * self.salaries_percent_FCI
            + self.maintenance_costs_percent_FCI
            + self.laboratory_fees_percent_FCI
            + self.insurance_and_taxes_percent_FCI,
            doc="Maintenance-labor-chemical factor [fraction of equipment cost/year]",
        )

        super().build_global_params()

    def build_process_costs(self):
        """
        Calculating process wide costs.
        """
        # add the "basic" process costs
        super().build_process_costs()

        # now the detailed breakdowns
        # Other capital costs
        self.land_cost = pyo.Expression(
            expr=self.land_cost_percent_FCI * self.aggregate_capital_cost,
            doc="Land costs - based on aggregate capital costs",
        )
        self.working_capital = pyo.Expression(
            expr=self.working_capital_percent_FCI * self.aggregate_capital_cost,
            doc="Working capital - based on aggregate capital costs",
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
