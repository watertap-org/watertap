#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
This module contains the functions needed for the construction of flexible
desalination flowsheet
"""

from idaes.apps.grid_integration import OperationModel
from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Param,
    Var,
    units as pyunits,
)
from watertap.flowsheets.flex_desal import params as um_params
from watertap.flowsheets.flex_desal import unit_models as um


def add_operational_cost_expressions(blk, params: um_params.FlexDesalParams):
    """
    Adds cost expressions to the flowsheet
    """
    # Water revenue
    blk.water_revenue = Expression(
        expr=(
            params.product_water_price
            * blk.posttreatment.product_flowrate
            * params.timestep_hours
        ),
        doc="Revenue generated from product water",
    )

    # Customer cost
    blk.customer_cost = Param(
        initialize=0,
        mutable=True,
        doc="Fixed customer cost",
    )

    # Demand response revenue
    blk.demand_response_price = Param(
        initialize=0, mutable=True, doc="Demand-response prices"
    )
    blk.baseline_power = Param(
        initialize=100, mutable=True, doc="Baseline power requirement"
    )
    blk.demand_response_revenue = Expression(
        expr=blk.demand_response_price
        * (blk.baseline_power - blk.power_from_grid)
        * params.timestep_hours,
        doc="Revenue generated from demand response",
    )

    # Cost of emissions
    blk.emissions_intensity = Param(
        initialize=0, mutable=True, units=pyunits.kg / pyunits.kWh
    )
    blk.emissions_cost = Expression(
        expr=(
            blk.emissions_intensity
            * blk.power_from_grid
            * params.timestep_hours
            * params.emissions_cost
            / 907.185  # Conversion factor: $/ton to $/kg
        ),
        doc="Cost associated with carbon emissions",
    )

    # Cost of energy
    blk.LMP = Param(
        initialize=0,
        mutable=True,
        doc="Locational marginal price of electricity [$/kWh]",
    )
    blk.energy_cost = Expression(
        expr=blk.LMP * blk.power_from_grid * params.timestep_hours,
        doc="Cost of electricity purchased from the grid",
    )

    # Demand cost parameters
    blk.fixed_demand_rate = Param(
        initialize=0,
        mutable=True,
        doc="Constant demand tariff",
    )
    blk.variable_demand_rate = Param(
        initialize=0,
        mutable=True,
        doc="Variable demand tariff",
    )


def build_desal_flowsheet(blk, params: um_params.FlexDesalParams):
    """
    Builds a flowsheet instance of the entire desalination process

    Parameters
    ----------
    blk : Block
        Pyomo Block instance

    params : object
        Object containing model parameters
    """
    # Build units
    blk.intake = OperationModel(
        model_func=um.intake_operation_model,
        model_args={"params": params.intake},
    )
    blk.bypass_pretreatment_flow = Var(
        within=NonNegativeReals,
        doc="Flowrate bypassed to brine discharge due to pretreatment shutdown",
    )
    blk.pretreatment = OperationModel(
        model_func=um.pretreatment_operation_model,
        model_args={"params": params.pretreatment},
    )
    blk.reverse_osmosis = OperationModel(
        model_func=um.reverse_osmosis_operation_model,
        model_args={"params": params.ro},
    )
    blk.posttreatment = OperationModel(
        model_func=um.posttreatment_operation_model,
        model_args={"params": params},
    )
    blk.brine_discharge = OperationModel(
        model_func=um.brine_discharge_operation_model,
        model_args={"params": params},
    )

    # Flowsheet connections
    blk.arc_intake_pretreatment = Constraint(
        expr=blk.intake.product_flowrate
        == blk.pretreatment.feed_flowrate + blk.bypass_pretreatment_flow,
        doc="intake-pretreatment mass balance",
    )
    blk.suppress_pretreatment_bypass = Constraint(
        expr=blk.bypass_pretreatment_flow
        <= (1 - blk.pretreatment.op_mode) * params.intake.nominal_flowrate
    )
    blk.arc_pretreatment_ro = Constraint(
        expr=blk.pretreatment.product_flowrate == blk.reverse_osmosis.feed_flowrate,
        doc="pretreatment-reverse_osmosis mass balance",
    )
    blk.arc_ro_posttreatment = Constraint(
        expr=blk.reverse_osmosis.product_flowrate == blk.posttreatment.feed_flowrate,
        doc="reverse_osmosis-posttreatment mass balance",
    )
    blk.calculate_brine_discharge = Constraint(
        expr=blk.brine_discharge.feed_flowrate
        == (
            blk.intake.reject_flowrate
            + blk.bypass_pretreatment_flow
            + blk.pretreatment.reject_flowrate
            + blk.reverse_osmosis.reject_flowrate
            + blk.reverse_osmosis.leftover_flow
            + blk.posttreatment.reject_flowrate
        ),
        doc="Computes the total inflow to brine discharge",
    )

    blk.num_skids_online = Expression(
        expr=sum(blk.reverse_osmosis.ro_skid[:].op_mode),
        doc="Calculates the number of skids operating at time t",
    )

    blk.net_power_consumption = Expression(
        expr=blk.intake.power_consumption
        + blk.pretreatment.power_consumption
        + blk.reverse_osmosis.power_consumption
        + blk.posttreatment.power_consumption
        + blk.brine_discharge.power_consumption,
        doc="Net power consumed from the grid",
    )

    if params.include_onsite_solar:
        blk.power_generation = OperationModel(
            model_func=um.power_generation_operation_model,
            model_args={"params": params},
        )
        blk.net_power_consumption += -blk.power_generation.power_utilized

    if params.include_battery:
        blk.battery = OperationModel()
        blk.net_power_consumption += blk.battery.power_charge - blk.battery.discharge

    # Power purchased from the grid
    blk.power_from_grid = Var(
        within=NonNegativeReals,
        units=pyunits.kW,
        doc="Total power purchased from the grid",
    )
    blk.overall_power_balance = Constraint(
        expr=blk.power_from_grid == blk.net_power_consumption,
        doc="Computes the total power purchased from the grid",
    )

    # Add cost expressions
    add_operational_cost_expressions(blk, params)


def add_delayed_startup_constraints(m):
    """Adds the delayed startup constraints to the model"""
    params: um_params.FlexDesalParams = m.params

    # "Shutdown" post-treatment unit if RO startup is initiated
    @m.Constraint(m.period.index_set())
    def posttreatment_unit_commitment(blk, d, t):
        indices = [(d, t - i) for i in range(params.ro.startup_delay) if t - i > 0]
        return (1 - blk.period[d, t].posttreatment.op_mode) == sum(
            blk.period[p].reverse_osmosis.ro_skid[1].startup for p in indices
        )

    # Brine pump must operate if RO startup is initiated
    @m.Constraint(m.period.index_set())
    def brine_pump_unit_commitment(blk, d, t):
        indices = [(d, t - i) for i in range(params.ro.startup_delay) if t - i > 0]
        return blk.period[d, t].brine_discharge.op_mode == sum(
            blk.period[p].reverse_osmosis.ro_skid[1].startup for p in indices
        )


def add_demand_and_fixed_costs(m):
    """Adds variables and expressions/constraints for demand and fixed costs"""

    params: um_params.FlexDesalParams = m.params
    m.fixed_demand_cost = Var(
        within=NonNegativeReals,
        doc="Total fixed demand charge value for the entire horizon",
    )
    m.variable_demand_cost = Var(
        within=NonNegativeReals,
        doc="Total variable demand charge value for the entire time horizon",
    )
    m.fixed_monthly_cost = Var(
        within=NonNegativeReals,
        doc="Total customer cost for the entire time horizon",
    )

    @m.Constraint(m.period.index_set())
    def calculate_fixed_demand_cost(blk, d, t):
        return (
            blk.fixed_demand_cost
            >= blk.period[d, t].fixed_demand_rate
            * blk.period[d, t].power_from_grid
            * params.num_months
        )

    @m.Constraint(m.period.index_set())
    def calculate_variable_demand_cost(blk, d, t):
        return (
            blk.variable_demand_cost
            >= blk.period[d, t].variable_demand_rate
            * blk.period[d, t].power_from_grid
            * params.num_months
        )

    m.calculate_fixed_monthly_cost = Constraint(
        expr=m.fixed_monthly_cost == params.fixed_monthly_cost * params.num_months
    )


def add_useful_expressions(m):
    """Defines useful expressions for custom objective functions"""

    m.total_water_revenue = Expression(expr=sum(m.period[:, :].water_revenue))
    m.total_demand_response_revenue = Expression(
        expr=sum(m.period[:, :].demand_response_revenue)
    )
    m.total_emissions_cost = Expression(expr=sum(m.period[:, :].emissions_cost))


def constrain_water_production(m, baseline_production: float = None):
    """Constrains the total water production rate"""

    params: um_params.FlexDesalParams = m.params
    if baseline_production is not None:
        m.curtailment_fraction = Param(
            initialize=params.curtailment_fraction,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Fraction of water production that is curtailed",
        )

        m.baseline_production = Param(
            initialize=baseline_production,
            mutable=True,
            units=pyunits.m**3,
            doc="Baseline water production",
        )

        m.water_production_target = Constraint(
            expr=m.total_water_production
            >= m.baseline_production * (1 - m.curtailment_fraction)
        )

    elif params.annual_production_AF is not None:
        # Convert production rate from acre-ft/year to m^3/year
        annual_production_m3 = params.annual_production_AF * 1233.48
        m.production_target_abs = Param(
            initialize=annual_production_m3 / 365 * params.num_days,
            mutable=True,
            units=pyunits.m**3,
            doc="Absolute water production target",
        )

        m.water_production_target = Constraint(
            expr=m.total_water_production >= m.production_target_abs
        )

    else:
        raise ValueError("Water production targets not specified in params")
