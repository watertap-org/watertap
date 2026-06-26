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
Copied from the flexible desalination flowsheet, but altered ro surrogate.
"""

from idaes.apps.grid_integration import OperationModel
from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Param,
    Var,
    Binary,
    units as pyunits,
    value,
)

from watertap.flowsheets.flex_desal import params as um_params
from watertap.flowsheets.flex_desal import unit_models as um
from watertap.flowsheets.flex_desal.wrd_unit_models import (
    wrd_reverse_osmosis_operation_model,
    wrd_uf_operation_model,
)


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
        initialize=1062, mutable=True, doc="Baseline power requirement"
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


def build_wrd_desal_flowsheet(blk, params: um_params.FlexDesalParams):
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
    # pretreatment in this case refers to the UF
    blk.pretreatment = OperationModel(
        model_func=wrd_uf_operation_model,
        model_args={"params": params.wrd_uf},
    )
    blk.reverse_osmosis = OperationModel(
        model_func=wrd_reverse_osmosis_operation_model,
        model_args={"params": params.wrd_ro},
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

        blk.excess_solar_power = Var(
            within=NonNegativeReals,
            units=pyunits.kW,
            doc="Excess solar power that is generated but not utilized",
        )
        blk.net_power_consumption += (
            blk.excess_solar_power
        )  # Slack variable so that when plant power is zero, the solar has somewhere to go.

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
        indices = [(d, t - i) for i in range(params.wrd_ro.startup_delay) if t - i > 0]
        return (1 - blk.period[d, t].posttreatment.op_mode) == sum(
            blk.period[p].reverse_osmosis.ro_skid[1].startup for p in indices
        )

    # Brine pump must operate if RO startup is initiated
    @m.Constraint(m.period.index_set())
    def brine_pump_unit_commitment(blk, d, t):
        indices = [(d, t - i) for i in range(params.wrd_ro.startup_delay) if t - i > 0]
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


def add_flow_costs(m):
    """Adds expressions for feed and brine costs"""

    # If the brine discharge in being operated, then plant start-up is occurring, and the feed is being recirculated, so there is no cost.
    # Arguably this could be applied to the brine as well but that's tbd
    m.total_feed_cost = Expression(
        expr=sum(
            m.period[d, t].intake.feed_cost
            * (1 - m.period[d, t].brine_discharge.op_mode)
            for d, t in m.period.index_set()
        )
        * m.params.timestep_hours,
        doc="Total cost of feed water over the time horizon ($)",
    )

    m.total_brine_cost = Expression(
        expr=sum(m.period[:, :].brine_discharge.brine_cost) * m.params.timestep_hours,
        doc="Total cost of brine discharge over the time horizon ($)",
    )

    m.total_chemical_cost = Expression(
        expr=m.params.timestep_hours
        * (
            sum(m.period[:, :].intake.chemical_cost)
            + sum(m.period[:, :].posttreatment.chemical_cost)
        ),
        doc="Total cost of chemicals over the time horizon ($)",
    )


def add_flow_changes_penalty(m):
    # Add binary variables to track flowrate changes between consecutive time steps
    m.flow_changed = Var(
        m.set_days,
        m.set_time,
        range(1, m.params.wrd_ro.num_ro_skids + 1),
        within=Binary,
        doc="Binary variable: 1 if RO skid flowrate changes from previous time step, 0 otherwise",
    )

    # Add binary variables to track UF pump flowrate changes between consecutive time steps
    m.uf_flow_changed = Var(
        m.set_days,
        m.set_time,
        range(1, m.params.wrd_uf.num_uf_pumps + 1),
        within=Binary,
        doc="Binary variable: 1 if UF pump flowrate changes from previous time step, 0 otherwise",
    )

    # Add constraints to detect flowrate changes
    # Big-M value used for the constraint, based on max flowrate
    big_M = 2 * m.params.wrd_ro.maximum_flowrate
    big_M_uf = 2 * m.params.wrd_uf.maximum_flowrate

    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def track_flow_changes(m_blk, d, t, i):
        # Skip first time step (no previous time step to compare)
        if t == 1:
            return Constraint.Skip

        # Get current and previous period flowrates
        current_flow = m_blk.period[d, t].reverse_osmosis.ro_skid[i].feed_flowrate
        prev_flow = m_blk.period[d, t - 1].reverse_osmosis.ro_skid[i].feed_flowrate

        flow_diff = current_flow - prev_flow

        return m_blk.flow_changed[d, t, i] * big_M >= flow_diff

    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def track_flow_changes_neg(m_blk, d, t, i):
        # Skip first time step
        if t == 1:
            return Constraint.Skip

        current_flow = m_blk.period[d, t].reverse_osmosis.ro_skid[i].feed_flowrate
        prev_flow = m_blk.period[d, t - 1].reverse_osmosis.ro_skid[i].feed_flowrate

        flow_diff = prev_flow - current_flow

        # Capture negative direction of flow change
        return m_blk.flow_changed[d, t, i] * big_M >= flow_diff

    # Track UF pump flowrate changes
    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_uf.num_uf_pumps + 1))
    def track_uf_flow_changes(m_blk, d, t, i):
        # Skip first time period of first day (no previous period to compare)
        if t == 1:
            return Constraint.Skip

        # Get current and previous period flowrates
        current_flow = m_blk.period[d, t].pretreatment.uf_pumps[i].feed_flowrate
        prev_flow = m_blk.period[d, t - 1].pretreatment.uf_pumps[i].feed_flowrate

        # If flows are different, uf_flow_changed must be 1
        flow_diff = current_flow - prev_flow

        return m_blk.uf_flow_changed[d, t, i] * big_M_uf >= flow_diff

    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_uf.num_uf_pumps + 1))
    def track_uf_flow_changes_neg(m_blk, d, t, i):
        # Skip first time period of first day
        if t == 1:
            return Constraint.Skip

        current_flow = m_blk.period[d, t].pretreatment.uf_pumps[i].feed_flowrate
        prev_flow = m_blk.period[d, t - 1].pretreatment.uf_pumps[i].feed_flowrate

        flow_diff = prev_flow - current_flow

        # Capture negative direction of flow change
        return m_blk.uf_flow_changed[d, t, i] * big_M_uf >= flow_diff

    # The penalty is simply the number of flow changes multiplied by a scaling factor
    m.flow_changes_penalty = Expression(
        expr=50  # Scaling factor (adjust as needed). This is equivalent to a cost of $50 per flowrate change
        * (
            sum(
                sum(
                    sum(m.flow_changed[d, t, i] for t in m.set_time) for d in m.set_days
                )
                for i in range(1, m.params.wrd_ro.num_ro_skids + 1)
            )
            + sum(
                sum(
                    sum(m.uf_flow_changed[d, t, i] for t in m.set_time)
                    for d in m.set_days
                )
                for i in range(1, m.params.wrd_uf.num_uf_pumps + 1)
            )
        )
    )


def calculate_replacement_costs(m):
    """Adds expression for the replacement cost"""
    params: um_params.WRD_ROParams = m.params.wrd_ro

    m.degree_of_flex = Expression(
        expr=sum(
            m.period[d, t].reverse_osmosis.ro_skid[i].shutdown
            for d in m.set_days
            for t in m.set_time
            for i in range(1, params.num_ro_skids + 1)
        )
        / (
            2 * m.params.num_days * params.num_ro_skids
        ),  # 2 is arbitrary. Means that 2 shutdowns per day per skid would yield a raw_degree_of_flex of 1
        doc="Constraint to compute raw flexibility metric",
    )

    if params.replacement_types:
        for i, replacement_type in enumerate(params.replacement_types):
            # Create a variable for that replacement type
            setattr(
                m,
                f"replacement_cost_{replacement_type}",
                Param(
                    within=NonNegativeReals,
                    initialize=params.replacement_costs[i],
                    doc=f"Replacement cost for {replacement_type}",
                ),
            )

        m.total_replacement_cost = Expression(
            expr=(
                sum(
                    getattr(m, f"replacement_cost_{replacement_type}")
                    / (
                        params.replacement_lifetimes[i]
                        * (
                            1
                            - params.replacement_max_flex_penalty[i] * m.degree_of_flex
                        )
                    )
                    * m.params.num_months
                    / 12
                    for i, replacement_type in enumerate(params.replacement_types)
                )
            ),
            doc="Total replacement costs annualized over the time horizon",
        )


def calculate_flexibility_metrics(
    m,
    baseline_power=1000,
    baseline_electricity_cost=100000,
    baseline_replacement_cost=993,
):
    """Calculates flexibility metrics based on model results. Should be called after solving the model."""

    maximum_power = max(
        value(m.period[d, t].power_from_grid) for d, t in m.period.index_set()
    )
    print(f"Maximum power (kW): {maximum_power:.2f}")

    # Discharging energy capacity
    discharge_energy_capacity = (
        sum(
            max(0, baseline_power - value(m.period[d, t].power_from_grid))
            for d, t in m.period.index_set()
        )
        * m.params.timestep_hours
    )
    print(f"Discharge energy capacity (kWh): {discharge_energy_capacity:.2f}")
    discharge_time = sum(
        (value(m.period[d, t].power_from_grid) < baseline_power)
        * m.params.timestep_hours
        for d, t in m.period.index_set()
    )

    # Charging energy capacity
    charge_energy_capacity = (
        sum(
            max(0, value(m.period[d, t].power_from_grid) - baseline_power)
            for d, t in m.period.index_set()
        )
        * m.params.timestep_hours
    )
    print(f"Charge energy capacity (kWh): {charge_energy_capacity:.2f}")
    charge_time = sum(
        (value(m.period[d, t].power_from_grid) > baseline_power)
        * m.params.timestep_hours
        for d, t in m.period.index_set()
    )

    # Power Capacities
    if discharge_time > 0:
        discharge_power_capacity = discharge_energy_capacity / discharge_time
    else:
        discharge_power_capacity = float("nan")

    if charge_time > 0:
        charge_power_capacity = charge_energy_capacity / charge_time
    else:
        charge_power_capacity = float("nan")

    if discharge_energy_capacity == 0:
        LVOF = float("nan")
    else:
        LVOF = (
            (
                baseline_electricity_cost
                - value(
                    m.total_energy_cost
                    + m.total_demand_cost
                    + m.total_customer_cost
                    - m.total_demand_response_revenue
                )
            )
            - (baseline_replacement_cost - value(m.total_replacement_cost))
        ) / (discharge_energy_capacity)
    print(f"Levelized Cost of Flexibility ($/kWh): {LVOF:.2f}")
    # Levelized Cost of Flexibility. Only flexibility costs are the replacement costs

    m.maximum_power = Var(initialize=maximum_power)
    m.energy_capacity = Var(initialize=discharge_energy_capacity)
    m.power_capacity = Var(initialize=discharge_power_capacity)
    m.LVOF = Var(initialize=LVOF)


def add_useful_expressions(m):
    """Defines useful expressions for custom objective functions"""

    m.total_water_revenue = Expression(expr=sum(m.period[:, :].water_revenue))
    m.total_demand_response_revenue = Expression(
        expr=sum(m.period[:, :].demand_response_revenue)
    )
    m.total_emissions_cost = Expression(expr=sum(m.period[:, :].emissions_cost))


def begin_and_end_constraint(m):
    """Force RO train 1 op_mode to match between first and last timesteps."""
    period_points = list(m.period.index_set())
    if not period_points:
        return

    first_point = period_points[0]
    last_point = period_points[-1]

    @m.Constraint()
    def match_train_1_at_start_and_end(blk):
        return (
            blk.period[first_point].reverse_osmosis.ro_skid[1].op_mode
            == blk.period[last_point].reverse_osmosis.ro_skid[1].op_mode
        )


### NOT USED IN THE TUTORIAL EXAMPLE - That might mean they aren't being tested? ###
def fix_operations_for_first_four_days(m, peak_hours=None):
    """Fix all RO trains to expected behavior for first four days. This could be some part of an initialization that improve solve times."""
    for d, p in m.period:
        if p <= 4 * 24:  # Assuming hourly time steps
            if p <= 2:
                # Avoiding constraint that plant has to be on at first time step.
                pass
            elif peak_hours is not None and peak_hours[p]:
                m.period[d, p].reverse_osmosis.ro_skid[4].op_mode.fix(
                    0
                )  # 4th skid off during peak hours
            else:
                # Ensure plant is on during non-peak hours
                m.period[d, p].reverse_osmosis.ro_skid[1].op_mode.fix(
                    1
                )  # Plant must be on


def add_delayed_shutdown_constraints(m):
    # Consider implmenting with the add_ramping_limits from IDAES price_taker_model
    """Adds the delayed shutdown constraints to the model"""
    params: um_params.FlexDesalParams = m.params

    """Specific to WRD where the planned shutdowns seem to occur over a period 60-100 minutes, meaning it's 
    not realistic to have all trains go from on to off in same hour. This says 30 mins per train about right"""

    @m.Constraint(m.period.index_set())
    def posttreatment_unit_commitment_shutdown(blk, d, t):
        return (
            sum(
                blk.period[d, t].reverse_osmosis.ro_skid[i].shutdown
                for i in range(1, params.wrd_ro.num_ro_skids + 1)
            )
            <= 2
        )


def add_working_hours_constraint(m):
    "Prevents shutdowns and startups duing nonworking hours (e.g., 6pm-8am)"

    @m.Constraint(m.period.index_set())
    def prevent_shutdowns(blk, d, t):
        if (
            t - 1
        ) % 24 in m.params.nonworking_hours:  # Assuming time index is in hours and starts at 1
            return (
                sum(
                    blk.period[d, t].reverse_osmosis.ro_skid[i].shutdown
                    for i in range(1, blk.params.wrd_ro.num_ro_skids + 1)
                )
                == 0
            )
        else:
            return Constraint.Skip

    @m.Constraint(m.period.index_set())
    def prevent_startups(blk, d, t):
        if (
            t - 1
        ) % 24 in m.params.nonworking_hours:  # Assuming time index is in hours and starts at 1
            return (
                sum(
                    blk.period[d, t].reverse_osmosis.ro_skid[i].startup
                    for i in range(1, blk.params.wrd_ro.num_ro_skids + 1)
                )
                == 0
            )
        else:
            return Constraint.Skip


def add_maximum_shutdowns(m):
    """Adds rolling 24-hour shutdown limits over the full period index."""
    params: um_params.FlexDesalParams = m.params

    window_steps = max(1, int(round(24 / params.timestep_hours)))
    period_points = list(m.period.index_set())
    num_windows = max(0, len(period_points) - window_steps + 1)

    @m.Constraint(range(num_windows))
    def max_shutdowns_per_24h_window(blk, w):
        return (
            sum(
                blk.period[period_points[k]].reverse_osmosis.ro_skid[1].shutdown
                for k in range(w, w + window_steps)
            )
            <= params.max_daily_shutdowns
        )


def restrict_flexible_trains(m, num_flexible_trains):
    ro_skids = sorted(list(m.period[1, 1].reverse_osmosis.set_ro_skids))
    n_ro_skids = len(ro_skids)

    if num_flexible_trains < 0 or num_flexible_trains > n_ro_skids:
        raise ValueError(
            "Invalid num_flexible_trains "
            f"'{num_flexible_trains}'. Valid range is [0, {n_ro_skids}]."
        )

    non_flexible_skids = ro_skids[: n_ro_skids - num_flexible_trains]

    for p in m.period:
        for skid in non_flexible_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.startup.fix(0)
            ro_skid.shutdown.fix(0)


# These ones might not need to be included at all
def repeat_weekdays(m):
    """Ensures operations during first four days are repeated"""

    detla_time = (
        24 / m.params.timestep_hours
    )  # Assuming time index is in hours and starts at 1

    @m.Constraint(m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def repeat_weekday_flowrates(blk, t, i):
        if t >= detla_time + 1 and t <= 4 * detla_time:  # Compare day 2-4 to day 1
            return (
                blk.period[1, t].reverse_osmosis.ro_skid[i].feed_flowrate
                == blk.period[1, t].reverse_osmosis.ro_skid[i].feed_flowrate
            )
        else:
            return Constraint.Skip

    @m.Constraint(m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def repeat_weekday_recovery(blk, t, i):
        if t >= detla_time + 1 and t <= 4 * detla_time:  # Compare day 2-4 to day 1
            return (
                blk.period[1, t].reverse_osmosis.ro_skid[i].recovery
                == blk.period[1, t].reverse_osmosis.ro_skid[i].recovery
            )
        else:
            return Constraint.Skip

    @m.Constraint(m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def repeat_weekend_flowrate(blk, t, i):
        if t >= 5 * detla_time + 1 and t <= 6 * detla_time:  # Compare day 6 to day 7
            return (
                blk.period[1, t].reverse_osmosis.ro_skid[i].feed_flowrate
                == blk.period[1, t].reverse_osmosis.ro_skid[i].feed_flowrate
            )
        else:
            return Constraint.Skip

    @m.Constraint(m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def repeat_weekend_recovery(blk, t, i):
        if t >= 5 * detla_time + 1 and t <= 6 * detla_time:  # Compare day 6 to day 7
            return (
                blk.period[1, t].reverse_osmosis.ro_skid[i].recovery
                == blk.period[1, t].reverse_osmosis.ro_skid[i].recovery
            )
        else:
            return Constraint.Skip


def prevent_consecutive_flow_changes(m):
    # Prevent back-to-back flow-change events unless shutdown is occurring.
    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_ro.num_ro_skids + 1))
    def no_consecutive_ro_flow_changes(m_blk, d, t, i):
        if d == 1 and t <= 2:
            return Constraint.Skip

        shutdown_now = m_blk.period[d, t].reverse_osmosis.ro_skid[i].shutdown
        shutdown_prev = m_blk.period[d, t - 1].reverse_osmosis.ro_skid[i].shutdown

        return (
            m_blk.flow_changed[d, t, i] + m_blk.flow_changed[d, t - 1, i]
            <= 1 + shutdown_now + shutdown_prev
        )

    @m.Constraint(m.set_days, m.set_time, range(1, m.params.wrd_uf.num_uf_pumps + 1))
    def no_consecutive_uf_flow_changes(m_blk, d, t, i):
        if d == 1 and t <= 2:
            return Constraint.Skip

        shutdown_now = m_blk.period[d, t].pretreatment.uf_pumps[i].shutdown
        shutdown_prev = m_blk.period[d, t - 1].pretreatment.uf_pumps[i].shutdown

        return (
            m_blk.uf_flow_changed[d, t, i] + m_blk.uf_flow_changed[d, t - 1, i]
            <= 1 + shutdown_now + shutdown_prev
        )
