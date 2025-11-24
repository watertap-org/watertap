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
This module contains unit models needed for flexible desalination
analysis.
"""

from idaes.apps.grid_integration import OperationModel
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    RangeSet,
    Var,
    exp,
    units as pyunits,
)
from watertap.flowsheets.flex_desal import params as um_params


# NOTE: OperationModel class automatically adds startup, shutdown,
# and op_mode binary variables. So, no need to define these variables
# explicitly.
def _add_required_variables(blk):
    """Function for defining common variables and constraints"""
    # Declare variables
    blk.energy_intensity = Var(
        within=NonNegativeReals, units=pyunits.kWh / pyunits.m**3
    )
    blk.power_consumption = Var(within=NonNegativeReals, units=pyunits.kW)

    blk.recovery = Var(
        within=NonNegativeReals, bounds=(0, 1), units=pyunits.dimensionless
    )
    blk.feed_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    blk.product_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    blk.reject_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)

    # Declare model constraints
    blk.mass_balance = Constraint(
        expr=blk.feed_flowrate == blk.product_flowrate + blk.reject_flowrate,
        doc="Overall mass balance",
    )
    blk.calculate_product_flowrate = Constraint(
        expr=blk.product_flowrate == blk.feed_flowrate * blk.recovery,
        doc="Computes product flowrate",
    )
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption == blk.energy_intensity * blk.product_flowrate,
        doc="Power requirement for the unit",
    )


def intake_operation_model(blk, params: um_params.IntakeParams):
    """
    Builds operation model for the intake

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.recovery.fix(params.get_recovery)
    blk.energy_intensity.fix(params.energy_intensity)


def pretreatment_operation_model(blk, params: um_params.PretreatmentParams):
    """
    Builds operation model for the pretreatment unit

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.recovery.fix(params.get_recovery)
    blk.energy_intensity.fix(params.energy_intensity)


def ro_skid_operation_model(blk, params: um_params.ROParams):
    """
    Builds operation model for an RO skid

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    # # Remove the constant if its magnitude is too small
    # if abs(params.surrogate_d) < 1e-3:
    #     params.surrogate_d = 0

    _add_required_variables(blk)
    blk.coeffs = Param(["a", "b", "c", "d"], initialize=params.surrogate_coeffs)

    blk.operational_limits = Constraint(
        expr=blk.feed_flowrate == blk.op_mode * params.nominal_flowrate
    )

    if params.surrogate_type == "exponential_quadratic":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (
                blk.coeffs["a"] * exp(-blk.coeffs["b"] * blk.recovery)
                + blk.coeffs["c"] * blk.recovery**2
                + blk.coeffs["d"]
            ),
            doc="Calculates the specific energy requirement",
        )

    elif params.surrogate_type == "quadratic_surrogate":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (
                blk.coeffs["a"] * blk.recovery**2
                + blk.coeffs["b"] * blk.recovery
                + blk.coeffs["c"]
            )
        )


def reverse_osmosis_operation_model(blk, params: um_params.ROParams):
    """
    Builds operation model for the reverse osmosis unit

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    # Declare required variables
    _add_required_variables(blk)
    blk.inlet_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    # Defining a slack variable for flowrate that is not accounted
    # for by the sum of RO intake pumps
    blk.leftover_flow = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)

    # Build RO skid models
    blk.set_ro_skids = RangeSet(params.num_ro_skids)
    blk.ro_skid = OperationModel(
        blk.set_ro_skids,
        model_func=ro_skid_operation_model,
        model_args={"params": params},
        minimum_up_time=params.minimum_uptime,
        minimum_down_time=params.minimum_downtime,
    )

    # Remove overall mass balance and power consumption calculation
    blk.del_component(blk.recovery)
    blk.del_component(blk.energy_intensity)
    blk.del_component(blk.mass_balance)
    blk.del_component(blk.calculate_product_flowrate)
    blk.del_component(blk.calculate_power_consumption)

    # Declare required constraints
    blk.calculate_leftover_flow = Constraint(
        expr=blk.feed_flowrate == blk.inlet_flowrate + blk.leftover_flow,
        doc="Calculates leftover flowrate",
    )
    blk.feed_mass_balance = Constraint(
        expr=blk.inlet_flowrate
        == sum(blk.ro_skid[i].feed_flowrate for i in blk.set_ro_skids),
        doc="Mass balance at the feed",
    )
    blk.product_mass_balance = Constraint(
        expr=blk.product_flowrate
        == sum(blk.ro_skid[i].product_flowrate for i in blk.set_ro_skids),
        doc="Mass balance on permeate side",
    )
    blk.reject_mass_balance = Constraint(
        expr=blk.reject_flowrate
        == sum(blk.ro_skid[i].reject_flowrate for i in blk.set_ro_skids),
        doc="Mass balance on brine side",
    )
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption
        == sum(blk.ro_skid[i].power_consumption for i in blk.set_ro_skids),
        doc="Calculates the total power requirement for RO",
    )

    # symmetry breaking for >1 skid. Skids can only operate if the previous skid is on
    @blk.Constraint(blk.set_ro_skids)
    def symmetry_breaking_cuts(b, index):
        if index == 1:
            return Constraint.Skip
        return b.ro_skid[index].op_mode <= b.ro_skid[index - 1].op_mode

    # Ensure that the operation of minimum number of skids is identical
    blk.set_min_operating_skids = RangeSet(2, params.minimum_operating_skids)

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_startup(b, index):
        return b.ro_skid[index].startup == b.ro_skid[1].startup

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_op_mode(b, index):
        return b.ro_skid[index].op_mode == b.ro_skid[1].op_mode

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_shutdown(b, index):
        return b.ro_skid[index].shutdown == b.ro_skid[1].shutdown

    # Update bounds on recovery and energy intensity for all skids
    ei_lb, ei_ub = params.get_energy_intensity_bounds()
    for skid in blk.set_ro_skids:
        blk.ro_skid[skid].recovery.setlb(params.minimum_recovery)
        blk.ro_skid[skid].recovery.setub(params.maximum_recovery)
        blk.ro_skid[skid].energy_intensity.setlb(ei_lb)
        blk.ro_skid[skid].energy_intensity.setub(ei_ub)


def posttreatment_operation_model(blk, params: um_params.FlexDesalParams):
    """
    Builds operation model for the posttreatment unit

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)

    blk.energy_intensity.fix(params.posttreatment.energy_intensity)
    blk.del_component(blk.recovery)
    blk.del_component(blk.calculate_product_flowrate)

    # If the posttreatment unit is not operating, then set product flowrate to zero
    # Connect post-treatment operation with the RO startup variables
    blk.suppress_product_flowrate = Constraint(
        expr=blk.product_flowrate <= params.intake.nominal_flowrate * blk.op_mode
    )
    blk.suppress_reject_flowrate = Constraint(
        expr=blk.reject_flowrate <= params.intake.nominal_flowrate * (1 - blk.op_mode)
    )

    # Update the power consumption calculation constraint
    blk.calculate_power_consumption.set_value(
        blk.power_consumption == blk.energy_intensity * blk.feed_flowrate,
    )


def brine_discharge_operation_model(blk, params: um_params.FlexDesalParams):
    """
    Builds operation model for the brine discharge unit

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    # Declare model parameters
    blk.energy_intensity = Param(
        initialize=params.brinedischarge.energy_intensity,
        units=pyunits.kWh / pyunits.m**3,
        mutable=True,
    )

    # Declare essential variables
    blk.feed_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    blk.power_consumption = Var(within=NonNegativeReals, units=pyunits.kW)

    # the brine sump only consumes power if the RO is off,
    # otherwise brine is pushed out by the leftover RO pressure
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption
        >= blk.energy_intensity
        * (blk.feed_flowrate + params.intake.nominal_flowrate * (blk.op_mode - 1)),
        doc="Power requirement for brine discharge",
    )


def power_generation_operation_model(blk, params: um_params.FlexDesalParams):
    """
    Builds the operation model for onsite power generation

    Parameters
    ----------
    blk : OperationModel
        Instance of IDAES OperationModel

    design_blk : DesignModel
        Design model containing information on the peak capacity
    """
    # Declare model parameters
    blk.capacity_factor = Param(
        initialize=0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="capacity factor of onsite power generation",
    )

    # Declare essential variables
    # The installed capacity could exceed the power requirement
    # We will assume that the excess power is curtailed.
    blk.power_utilized = Var(
        initialize=0,
        within=NonNegativeReals,
        units=pyunits.kW,
        doc="Power utilized by the system",
    )
    blk.power_curtailed = Var(
        initialize=0,
        within=NonNegativeReals,
        units=pyunits.kW,
        doc="Total power curtailed by the system",
    )

    # Energy balance: Sum of power utilized and the power
    # curtailed must be equal to the total power generated.
    blk.calculate_power_generation = Constraint(
        expr=(
            blk.power_utilized + blk.power_curtailed
            == params.onsite_capacity * blk.capacity_factor
        ),
        doc="Computes the total power generated onsite",
    )
