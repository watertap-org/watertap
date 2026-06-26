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


from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    RangeSet,
    Var,
    units as pyunits,
)

from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.apps.grid_integration import OperationModel

from watertap.flowsheets.flex_desal import params as um_params
from watertap.flowsheets.flex_desal.unit_models import _add_required_variables


def ro_skid_operation_model(blk, params: um_params.WRD_ROParams):
    """
    Builds operation model for an RO skid

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.coeffs = Param(["a", "b", "c"], initialize=params.surrogate_coeffs)

    blk.op_flow_limits_lower = Constraint(
        expr=blk.feed_flowrate >= blk.op_mode * params.minimum_flowrate,
        doc="Enforce minimum flowrate when operating",
    )
    blk.op_flow_limits_upper = Constraint(
        expr=blk.feed_flowrate <= blk.op_mode * params.maximum_flowrate,
        doc="Enforce maximum flowrate when operating",
    )

    if params.surrogate_type == "quadratic_energy_intensity":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (
                blk.coeffs["a"]
                + blk.coeffs["b"] * blk.feed_flowrate
                + blk.coeffs["c"] * blk.feed_flowrate**2
            ),
            doc="Calculates the specific energy requirement",
        )

    elif params.surrogate_type == "PySMO_polyfit":
        energy_surrogate = PysmoSurrogate.load_from_file(params.surrogate_file)
        if energy_surrogate._input_bounds["Feed Flow m3/hr"][0] != 0:
            raise ValueError(
                "Surrogate input bounds are not correct. Lower bound should be 0."
            )
        blk.energy_surrogate = SurrogateBlock()
        blk.energy_surrogate.build_model(
            energy_surrogate,
            input_vars=[blk.recovery, blk.feed_flowrate],  # RR,
            output_vars=[blk.energy_intensity],
        )

        # reject_flowrate = feed - product = feed*(1-recovery), already a linear
        # variable via mass_balance.  Using it here avoids the bilinear product
        # feed*(1-recovery)
        blk.flow_limit_from_RR = Constraint(
            expr=blk.reject_flowrate
            >= blk.op_mode * (3 * pyunits.m**3 / pyunits.hr * 15),
            doc="Minimum reject flowrate when operating (linearised via reject_flowrate variable)",
        )

    else:
        raise ValueError("Unrecognized surrogate type")


def wrd_reverse_osmosis_operation_model(blk, params: um_params.WRD_ROParams):
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
    # Defining a slack variable for flowrate that is not accounted for by the sum of RO intake pumps
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
    def symmetry_breaking_cuts_ro(b, index):
        if index == 1:
            return Constraint.Skip
        return b.ro_skid[index].op_mode <= b.ro_skid[index - 1].op_mode

    # Also add one for the flowrate itself
    @blk.Constraint(blk.set_ro_skids)
    def symmetry_breaking_cuts_ro_flowrate(b, index):
        if index == 1:
            return Constraint.Skip
        return b.ro_skid[index].feed_flowrate <= b.ro_skid[index - 1].feed_flowrate

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

    for skid in blk.set_ro_skids:
        # Bounds on recovery
        blk.ro_skid[skid].recovery.setlb(params.minimum_recovery)
        blk.ro_skid[skid].recovery.setub(params.maximum_recovery)
        # Add an upper bound for flow to help solver
        blk.ro_skid[skid].feed_flowrate.setub(params.maximum_flowrate)


def uf_pump_operation_model(blk, params: um_params.WRD_UFParams):
    """
    Builds operation model for a UF pump

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.coeffs = Param(["a", "b", "c"], initialize=params.surrogate_coeffs)

    blk.operational_limits_lower = Constraint(
        expr=blk.feed_flowrate >= blk.op_mode * params.minimum_flowrate,
        doc="Enforce minimum flowrate when operating",
    )
    blk.operational_limits_upper = Constraint(
        expr=blk.feed_flowrate <= blk.op_mode * params.maximum_flowrate,
        doc="Enforce maximum flowrate when operating",
    )

    if params.surrogate_type == "linear_energy_intensity":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (blk.coeffs["a"] + blk.coeffs["b"] * blk.feed_flowrate),
            doc="Calculates the specific energy requirement",
        )  # This shouldn't be needed tbh. Linear is quadratic with c=0
    elif params.surrogate_type == "quadratic_energy_intensity":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (
                blk.coeffs["a"]
                + blk.coeffs["b"] * blk.feed_flowrate
                + blk.coeffs["c"] * blk.feed_flowrate**2
            ),
            doc="Calculates the specific energy requirement",
        )
    else:
        raise ValueError("Unrecognized surrogate type")


def wrd_uf_operation_model(blk, params: um_params.WRD_UFParams):
    """
    Builds operation model for UF unit in WRD case

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
    blk.set_uf_pumps = RangeSet(params.num_uf_pumps)
    blk.uf_pumps = OperationModel(
        blk.set_uf_pumps,
        model_func=uf_pump_operation_model,
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
        == sum(blk.uf_pumps[i].feed_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance at the feed",
    )
    blk.product_mass_balance = Constraint(
        expr=blk.product_flowrate
        == sum(blk.uf_pumps[i].product_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance on permeate side",
    )
    blk.reject_mass_balance = Constraint(
        expr=blk.reject_flowrate
        == sum(blk.uf_pumps[i].reject_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance on brine side",
    )
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption
        == sum(blk.uf_pumps[i].power_consumption for i in blk.set_uf_pumps),
        doc="Calculates the total power requirement for RO",
    )

    # symmetry breaking for >1 skid. Skids can only operate if the previous skid is on
    @blk.Constraint(blk.set_uf_pumps)
    def symmetry_breaking_cuts_uf(b, index):
        if index == 1:
            return Constraint.Skip
        return b.uf_pumps[index].op_mode <= b.uf_pumps[index - 1].op_mode

    # Also add one for the flowrate itself
    @blk.Constraint(blk.set_uf_pumps)
    def symmetry_breaking_cuts_uf_flowrate(b, index):
        if index == 1:
            return Constraint.Skip
        return b.uf_pumps[index].feed_flowrate <= b.uf_pumps[index - 1].feed_flowrate

    # Ensure that the operation of minimum number of skids is identical
    blk.set_min_operating_pumps = RangeSet(2, params.minimum_operating_pumps)

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_startup(b, index):
        return b.uf_pumps[index].startup == b.uf_pumps[1].startup

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_op_mode(b, index):
        return b.uf_pumps[index].op_mode == b.uf_pumps[1].op_mode

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_shutdown(b, index):
        return b.uf_pumps[index].shutdown == b.uf_pumps[1].shutdown

    # Set a maximum flowrate to reduce search space
    for pump in blk.set_uf_pumps:
        blk.uf_pumps[pump].feed_flowrate.setub(params.maximum_flowrate)
