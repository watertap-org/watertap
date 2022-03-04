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
import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.generic_models.costing.costing_base import (
    FlowsheetCostingBlockData, register_idaes_currency_units)

from watertap.core.zero_order_base import ZeroOrderBase


@declare_process_block_class("ZeroOrderCosting")
class ZeroOrderCostingData(FlowsheetCostingBlockData):

    # Register currency and conversion rates based on CE Index
    # TODO : Consider way to do this from data
    register_idaes_currency_units()

    def build_global_params(self):
        """
        To minimize overhead, only create global parameters for now.

        Unit-specific parameters will be added as sub-Blocks on a case-by-case
        basis as a unit of that type is costed.
        """
        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Define expected flows
        self.defined_flows = {
            "electricity":
                0.0595*pyo.units.USD_2019/pyo.units.kW/pyo.units.hour}

        # Costing factors
        self.land_cost_percent_FCI = pyo.Var(initialize=0.07,
                                             units=pyo.units.dimensionless,
                                             doc="Land cost as % FCI")
        self.working_capital_percent_FCI = pyo.Var(
            initialize=0.008,
            units=pyo.units.dimensionless,
            doc="Working capital as % FCI")

        # TODO: Load values from database
        # Fix all Vars
        for v in self.component_objects(pyo.Var, descend_into=True):
            v.fix()

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

        # Other variable costs

    def initialize(self):
        """
        Here we can add intialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now, no additional process level costs to initialize
        pass

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
            # If not present, call emthod to create parameter Block
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
