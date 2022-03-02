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
from idaes.generic_models.costing.costing_base import FlowsheetCostingBlockData

from watertap.core.zero_order_base import ZeroOrderBase


@declare_process_block_class("ZeroOrderCosting")
class ZeroOrderCostingData(FlowsheetCostingBlockData):

    # Register currency and conversion rates based on CE Index
    # TODO : Consider way to do this from data
    pyo.units.load_definitions_from_strings([
        "USD_500 = [currency]",
        "USD_1990 = 500/357.6 * USD_500",
        "USD_1991 = 500/361.3 * USD_500",
        "USD_1992 = 500/358.2 * USD_500",
        "USD_1993 = 500/359.2 * USD_500",
        "USD_1994 = 500/368.1 * USD_500",
        "USD_1995 = 500/381.1 * USD_500",
        "USD_1996 = 500/381.7 * USD_500",
        "USD_1997 = 500/386.5 * USD_500",
        "USD_1998 = 500/389.5 * USD_500",
        "USD_1999 = 500/390.6 * USD_500",
        "USD_2000 = 500/394.1 * USD_500",
        "USD_2001 = 500/394.3 * USD_500",
        "USD_2002 = 500/395.6 * USD_500",
        "USD_2003 = 500/402.0 * USD_500",
        "USD_2004 = 500/444.2 * USD_500",
        "USD_2005 = 500/468.2 * USD_500",
        "USD_2006 = 500/499.6 * USD_500",
        "USD_2007 = 500/525.4 * USD_500",
        "USD_2008 = 500/575.4 * USD_500",
        "USD_2009 = 500/521.9 * USD_500",
        "USD_2010 = 500/550.8 * USD_500",
        "USD_2011 = 500/585.7 * USD_500",
        "USD_2012 = 500/584.6 * USD_500",
        "USD_2013 = 500/567.3 * USD_500",
        "USD_2014 = 500/576.1 * USD_500",
        "USD_2015 = 500/556.8 * USD_500",
        "USD_2016 = 500/541.7 * USD_500",
        "USD_2017 = 500/567.5 * USD_500",
        "USD_2018 = 500/671.1 * USD_500",
        "USD_2019 = 500/680.0 * USD_500"])

    # TODO: Define any default flow costs of interest

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Define expected flows
        self.defined_flows = {
            "electricity": 0.0595*pyo.units.USD_2019/pyo.units.kW/pyo.units.hour}

    def build_process_costs(self):
        """
        This is where you do all your process wide costing.
        This is completely up to you, but you will have access to the
        following aggregate costs:

            1. blk.aggregate_capital_cost
            2. blk.aggregate_fixed_operating_cost
            3. blk.aggregate_variable_operating_cost
            4. blk.aggregate_flow_costs (indexed by flow type)
        """
        # TODO: Do we have any process level methods to add here?
        pass

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
    def exponential_form(blk):
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

        try:
            sblock = blk.unit_model.properties_in[t0]
        except AttributeError:
            # Pass-through case
            sblock = blk.unit_model.properties[t0]

        # TODO: More bases
        if basis == "flow_vol_inlet":
            state = sblock.flow_vol
        else:
            raise ValueError(
                f"{blk.name} - unrecognized basis in parameter declaration: "
                f"{basis}.")

        # Get reference flow and parameters
        try:
            pblock = getattr(blk.config.flowsheet_costing_block,
                             blk.unit_model._tech_type)
        except AttributeError:
            # Parameters for this tehcnology haven't been added yet
            pblock = pyo.Block()

            # Add block to FlowsheetCostingBlock
            blk.config.flowsheet_costing_block.add_component(
                blk.unit_model._tech_type, pblock)

            # Add requried parameters
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
            pblock.reference_state = pyo.Var(
                initialize=float(
                    parameter_dict[
                        "capital_cost"]["reference_state"]["value"]),
                units=getattr(
                    pyo.units,
                    parameter_dict[
                        "capital_cost"]["reference_state"]["units"]),
                doc="Reference state for capital cost expression")

            pblock.capital_a_parameter.fix()
            pblock.capital_b_parameter.fix()
            pblock.reference_state.fix()

        A = pblock.capital_a_parameter
        B = pblock.capital_b_parameter
        state_ref = pblock.reference_state

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost ==
            A*pyo.units.convert(state/state_ref,
                                to_units=pyo.units.dimensionless)**B)

        # Register flows if present
        if hasattr(blk.unit_model, "electricity"):
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.electricity[t0], "electricity")

    # -------------------------------------------------------------------------
    # Map costing methods to unit model classes
    unit_mapping = {ZeroOrderBase: exponential_form}
