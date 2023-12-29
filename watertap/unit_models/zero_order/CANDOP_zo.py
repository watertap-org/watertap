#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
This module contains a zero-order representation of a CANDO+P reactor unit.
"""

import pyomo.environ as pyo
from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData

# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("CANDOPZO")
class CANDOPData(ZeroOrderBaseData):
    """
    Zero-Order model for a CANDO+P reactor unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "CANDO_P"

        build_sido_reactive(self)

        # Create electricity variable
        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )
        self._perf_var_dict["Electricity Demand"] = self.electricity

        # Create electricity intensity variable and constraint. For this
        # model, electricity demand is calculated based on the amount of
        # nitrogen reacted.
        # TODO The information I have says that this electricity consumption
        #      accounts pumping, oxygenation, and stirring. At some point,
        #      perhaps we should come back and adjust the model to account for
        #      pumping costs separately.
        self.electricity_intensity_N = Var(
            units=pyunits.kWh / pyunits.kg,
            bounds=(0, None),
            doc="Electricity demand per kg N reacted",
        )
        self._fixed_perf_vars.append(self.electricity_intensity_N)
        self._perf_var_dict["Electricity Intensity"] = self.electricity_intensity_N

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "nitrogen consumed.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                pyunits.convert(
                    b.extent_of_reaction[t, "n_reaction"] * b.electricity_intensity_N,
                    to_units=pyunits.kW,
                )
            )

        # Create oxygen demand variables and constraint. The amount of oxygen
        # consumed is assumed to be a linear function of the amount of
        # nitrogen reacted.
        self.O2_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Oxygen demand",
        )
        self._perf_var_dict["Oxygen Demand"] = self.O2_demand
        self.oxygen_nitrogen_ratio = Var(
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Oxygen consumed - nitrogen reacted ratio",
        )
        self._fixed_perf_vars.append(self.oxygen_nitrogen_ratio)
        self._perf_var_dict[
            "Oxygen consumed / nitrogen reacted ratio (mass basis)"
        ] = self.oxygen_nitrogen_ratio

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for oxygen consumption."
        )
        def oxygen_consumption(b, t):
            return b.O2_demand[t] == (
                b.extent_of_reaction[t, "n_reaction"] * b.oxygen_nitrogen_ratio
            )

    @property
    def default_costing_method(self):
        return self.cost_CANDOP

    @staticmethod
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
        size_param, size_cost = blk.unit_model._get_tech_parameters(
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
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
