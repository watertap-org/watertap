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
This module contains a zero-order representation of a deep well injection unit.
"""

import pyomo.environ as pyo
from pyomo.environ import Reference, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("DeepWellInjectionZO")
class DeepWellInjectionZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a deep well injection unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "deep_well_injection"
        build_pt(self)
        self._Q = Reference(self.properties[:].flow_vol)

        pump_electricity(self, self._Q)

        self.pipe_distance = Var(
            self.flowsheet().config.time, units=pyunits.miles, doc="Piping distance"
        )

        self.pipe_diameter = Var(
            self.flowsheet().config.time, units=pyunits.inches, doc="Pipe diameter"
        )

        self.flow_basis = Var(
            self.flowsheet().time, units=pyunits.m**3 / pyunits.hour, doc="flow basis"
        )

        self._fixed_perf_vars.append(self.pipe_distance)
        self._fixed_perf_vars.append(self.pipe_diameter)
        self._fixed_perf_vars.append(self.flow_basis)

        self._perf_var_dict["Pipe Distance (miles)"] = self.pipe_distance
        self._perf_var_dict["Pipe Diameter (inches)"] = self.pipe_diameter

    @property
    def default_costing_method(self):
        return self.cost_deep_well_injection

    @staticmethod
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
        A, B, C = blk.unit_model._get_tech_parameters(
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
        blk.unit_model._general_power_law_form(
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
