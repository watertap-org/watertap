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
This module contains a zero-order representation of a storage tank unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("StorageTankZO")
class StorageTankZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a storage tank unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        build_pt(self)
        constant_intensity(self)

        # self._has_recovery_removal = False

        self._tech_type = "storage_tank"

        self.storage_time = Var(
            self.flowsheet().time, units=pyunits.hours, doc="Storage time needed"
        )

        self.surge_capacity = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Additional capacity needed for surge flow",
        )

        self._fixed_perf_vars.append(self.storage_time)
        self._fixed_perf_vars.append(self.surge_capacity)

        self.tank_volume = Var(
            self.flowsheet().time, units=pyunits.m**3, doc="Storage tank volume"
        )

        @self.Constraint(self.flowsheet().time, doc="Tank volume constraint")
        def tank_volume_constraint(b, t):
            return b.tank_volume[t] == pyunits.convert(
                b.properties[t].flow_vol, to_units=pyunits.m**3 / pyunits.hr
            ) * b.storage_time[t] * (1 + b.surge_capacity[t])

        self._perf_var_dict["Storage Time (hr)"] = self.storage_time
        self._perf_var_dict["Surge Capacity (%)"] = self.surge_capacity
        self._perf_var_dict["Tank Volume (m3)"] = self.tank_volume

    @property
    def default_costing_method(self):
        return self.cost_storage_tank

    @staticmethod
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
        A, B = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        blk.unit_model._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
