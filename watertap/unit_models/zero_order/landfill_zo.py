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
This module contains a zero-order representation of a landfill unit.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

# Some more inforation about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("LandfillZO")
class LandfillZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a landfill unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "landfill"

        build_pt(self)
        constant_intensity(self)

        self.capacity_basis = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            doc="capacity basis for capital cost",
        )

        self.total_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            doc="total mass flow rate",
        )

        self._fixed_perf_vars.append(self.capacity_basis)

        @self.Constraint(self.flowsheet().time, doc="Total mass constraint")
        def total_mass_constraint(b, t):
            return b.total_mass[t] == sum(
                pyunits.convert(
                    b.inlet.flow_mass_comp[t, m], to_units=pyunits.kg / pyunits.hr
                )
                for m in b.config.property_package.component_list
            )

        self._perf_var_dict["Capacity Basis (kg/hr)"] = self.capacity_basis
        self._perf_var_dict["Total Mass (kg/hr)"] = self.total_mass

    @property
    def default_costing_method(self):
        return self.cost_landfill

    @staticmethod
    def cost_landfill(blk, number_of_parallel_units=1):
        """
        General method for costing landfill. Capital cost is based on the total mass and
        capacity basis.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """

        t0 = blk.flowsheet().time.first()
        sizing_term = blk.unit_model.total_mass[t0] / blk.unit_model.capacity_basis[t0]

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
