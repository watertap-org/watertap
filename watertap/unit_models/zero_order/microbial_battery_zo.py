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
This module contains a zero-order representation of a microbial battery water
treatment unit.
"""

import pyomo.environ as pyo
from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData


# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("MicrobialBatteryZO")
class MicrobialBatteryData(ZeroOrderBaseData):
    """
    Zero-Order model for a microbial battery water treatment unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "microbial_battery"

        build_sido_reactive(self)
        constant_intensity(self)

        # Create hydraulic retention time variable
        self.HRT = Var(
            units=pyunits.hr,
            bounds=(0, None),
            doc="Hydraulic retention time of water treatment unit",
        )
        self._perf_var_dict["Hydraulic Retention Time"] = self.HRT
        self._fixed_perf_vars.append(self.HRT)

        # Create reactor volume variable
        self.reactor_volume = Var(
            units=pyunits.m**3,
            bounds=(0, None),
            doc="Volume of water treatment unit",
        )
        self._perf_var_dict["Reactor Volume"] = self.reactor_volume

        @self.Constraint(self.flowsheet().time, doc="Constraint for reactor volume.")
        def reactor_volume_rule(b, t):
            return b.reactor_volume == (
                pyunits.convert(
                    b.HRT * b.properties_in[t].flow_vol, to_units=pyunits.m**3
                )
            )

    @property
    def default_costing_method(self):
        return self.cost_microbial_battery

    @staticmethod
    def cost_microbial_battery(blk):
        """
        General method for costing microbial battery treatment unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        sizing_cost = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * sizing_cost,
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
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_in[t0].flow_mass_comp["filtration_media"],
            "filtration_media",
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp["filtration_media"],
            "filtration_media_disposal",
        )
