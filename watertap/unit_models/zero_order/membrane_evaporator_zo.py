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
This module contains a zero-order representation of a membrane evaporation.
"""

import pyomo.environ as pyo
from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData


# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("MembraneEvaporatorZO")
class MembraneEvaporatorData(ZeroOrderBaseData):
    """
    Zero-Order model for a membrane evaporator.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "membrane_evaporator"

        build_sido(self)
        constant_intensity(self)

        # Create water flux variable
        self.water_flux = Var(
            units=pyunits.m / pyunits.hr,
            bounds=(0, None),
            doc="Water flux through membrane",
        )
        self._perf_var_dict["Water Flux"] = self.water_flux
        self._fixed_perf_vars.append(self.water_flux)

        # Create membrane area variable
        self.membrane_area = Var(
            units=pyunits.m**2, bounds=(0, None), doc="Membrane area"
        )
        self._perf_var_dict["Membrane Area"] = self.membrane_area

        @self.Constraint(self.flowsheet().time, doc="Constraint for water flux.")
        def wat_flux(b, t):
            return b.properties_byproduct[t].flow_vol == (
                pyunits.convert(
                    b.water_flux * b.membrane_area, to_units=pyunits.m**3 / pyunits.s
                )
            )

    @property
    def default_costing_method(self):
        return self.cost_membrane_evaporator

    @staticmethod
    def cost_membrane_evaporator(blk):
        """
        General method for costing membrane evaporator unit.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        memb_cost = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["membrane_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.membrane_area * memb_cost,
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
        # TODO: Consider adding heat as a registered flow, since the inlet
        #       stream to the membrane evaporator must be heated to ~37 C.
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
