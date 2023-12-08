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
This module contains a zero-order representation of a UV reactor unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVZO")
class UVZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a UV unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "uv"

        build_siso(self)
        constant_intensity(self)

        self.uv_reduced_equivalent_dose = Var(
            self.flowsheet().time,
            units=pyunits.mJ / pyunits.cm**2,
            doc="Reduced equivalent dosage",
        )
        self.uv_transmittance_in = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="UV transmittance of solution at UV reactor inlet",
        )

        self.recovery_frac_mass_H2O.fix(1)
        self._fixed_perf_vars.append(self.uv_reduced_equivalent_dose)
        self._fixed_perf_vars.append(self.uv_transmittance_in)

        self._perf_var_dict[
            "UV Reduced Equivalent Dosage (mJ/cm^2)"
        ] = self.uv_reduced_equivalent_dose
        self._perf_var_dict["UV Transmittance of Feed"] = self.uv_transmittance_in

    @property
    def default_costing_method(self):
        return self.cost_uv

    @staticmethod
    def cost_uv(blk):
        """
        General method for costing UV reactor units. Capital cost is based on
        the inlet flow, UV reduced equivalent dosage, and UV transmittance at
        the inlet.
        """
        t0 = blk.flowsheet().time.first()
        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "reactor_cost",
                "lamp_cost",
            ],
        )

        expr = blk.unit_model._get_uv_capital_cost(blk, A, B)

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

    @staticmethod
    def _get_uv_capital_cost(blk, A, B):
        """
        Generate expression for capital cost of UV reactor.
        """
        t0 = blk.flowsheet().time.first()

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.m**3 / pyo.units.hr,
        )

        E = pyo.units.convert(blk.unit_model.electricity[t0], to_units=pyo.units.kW)

        expr = pyo.units.convert(
            A * Q + B * E,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        return expr
