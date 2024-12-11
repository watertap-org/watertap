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
This module contains a zero-order representation of a nanofiltration unit
operation.
"""

import pyomo.environ as pyo
from idaes.core import declare_process_block_class
from pyomo.environ import Var, units as pyunits
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Andrew Lee, Adam Atia"


@declare_process_block_class("NanofiltrationZO")
class NanofiltrationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Nanofiltration unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "nanofiltration"

        build_sido(self)

        if (
            self.config.process_subtype == "default"
            or self.config.process_subtype is None
        ):
            constant_intensity(self)
        else:
            self.rejection_comp = Var(
                self.flowsheet().time,
                self.config.property_package.config.solute_list,
                units=pyunits.dimensionless,
                doc="Component rejection",
            )

            self.water_permeability_coefficient = Var(
                self.flowsheet().time,
                units=pyunits.L / pyunits.m**2 / pyunits.hour / pyunits.bar,
                doc="Membrane water permeability coefficient, A",
            )

            self.applied_pressure = Var(
                self.flowsheet().time,
                units=pyunits.bar,
                doc="Net driving pressure across membrane",
            )

            self.area = Var(units=pyunits.m**2, doc="Membrane area")

            self._fixed_perf_vars.append(self.applied_pressure)
            self._fixed_perf_vars.append(self.water_permeability_coefficient)

            @self.Constraint(self.flowsheet().time, doc="Water permeance constraint")
            def water_permeance_constraint(b, t):
                return b.properties_treated[t].flow_vol == pyunits.convert(
                    b.water_permeability_coefficient[t]
                    * b.area
                    * b.applied_pressure[t],
                    to_units=pyunits.m**3 / pyunits.s,
                )

            @self.Constraint(
                self.flowsheet().time,
                self.config.property_package.config.solute_list,
                doc="Solute [observed] rejection constraint",
            )
            def rejection_constraint(b, t, j):
                return (
                    b.rejection_comp[t, j]
                    == 1
                    - b.properties_treated[t].conc_mass_comp[j]
                    / b.properties_in[t].conc_mass_comp[j]
                )

            self._perf_var_dict["Membrane Area (m^2)"] = self.area
            self._perf_var_dict["Net Driving Pressure (bar)"] = self.applied_pressure
            self._perf_var_dict["Water Permeability Coefficient (LMH/bar)"] = (
                self.water_permeability_coefficient
            )
            self._perf_var_dict[f"Rejection"] = self.rejection_comp

    @property
    def default_costing_method(self):
        return self.cost_nanofiltration

    @staticmethod
    def cost_nanofiltration(blk, number_of_parallel_units=1):
        """
        General method for costing nanofiltration. Costing is carried out
        using either the general_power_law form or the standard form which
        computes membrane cost and replacement rate.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        # Get cost method for this technology
        cost_method = blk.unit_model._get_unit_cost_method(blk)
        valid_methods = ["cost_power_law_flow", "cost_membrane"]
        if cost_method == "cost_power_law_flow":
            blk.unit_model.cost_power_law_flow(blk, number_of_parallel_units)
        elif cost_method == "cost_membrane":
            # NOTE: number of units does not matter for cost_membrane
            #       as its a linear function of membrane area
            blk.unit_model.cost_membrane(blk)
        else:
            raise KeyError(
                f"{cost_method} is not a relevant cost method for "
                f"{blk.unit_model._tech_type}. Specify one of the following "
                f"cost methods in the unit's YAML file: {valid_methods}"
            )

    @staticmethod
    def cost_membrane(blk):
        """
        Get membrane cost based on membrane area and unit membrane costs
        as well as fixed operating cost for membrane replacement.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )
        # Get costing parameter sub-block for this technology
        mem_cost, rep_rate = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["membrane_cost", "membrane_replacement_rate"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit operation",
        )

        capex_expr = pyo.units.convert(
            mem_cost * pyo.units.convert(blk.unit_model.area, to_units=pyo.units.m**2),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * capex_expr
        )

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                rep_rate
                * mem_cost
                * pyo.units.convert(blk.unit_model.area, to_units=pyo.units.m**2),
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )
