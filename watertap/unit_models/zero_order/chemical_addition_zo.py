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
This module contains a zero-order representation of a chemical addition unit
operation.
"""
import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError

from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData

# Some more inforation about this module
__author__ = "Andrew Lee"


@declare_process_block_class("ChemicalAdditionZO")
class ChemicalAdditionZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a chemical addition unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "chemical_addition"

        if self.config.process_subtype is None:
            raise ConfigurationError(
                f"{self.name} - zero-order chemical addition operations "
                "require the process_subtype configuration argument to be set"
            )

        build_pt(self)

        self.chemical_dosage = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.mg / pyo.units.L,
            bounds=(0, None),
            doc="Dosing rate of chemical",
        )

        self.solution_density = pyo.Var(
            bounds=(0, None),
            units=pyo.units.kg / pyo.units.m**3,
            doc="Mass density of chemical solution",
        )
        self.ratio_in_solution = pyo.Var(
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Mass fraction of chemical in solution",
        )

        self.chemical_flow_vol = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.m**3 / pyo.units.s,
            bounds=(0, None),
            doc="Volumetric flow rate of chemical solution",
        )

        self._fixed_perf_vars.append(self.chemical_dosage)
        self._fixed_perf_vars.append(self.solution_density)
        self._fixed_perf_vars.append(self.ratio_in_solution)

        self._perf_var_dict["Chemical Dosage"] = self.chemical_dosage
        self._perf_var_dict["Chemical Flow"] = self.chemical_flow_vol

        def rule_chem_flow(blk, t):
            return blk.chemical_flow_vol[t] == pyo.units.convert(
                blk.chemical_dosage[t]
                * blk.properties[t].flow_vol
                / (blk.solution_density * blk.ratio_in_solution),
                to_units=pyo.units.m**3 / pyo.units.s,
            )

        self.chemical_flow_constraint = pyo.Constraint(
            self.flowsheet().time, rule=rule_chem_flow
        )

        pump_electricity(self, self.chemical_flow_vol)

    @property
    def default_costing_method(self):
        return self.cost_chemical_addition

    @staticmethod
    def cost_chemical_addition(blk, number_of_parallel_units=1):
        """
        General method for costing chemical addition processes. Capital cost is
        based on the mass flow rate of chemical added.
        This method also registers the chemical flow and electricity demand as
        costed flows.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        chem_name = blk.unit_model.config.process_subtype

        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (
            blk.unit_model.chemical_dosage[t0]
            * blk.unit_model.properties[t0].flow_vol
            / blk.unit_model.ratio_in_solution
        )
        sizing_term = blk.unit_model.chemical_flow_vol[t0] / (
            pyo.units.gal / pyo.units.day
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
        blk.config.flowsheet_costing_block.cost_flow(chem_flow_mass, chem_name)
