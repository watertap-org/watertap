#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of a chemical addition unit
operation.
"""

import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError

from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData

__author__ = "Andrew Lee, Kurban Sitterley"


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

        self.chemical_flow_mass = pyo.Var(
            initialize=1,
            units=pyo.units.kg / pyo.units.s,
            bounds=(0, None),
            doc="Mass flow rate of chemical solution",
        )

        self._fixed_perf_vars.append(self.chemical_dosage)
        self._fixed_perf_vars.append(self.solution_density)
        self._fixed_perf_vars.append(self.ratio_in_solution)

        self._perf_var_dict["Chemical Dosage"] = self.chemical_dosage
        self._perf_var_dict["Chemical Flow"] = self.chemical_flow_vol

        def rule_chem_flow_mass(blk):
            return blk.chemical_flow_mass == pyo.units.convert(
                blk.chemical_dosage
                * blk.properties[0].flow_vol
                / blk.ratio_in_solution,
                to_units=pyo.units.kg / pyo.units.s,
            )

        self.chemical_flow_mass_constraint = pyo.Constraint(rule=rule_chem_flow_mass)

        def rule_chem_flow_vol(blk):
            return blk.chemical_flow_vol[0] == pyo.units.convert(
                blk.chemical_flow_mass / blk.solution_density,
                to_units=pyo.units.m**3 / pyo.units.s,
            )

        self.chemical_flow_vol_constraint = pyo.Constraint(rule=rule_chem_flow_vol)

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

        if chem_name in ["alum"]:
            basis_units = pyo.units.gal / pyo.units.hr
            sizing_term = blk.unit_model.chemical_flow_vol[0] / basis_units
        elif chem_name in ["lime", "anhydrous_ammonia", "chlorine"]:
            basis_units = pyo.units.lb / pyo.units.day
            sizing_term = blk.unit_model.chemical_flow_mass / basis_units
        else:
            basis_units = pyo.units.gal / pyo.units.day
            sizing_term = blk.unit_model.chemical_flow_vol[0] / basis_units

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
            blk.unit_model.electricity[0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass, chem_name
        )
