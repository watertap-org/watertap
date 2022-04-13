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
