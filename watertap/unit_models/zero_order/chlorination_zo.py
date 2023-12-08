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
This module contains a zero-order representation of a Chlorination unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("ChlorinationZO")
class ChlorinationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Chlorination unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "chlorination"

        build_siso(self)
        constant_intensity(self)

        self.initial_chlorine_demand = Var(
            self.flowsheet().time,
            units=pyunits.mg / pyunits.liter,
            doc="Initial chlorine demand",
        )

        self.contact_time = Var(
            self.flowsheet().time, units=pyunits.hour, doc="Chlorine contact time"
        )
        self.concentration_time = Var(
            self.flowsheet().time,
            units=(pyunits.mg * pyunits.minute) / pyunits.liter,
            doc="CT value for chlorination",
        )
        self.chlorine_decay_rate = Var(
            self.flowsheet().time,
            units=pyunits.mg / (pyunits.L * pyunits.hour),
            doc="Chlorine decay rate",
        )

        self.recovery_frac_mass_H2O.fix(1)
        self._fixed_perf_vars.append(self.initial_chlorine_demand)
        self._fixed_perf_vars.append(self.contact_time)
        self._fixed_perf_vars.append(self.concentration_time)
        self._fixed_perf_vars.append(self.chlorine_decay_rate)

        self.chlorine_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Chlorine dose"
        )

        @self.Constraint(self.flowsheet().time, doc="Chlorine dose constraint")
        def chlorine_dose_constraint(b, t):
            return b.chlorine_dose[t] == self.initial_chlorine_demand[
                t
            ] + self.chlorine_decay_rate[t] * self.contact_time[t] + (
                self.concentration_time[t]
                / pyunits.convert(self.contact_time[t], to_units=pyunits.minute)
            )

        self._perf_var_dict["Chlorine Dose (mg/L)"] = self.chlorine_dose
        self._perf_var_dict[
            "Initial Chlorine Demand (mg/L)"
        ] = self.initial_chlorine_demand
        self._perf_var_dict["Contact Time (hr)"] = self.contact_time
        self._perf_var_dict["CT Value ((mg*min)/L)"] = self.concentration_time
        self._perf_var_dict[
            "Chlorine Decay Rate (mg/(L*hr))"
        ] = self.chlorine_decay_rate

    @property
    def default_costing_method(self):
        return self.cost_chlorination

    @staticmethod
    def cost_chlorination(blk):
        """
        General method for costing chlorination units. Capital cost is based on
        the both inlet flow and dosage of chlorine.
        This method also registers the chemical flow and electricity demand as
        costed flows.
        """
        t0 = blk.flowsheet().time.first()
        chem_flow_mass = (
            blk.unit_model.chlorine_dose[t0] * blk.unit_model.properties_in[t0].flow_vol
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "capital_c_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        ln_Q = pyo.log(
            pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol
                / (pyo.units.m**3 / pyo.units.hour),
                to_units=pyo.units.dimensionless,
            )
        )
        ln_D = pyo.log(
            pyo.units.convert(
                blk.unit_model.chlorine_dose[t0] / (pyo.units.mg / pyo.units.liter),
                to_units=pyo.units.dimensionless,
            )
        )

        expr = pyo.units.convert(
            A * ln_Q + B * ln_D + C * ln_Q * ln_D,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

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
        blk.config.flowsheet_costing_block.cost_flow(chem_flow_mass, "chlorine")
