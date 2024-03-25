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
This module contains a zero-order representation of a Ozone-AOP unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.unit_models.zero_order.ozone_zo import OzoneZOData
from watertap.unit_models.zero_order.aop_addition_zo import AOPAdditionMixin

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("OzoneAOPZO")
class OzoneAOPZOData(OzoneZOData, AOPAdditionMixin):
    """
    Zero-Order model for a Ozone-AOP unit operation.
    """

    def build(self):
        super().build()

        self._tech_type = "ozone_aop"

        self.oxidant_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Oxidant dosage"
        )

        self.chemical_flow_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Mass flow rate of oxidant solution",
        )

        self.ozone_toc_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of ozone to total organic carbon",
        )

        self.oxidant_ozone_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of oxidant to ozone",
        )

        self._fixed_perf_vars.append(self.oxidant_ozone_ratio)

        @self.Constraint(self.flowsheet().time, doc="Ozone/TOC ratio constraint")
        def ozone_toc_ratio_constraint(b, t):
            return b.ozone_toc_ratio[t] == 1 + pyunits.convert(
                b.concentration_time[t]
                / b.contact_time[t]
                / b.properties_in[t].conc_mass_comp["toc"],
                to_units=pyunits.dimensionless,
            )

        @self.Constraint(self.flowsheet().time, doc="Oxidant dose constraint")
        def oxidant_dose_constraint(b, t):
            return b.oxidant_dose[t] == pyunits.convert(
                b.oxidant_ozone_ratio[t]
                * b.ozone_toc_ratio[t]
                * b.properties_in[t].conc_mass_comp["toc"],
                to_units=pyunits.mg / pyunits.L,
            )

        @self.Constraint(self.flowsheet().time, doc="Oxidant mass flow constraint")
        def chemical_flow_mass_constraint(b, t):
            return b.chemical_flow_mass[t] == pyunits.convert(
                b.oxidant_dose[t] * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        self._perf_var_dict["Oxidant Dosage (mg/L)"] = self.oxidant_dose
        self._perf_var_dict["Oxidant Flow (kg/s)"] = self.chemical_flow_mass
        self._perf_var_dict["Oxidant/Ozone Ratio"] = self.oxidant_ozone_ratio
        self._perf_var_dict["Ozone/TOC Ratio"] = self.ozone_toc_ratio

    @property
    def default_costing_method(self):
        return self.cost_ozonation_aop

    @staticmethod
    def cost_ozonation_aop(blk):
        """
        General method for costing ozonation with AOP. Capital cost is
        based on the inlet flowrate, dosage of ozone and flow rate of H2O2.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "ozone_capital_a_parameter",
                "ozone_capital_b_parameter",
                "ozone_capital_c_parameter",
                "ozone_capital_d_parameter",
                "aop_capital_a_parameter",
                "aop_capital_b_parameter",
            ],
        )

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        # Get costing term for ozone addition
        expr = blk.unit_model._get_ozone_capital_cost(blk, A, B, C, D)

        # Add costing term for AOP addition
        expr += blk.unit_model._get_aop_capital_cost(blk, E, F)

        blk.capital_cost_constraint = pyo.Constraint(expr=blk.capital_cost == expr)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.chemical_flow_mass[t0], "hydrogen_peroxide"
        )
