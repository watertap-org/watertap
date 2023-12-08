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
This module contains a zero-order representation of a Ozone reactor unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from watertap.core import build_siso, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("OzoneZO")
class OzoneZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Ozone unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "ozonation"

        build_siso(self)

        if "toc" not in self.config.property_package.config.solute_list:
            raise ConfigurationError(
                "toc must be in solute list for Ozonation or Ozone/AOP"
            )

        self.contact_time = Var(
            self.flowsheet().time, units=pyunits.minute, doc="Ozone contact time"
        )

        self.concentration_time = Var(
            self.flowsheet().time,
            units=(pyunits.mg * pyunits.minute) / pyunits.liter,
            doc="CT value for ozone contactor",
        )

        self.mass_transfer_efficiency = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ozone mass transfer efficiency",
        )

        self.specific_energy_coeff = Var(
            self.flowsheet().time,
            units=pyunits.kWh / pyunits.lb,
            bounds=(0, None),
            doc="Specific energy consumption for ozone generation",
        )

        self._fixed_perf_vars.append(self.contact_time)
        self._fixed_perf_vars.append(self.concentration_time)
        self._fixed_perf_vars.append(self.mass_transfer_efficiency)
        self._fixed_perf_vars.append(self.specific_energy_coeff)

        self.ozone_flow_mass = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.lb / pyunits.hr,
            doc="Mass flow rate of ozone",
        )

        self.ozone_consumption = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.mg / pyunits.liter,
            doc="Ozone consumption",
        )

        self.electricity = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Ozone generation power demand",
        )

        @self.Constraint(self.flowsheet().time, doc="Ozone consumption constraint")
        def ozone_consumption_constraint(b, t):
            return (
                b.ozone_consumption[t]
                == (
                    (
                        pyunits.convert(
                            b.properties_in[t].conc_mass_comp["toc"],
                            to_units=pyunits.mg / pyunits.liter,
                        )
                        + self.concentration_time[t] / self.contact_time[t]
                    )
                )
                / self.mass_transfer_efficiency[t]
            )

        @self.Constraint(self.flowsheet().time, doc="Ozone mass flow constraint")
        def ozone_flow_mass_constraint(b, t):
            return b.ozone_flow_mass[t] == pyunits.convert(
                b.properties_in[t].flow_vol * b.ozone_consumption[t],
                to_units=pyunits.lb / pyunits.hr,
            )

        @self.Constraint(self.flowsheet().time, doc="Ozone power constraint")
        def electricity_constraint(b, t):
            return b.electricity[t] == (
                b.specific_energy_coeff[t] * b.ozone_flow_mass[t]
            )

        self._perf_var_dict["Ozone Contact Time (min)"] = self.contact_time
        self._perf_var_dict["Ozone CT Value ((mg*min)/L)"] = self.concentration_time
        self._perf_var_dict[
            "Ozone Mass Transfer Efficiency"
        ] = self.mass_transfer_efficiency
        self._perf_var_dict["Ozone Mass Flow (lb/hr)"] = self.ozone_flow_mass
        self._perf_var_dict["Ozone Unit Power Demand (kW)"] = self.electricity

    @property
    def default_costing_method(self):
        return self.cost_ozonation

    @staticmethod
    def cost_ozonation(blk):
        """
        General method for costing ozone addition. Capital cost is
        based on the inlet flowrate and dosage of ozone.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "ozone_capital_a_parameter",
                "ozone_capital_b_parameter",
                "ozone_capital_c_parameter",
                "ozone_capital_d_parameter",
            ],
        )
        # Get costing term for ozone addition
        expr = blk.unit_model._get_ozone_capital_cost(blk, A, B, C, D)

        # Add cost variable
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
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

    @staticmethod
    def _get_ozone_capital_cost(blk, A, B, C, D):
        """
        Generate expressions for capital cost of ozonation system.
        """
        t0 = blk.flowsheet().time.first()

        ln_Q = pyo.log(
            pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol
                / (pyo.units.m**3 / pyo.units.hour),
                to_units=pyo.units.dimensionless,
            )
        )
        dosage = pyo.units.convert(
            blk.unit_model.ozone_consumption[t0] / (pyo.units.mg / pyo.units.liter),
            to_units=pyo.units.dimensionless,
        )

        expr = pyo.units.convert(
            A + B * dosage + C * ln_Q + D * dosage * ln_Q,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        return expr
