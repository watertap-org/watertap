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
This module contains a zero-order representation of an autothermal hydrothermal liquefaction unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("ATHTLZO")
class ATHTLZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an autothermal hydrothermal liquefaction (AT-HTL) unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "autothermal_hydrothermal_liquefaction"

        build_sido_reactive(self)

        self.flow_mass_in = Var(
            self.flowsheet().time,
            units=pyunits.t / pyunits.hour,
            bounds=(0, None),
            doc="Inlet mass flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for inlet mass flowrate.",
        )
        def cons_flow_mass(b, t):
            return b.flow_mass_in[t] == pyunits.convert(
                sum(
                    b.properties_in[t].flow_mass_comp[j]
                    for j in b.properties_in[t].component_list
                ),
                to_units=pyunits.t / pyunits.hour,
            )

        self._perf_var_dict["Inlet Mass Flowrate"] = self.flow_mass_in

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self._perf_var_dict["Electricity Demand"] = self.electricity

        self.energy_electric_flow_mass = Var(
            units=pyunits.kWh / pyunits.t,
            doc="Electricity intensity with respect to inlet flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on inlet flow rate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == pyunits.convert(
                b.energy_electric_flow_mass * b.flow_mass_in[t], to_units=pyunits.kW
            )

        self._fixed_perf_vars.append(self.energy_electric_flow_mass)
        self._perf_var_dict["Electricity Intensity"] = self.energy_electric_flow_mass

        self.catalyst_dosage = Var(
            units=pyunits.pound / pyunits.t,
            bounds=(0, None),
            doc="Dosage of catalyst per inlet flow",
        )

        self._fixed_perf_vars.append(self.catalyst_dosage)

        self._perf_var_dict["Dosage of catalyst per inlet flow"] = self.catalyst_dosage

        self.catalyst_flow = Var(
            self.flowsheet().time,
            units=pyunits.pound / pyunits.hr,
            bounds=(0, None),
            doc="Catalyst flow",
        )
        self._perf_var_dict["Catalyst flow"] = self.catalyst_flow

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for catalyst flow based on inlet flow rate.",
        )
        def eq_catalyst_flow(b, t):
            return b.catalyst_flow[t] == pyunits.convert(
                b.catalyst_dosage * b.flow_mass_in[t],
                to_units=pyunits.pound / pyunits.hr,
            )

    @property
    def default_costing_method(self):
        return self.cost_autothermal_hydrothermal_liquefaction

    @staticmethod
    def cost_autothermal_hydrothermal_liquefaction(blk):
        """
        General method for costing autothermal-hydrothermal liquefaction unit. Capital cost
        is based on the HTL reactor, booster pump, solid filter, other equipment, and
        heat oil system.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            A,
            B,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
            J,
            K,
            L,
            M,
            N,
            O,
            P,
            Q,
            R,
            S,
            T,
        ) = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "installation_factor_reactor",
                "equipment_cost_reactor",
                "base_flowrate_reactor",
                "scaling_exponent_reactor",
                "installation_factor_pump",
                "equipment_cost_pump",
                "base_flowrate_pump",
                "scaling_exponent_pump",
                "installation_factor_other",
                "equipment_cost_other",
                "base_flowrate_other",
                "scaling_exponent_other",
                "installation_factor_solid_filter",
                "equipment_cost_solid_filter",
                "base_flowrate_solid_filter",
                "scaling_exponent_solid_filter",
                "installation_factor_heat",
                "equipment_cost_heat",
                "base_flowrate_heat",
                "scaling_exponent_heat",
            ],
        )

        sizing_term_reactor = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / C),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_pump = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / G),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_other = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / K),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_solid_filter = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / O),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_heat = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / S),
            to_units=pyo.units.dimensionless,
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

        reactor_cost = pyo.units.convert(
            A * B * sizing_term_reactor**D,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        pump_cost = pyo.units.convert(
            E * F * sizing_term_pump**H,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        other_cost = pyo.units.convert(
            I * J * sizing_term_other**L,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        solid_filter_cost = pyo.units.convert(
            M * N * sizing_term_solid_filter**P,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        heat_cost = pyo.units.convert(
            Q * R * sizing_term_heat**T,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = reactor_cost + pump_cost + other_cost + solid_filter_cost + heat_cost

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
            blk.unit_model.catalyst_flow[t0], "catalyst_ATHTL"
        )
