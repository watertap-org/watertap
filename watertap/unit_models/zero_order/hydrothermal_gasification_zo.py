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
This module contains a zero-order representation of a hydrothermal gasification unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("HTGZO")
class HTGZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a hydrothermal gasification (HTG) unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "hydrothermal_gasification"

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
            doc="Constraint for electricity consumption based on inlet flowrate.",
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
        return self.cost_hydrothermal_gasification

    @staticmethod
    def cost_hydrothermal_gasification(blk):
        """
        General method for costing hydrothermal gasification unit. Capital cost
        is based on the CHG reactor and other wastewater treatment equipment including
        a feed pump, a booster pump, a feed/product exchanger, a fired heater,
        a hydrocyclone, and a product air fin cooler.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            IF_reactor,
            EP_reactor,
            F0_reactor,
            SE_reactor,
            IF_pump,
            EP_pump,
            F0_pump,
            SE_pump,
            IF_booster,
            EP_booster,
            F0_booster,
            SE_booster,
            IF_hydrocyclone,
            EP_hydrocyclone,
            F0_hydrocyclone,
            SE_hydrocyclone,
            IF_cooler,
            EP_cooler,
            F0_cooler,
            SE_cooler,
            IF_exchanger,
            EP_exchanger,
            F0_exchanger,
            Fnew_exchanger,
            SE_exchanger,
            IF_heater,
            EP_heater,
            F0_heater,
            Fnew_heater,
            SE_heater,
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
                "installation_factor_booster",
                "equipment_cost_booster",
                "base_flowrate_booster",
                "scaling_exponent_booster",
                "installation_factor_hydrocyclone",
                "equipment_cost_hydrocyclone",
                "base_flowrate_hydrocyclone",
                "scaling_exponent_hydrocyclone",
                "installation_factor_cooler",
                "equipment_cost_cooler",
                "base_flowrate_cooler",
                "scaling_exponent_cooler",
                "installation_factor_exchanger",
                "equipment_cost_exchanger",
                "base_area_exchanger",
                "new_area_exchanger",
                "scaling_exponent_exchanger",
                "installation_factor_heater",
                "equipment_cost_heater",
                "base_heat_duty_heater",
                "new_heat_duty_heater",
                "scaling_exponent_heater",
            ],
        )

        sizing_term_reactor = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_reactor),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_pump = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_pump),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_booster = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_booster),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_hydrocyclone = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_hydrocyclone),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_cooler = pyo.units.convert(
            (blk.unit_model.flow_mass_in[t0] / F0_cooler),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_exchanger = pyo.units.convert(
            (Fnew_exchanger / F0_exchanger),
            to_units=pyo.units.dimensionless,
        )

        sizing_term_heater = pyo.units.convert(
            (Fnew_heater / F0_heater),
            to_units=pyo.units.dimensionless,
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        reactor_cost = pyo.units.convert(
            IF_reactor * EP_reactor * sizing_term_reactor**SE_reactor,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        pump_cost = pyo.units.convert(
            IF_pump * EP_pump * sizing_term_pump**SE_pump,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        booster_cost = pyo.units.convert(
            IF_booster * EP_booster * sizing_term_booster**SE_booster,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        hydrocyclone_cost = pyo.units.convert(
            IF_hydrocyclone
            * EP_hydrocyclone
            * sizing_term_hydrocyclone**SE_hydrocyclone,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        cooler_cost = pyo.units.convert(
            IF_cooler * EP_cooler * sizing_term_cooler**SE_cooler,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        exchanger_cost = pyo.units.convert(
            IF_exchanger * EP_exchanger * sizing_term_exchanger**SE_exchanger,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        heater_cost = pyo.units.convert(
            IF_heater * EP_heater * sizing_term_heater**SE_heater,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = (
            reactor_cost
            + pump_cost
            + booster_cost
            + hydrocyclone_cost
            + cooler_cost
            + exchanger_cost
            + heater_cost
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
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.catalyst_flow[t0], "catalyst_HTG"
        )
