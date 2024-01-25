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
This module contains a zero-order representation of an ion exchange unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import Reference, units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("IonExchangeZO")
class IonExchangeZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an Ion exchange unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "ion_exchange"

        build_sido(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)

        # mutable parameter; default value found in WT3 for anion exchange
        if self.config.process_subtype == "clinoptilolite":
            pass
        else:
            self.eta_pump.set_value(0.8)
            # mutable parameter; default value of 2 bar converted to feet head
            self.lift_height.set_value(69.91052 * pyunits.feet)

        # Add variables and constraints for material requirements
        self.NaCl_flowrate = Var(
            self.flowsheet().time,
            initialize=1,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Flowrate of NaCl addition",
        )
        self.NaCl_dose = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Dosage of NaCl addition",
        )

        self._fixed_perf_vars.append(self.NaCl_dose)
        self._perf_var_dict["NaCl Addition"] = self.NaCl_flowrate

        @self.Constraint(self.flowsheet().time)
        def NaCl_constraint(blk, t):
            return blk.NaCl_flowrate[t] == blk.NaCl_dose * blk.properties_in[t].flow_vol

        self.resin_demand = Var(
            self.flowsheet().time,
            initialize=1,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Replacement rate of ion exchange resin",
        )
        self.resin_replacement = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Resin replacement as a function of flow",
        )

        self._fixed_perf_vars.append(self.resin_replacement)
        self._perf_var_dict["Resin Demand"] = self.resin_demand

        @self.Constraint(self.flowsheet().time)
        def resin_constraint(blk, t):
            return (
                blk.resin_demand[t]
                == blk.resin_replacement * blk.properties_in[t].flow_vol
            )

        if self.config.process_subtype == "clinoptilolite":
            if "ammonium_as_nitrogen" in self.config.property_package.solute_set:
                self.nitrogen_clay_ratio = Var(
                    self.flowsheet().config.time,
                    units=pyunits.dimensionless,
                    doc="Mass fraction of nitrogen in clay mixture",
                )

                self._fixed_perf_vars.append(self.nitrogen_clay_ratio)

                self.final_solids_mass = Var(
                    self.flowsheet().config.time,
                    units=pyunits.kg / pyunits.s,
                    doc="Solids mass flow in byproduct stream",
                )

                @self.Constraint(
                    self.flowsheet().time,
                    doc="Solids mass flow in byproduct stream constraint",
                )
                def solids_mass_flow_constraint(b, t):
                    return (
                        b.final_solids_mass[t]
                        == b.properties_byproduct[t].flow_mass_comp[
                            "ammonium_as_nitrogen"
                        ]
                        / b.nitrogen_clay_ratio[t]
                    )

                self._perf_var_dict[
                    "Nitrogen-Clay Mixture Ratio (kg/kg)"
                ] = self.nitrogen_clay_ratio
                self._perf_var_dict[
                    "Final mass flow of clay and nitrogen (kg/s)"
                ] = self.final_solids_mass
            else:
                raise KeyError(
                    "ammonium_as_nitrogen should be defined in solute_list for this subtype."
                )

    @property
    def default_costing_method(self):
        return self.cost_ion_exchange

    @staticmethod
    def cost_ion_exchange(blk):
        """
        Two methods for costing ion exchange:
        (1) General method for costing ion exchange units. Capital cost is based on
        the both inlet flow and TDS.
        This method also registers the NaCl demand, resin replacement and
        electricity demand as costed flows.
        (2) General method using unit capex and unit opex cost parameters, tailored for AMO
        wastewater resource recovery (process subtype: clinoptilolite)
        """
        t0 = blk.flowsheet().time.first()

        if blk.unit_model.config.process_subtype != "clinoptilolite":
            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )
            # Get costing parameter sub-block for this technology
            A, B, C, D = blk.unit_model._get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                [
                    "capital_a_parameter",
                    "capital_b_parameter",
                    "capital_c_parameter",
                    "capital_d_parameter",
                ],
            )

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
            T = pyo.units.convert(
                blk.unit_model.properties_in[t0].conc_mass_comp["tds"]
                / (pyo.units.mg / pyo.units.liter),
                to_units=pyo.units.dimensionless,
            )

            expr = pyo.units.convert(
                pyo.exp(A + B * ln_Q + C * T + D * ln_Q * T) * pyo.units.USD_2017,
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
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.NaCl_flowrate[t0], "sodium_chloride"
            )
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.resin_demand[t0], "ion_exchange_resin"
            )
        else:
            # Get parameter dict from database
            parameter_dict = (
                blk.unit_model.config.database.get_unit_operation_parameters(
                    blk.unit_model._tech_type,
                    subtype=blk.unit_model.config.process_subtype,
                )
            )
            # Get costing parameter sub-block for this technology
            unit_capex, unit_opex = blk.unit_model._get_tech_parameters(
                blk,
                parameter_dict,
                blk.unit_model.config.process_subtype,
                ["unit_capex", "unit_opex"],
            )

            # Add cost variable and constraint
            blk.capital_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency,
                bounds=(0, None),
                doc="Capital cost of unit operation",
            )

            capex_expr = pyo.units.convert(
                blk.unit_model.properties_in[t0].flow_vol * unit_capex,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )

            # Determine if a costing factor is required
            blk.costing_package.add_cost_factor(
                blk, parameter_dict["capital_cost"]["cost_factor"]
            )

            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost == blk.cost_factor * capex_expr
            )

            # Add fixed operating cost variable and constraint
            blk.fixed_operating_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
                bounds=(0, None),
                doc="Fixed operating cost of unit",
            )
            blk.fixed_operating_cost_constraint = pyo.Constraint(
                expr=blk.fixed_operating_cost
                == pyo.units.convert(
                    blk.unit_model.properties_in[t0].flow_vol * unit_opex,
                    to_units=blk.config.flowsheet_costing_block.base_currency
                    / blk.config.flowsheet_costing_block.base_period,
                )
            )

            # Register flows
            blk.config.flowsheet_costing_block.cost_flow(
                blk.unit_model.electricity[t0], "electricity"
            )
