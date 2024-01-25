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
This module contains a zero-order representation of a brine concentrator unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("BrineConcentratorZO")
class BrineConcentratorZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a brine concentrator unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "brine_concentrator"

        build_sido(self)

        if "tds" not in self.config.property_package.solute_set:
            raise KeyError(
                "TDS must be included in the solute list for determining"
                " electricity intensity and power consumption of the brine "
                "concentrator unit."
            )

        # Fitting parameters based on regressions for capital and electricity
        # developed from data in Table 5.1, Table A2.3 in:
        # Survey of High-Recovery and Zero Liquid Discharge Technologies for
        # Water Utilities (2008). WateReuse Foundation:
        # https://www.waterboards.ca.gov/water_issues/programs/grants_loans/water_recycling/research/02_006a_01.pdf
        # Capital = f(TDS, recovery, flow)
        # Electricity = f(TDS, recovery, flow)
        self.elec_coeff_1 = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Constant 1 in electricity intensity equation",
        )
        self.elec_coeff_2 = Var(
            units=pyunits.L / pyunits.mg * pyunits.kWh / pyunits.m**3,
            doc="Constant 2 in electricity intensity equation",
        )
        self.elec_coeff_3 = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Constant 3 in electricity intensity equation",
        )
        self.elec_coeff_4 = Var(
            units=pyunits.kWh / pyunits.m**6 * pyunits.hour,
            doc="Constant 4 in electricity intensity equation",
        )

        self._fixed_perf_vars.append(self.elec_coeff_1)
        self._fixed_perf_vars.append(self.elec_coeff_2)
        self._fixed_perf_vars.append(self.elec_coeff_3)
        self._fixed_perf_vars.append(self.elec_coeff_4)

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of brine concentrator",
        )
        self.electricity_intensity = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific energy consumption with respect to feed flowrate",
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Electricity intensity constraint"
        )
        def electricity_intensity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            tds_in = pyunits.convert(
                b.properties_in[t].conc_mass_comp["tds"],
                to_units=pyunits.mg / pyunits.L,
            )
            return (
                b.electricity_intensity[t]
                == b.elec_coeff_1
                + b.elec_coeff_2 * tds_in
                + b.elec_coeff_3 * b.recovery_frac_mass_H2O[t]
                + b.elec_coeff_4 * q_in
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Power consumption constraint"
        )
        def electricity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity[t] == b.electricity_intensity[t] * q_in

        self._perf_var_dict["Power Consumption (kW)"] = self.electricity
        self._perf_var_dict[
            "Electricity intensity per Inlet Flowrate  (kWh/m3)"
        ] = self.electricity_intensity

    @property
    def default_costing_method(self):
        return self.cost_brine_concentrator

    @staticmethod
    def cost_brine_concentrator(blk):
        """
        General method for costing brine concentration processes. Capital cost
        is based on the volumetirc flowrate and TDS of the incoming stream and
        the water recovery.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()
        inlet_state = blk.unit_model.properties_in[t0]

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
                "capital_a_parameter",
                "capital_b_parameter",
                "capital_c_parameter",
                "capital_d_parameter",
            ],
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

        expr = (
            pyo.units.convert(
                A, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                B * inlet_state.conc_mass_comp["tds"],
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                C * blk.unit_model.recovery_frac_mass_H2O[t0],
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                D * inlet_state.flow_vol,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
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
