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
This module contains a zero-order representation of a crystallizer unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData

__author__ = "Adam Atia, Kurban Sitterley"


@declare_process_block_class("CrystallizerZO")
class CrystallizerZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a crystallizer unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "crystallizer"

        build_sido(self)

        if "tds" not in self.config.property_package.solute_set:
            raise KeyError(
                "TDS must be included in the solute list for determining"
                " electricity intensity and power consumption of the crystallizer unit."
            )

        # Fitting parameters based on regressions for capital and electricity
        # developed from data in Table 4.4, Table A2.1, Table A2.3 in:
        # Survey of High-Recovery and Zero Liquid Discharge Technologies for
        # Water Utilities (2008). WateReuse Foundation:
        # https://www.waterboards.ca.gov/water_issues/programs/grants_loans/water_recycling/research/02_006a_01.pdf
        # Capital = f(TDS, purge fraction, flow)
        # Electricity = f(TDS, purge fraction, flow)
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
            doc="Power consumption of crystallizer",
        )
        self.electricity_intensity = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific energy consumption with respect to feed flowrate",
        )

        @self.Expression(
            doc="Purge fraction for crystallizer",
        )
        def purge_fraction(b):
            return 1 - b.recovery_frac_mass_H2O[0]

        @self.Constraint(doc="Electricity intensity constraint")
        def electricity_intensity_constraint(b):
            q_in = pyunits.convert(
                b.properties_in[0].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            tds_in = pyunits.convert(
                b.properties_in[0].conc_mass_comp["tds"],
                to_units=pyunits.mg / pyunits.L,
            )
            return (
                b.electricity_intensity
                == b.elec_coeff_1
                + b.elec_coeff_2 * tds_in
                + b.elec_coeff_3 * b.purge_fraction
                + b.elec_coeff_4 * q_in
            )

        @self.Constraint(doc="Power consumption constraint")
        def electricity_constraint(b):
            q_in = pyunits.convert(
                b.properties_in[0].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity[0] == b.electricity_intensity * q_in

        self._perf_var_dict["Power Consumption (kW)"] = self.electricity
        self._perf_var_dict["Electricity intensity per Inlet Flowrate  (kWh/m3)"] = (
            self.electricity_intensity
        )

    @property
    def default_costing_method(self):
        return self.cost_crystallizer

    @staticmethod
    def cost_crystallizer(blk):
        """
        General method for costing crystallizer processes. Capital cost
        is based on the volumetric flowrate and TDS of the incoming stream and
        the purge fraction.
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
                C * blk.unit_model.purge_fraction,
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
