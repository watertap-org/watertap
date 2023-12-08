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
This module contains a zero-order representation of a granular activated carbon unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from idaes.core.util.math import smooth_min
from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("GACZO")
class GACZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a granular activated carbon unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "gac"

        build_sido(self)

        # Empty Bed Contact Time
        self.empty_bed_contact_time = Var(
            units=pyunits.hour, bounds=(0, None), doc="Empty bed contact time of unit"
        )

        self._fixed_perf_vars.append(self.empty_bed_contact_time)
        self._perf_var_dict["Empty Bed Contact Time"] = self.empty_bed_contact_time

        # Electricity Demand
        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self.electricity_intensity_parameter = Var(
            units=pyunits.kW / pyunits.m**3,
            doc="Parameter for calculating electricity based on empty bed "
            "contact time",
        )

        self.energy_electric_flow_vol_inlet = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Electricity intensity with respect to inlet flowrate of unit",
        )

        @self.Constraint(doc="Electricity intensity based on empty bed contact time.")
        def electricity_intensity_constraint(b):
            return (
                b.energy_electric_flow_vol_inlet
                == b.electricity_intensity_parameter * b.empty_bed_contact_time
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "feed flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_flow_vol_inlet
                * pyunits.convert(
                    b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.hour
                )
            )

        self._fixed_perf_vars.append(self.electricity_intensity_parameter)

        self._perf_var_dict["Electricity Demand"] = self.electricity
        self._perf_var_dict[
            "Electricity Intensity"
        ] = self.energy_electric_flow_vol_inlet

        # Demand for activated carbon
        self.activated_carbon_replacement = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Replacement rate of activated carbon",
        )

        self.activated_carbon_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hour,
            bounds=(0, None),
            doc="Demand for activated carbon",
        )

        # Demand for activated carbon
        self.activated_carbon_bulk_density = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Bulk density, total mass of GAC per total bed volume",
        )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for activated carbon consumption."
        )
        def activated_carbon_equation(b, t):
            return b.activated_carbon_demand[t] == (
                b.activated_carbon_replacement
                * pyunits.convert(
                    b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.hour
                )
            )

        self._fixed_perf_vars.append(self.activated_carbon_replacement)
        self._fixed_perf_vars.append(self.activated_carbon_bulk_density)

        self._perf_var_dict["Activated Carbon Demand"] = self.activated_carbon_demand

    @property
    def default_costing_method(self):
        return self.cost_gac

    @staticmethod
    def cost_gac(blk, number_of_parallel_units=5):
        """
        Adapted from core GAC costing model initially released in v0.6.0

        3 equation capital cost estimation for GAC systems with: (i), contactor/pressure vessel cost by polynomial
        as a function of individual contactor volume; (ii), initial charge of GAC adsorbent cost by exponential as a
        function of required mass of GAC adsorbent; and (iii), other process costs (vessels, pipes, instrumentation, and
        controls) calculated by power law as a function of total contactor(s) volume. Operating costs calculated as the
        required makeup and regeneration of GAC adsorbent. Energy for backwash and booster pumps considered negligible
        compared to regeneration costs
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units in operation (default: 5)
        """
        t0 = blk.flowsheet().time.first()

        Q = blk.unit_model.properties_in[t0].flow_vol
        T = blk.unit_model.empty_bed_contact_time
        gac_dens = blk.unit_model.activated_carbon_bulk_density

        # total bed volume
        V = pyo.units.convert(
            Q, to_units=pyo.units.m**3 / pyo.units.seconds
        ) * pyo.units.convert(T, to_units=pyo.units.seconds)

        # mass of gac in bed
        bed_mass_gac = pyo.units.convert(
            V, to_units=pyo.units.m**3
        ) * pyo.units.convert(gac_dens, to_units=pyo.units.kg / pyo.units.m**3)

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Call contactor, adsorbent (GAC) initial charge cost, and other process cost parameters
        (
            A0,
            A1,
            A2,
            A3,
            B0,
            B1,
            C0,
            C1,
            bed_mass_max_ref,
        ) = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "contactor_cost_coeff_0",
                "contactor_cost_coeff_1",
                "contactor_cost_coeff_2",
                "contactor_cost_coeff_3",
                "adsorbent_unit_cost_coeff",
                "adsorbent_unit_cost_exp_coeff",
                "other_cost_coeff",
                "other_cost_exp",
                "bed_mass_max_ref",
            ],
        )

        contactor_cost = number_of_parallel_units * pyo.units.convert(
            (
                A3 * (V / number_of_parallel_units) ** 3
                + A2 * (V / number_of_parallel_units) ** 2
                + A1 * (V / number_of_parallel_units) ** 1
                + A0
            ),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        bed_mass_gac_ref = (
            smooth_min(
                bed_mass_max_ref / pyo.units.kg,
                pyo.units.convert(bed_mass_gac, to_units=pyo.units.kg) / pyo.units.kg,
            )
            * pyo.units.kg
        )

        adsorbent_unit_cost = pyo.units.convert(
            B0 * pyo.exp(bed_mass_gac_ref * B1),
            to_units=blk.config.flowsheet_costing_block.base_currency
            * pyo.units.kg**-1,
        )

        adsorbent_cost = adsorbent_unit_cost * bed_mass_gac

        other_process_cost = pyo.units.convert(
            (C0 * ((pyo.units.m**3) ** -C1) * V**C1),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = contactor_cost + adsorbent_cost + other_process_cost

        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        # electricity was not implemented in core GAC but retained in ZO update
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.activated_carbon_demand[t0], "activated_carbon"
        )
