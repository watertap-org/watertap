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
This module contains a zero-order representation of a coagulation/flocculation unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("CoagulationFlocculationZO")
class CoagulationFlocculationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Coagulation/Flocculation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "coag_and_floc"

        build_pt(self)

        self.alum_dose = Var(
            self.flowsheet().time,
            units=pyunits.mg / pyunits.L,
            bounds=(0, None),
            doc="Dosing rate of alum",
        )

        self.polymer_dose = Var(
            self.flowsheet().time,
            units=pyunits.mg / pyunits.L,
            bounds=(0, None),
            doc="Dosing rate of polymer",
        )

        self.anion_to_cation_polymer_ratio = Var(
            self.flowsheet().time,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Ratio of anionic to cationic polymer in dosage",
        )

        self.anionic_polymer_dose = Var(
            self.flowsheet().config.time,
            bounds=(0, None),
            units=pyunits.mg / pyunits.L,
            doc="Dosing rate of anionic polymer",
        )

        self.cationic_polymer_dose = Var(
            self.flowsheet().config.time,
            bounds=(0, None),
            units=pyunits.mg / pyunits.L,
            doc="Dosing rate of cationic polymer",
        )

        self.chemical_flow_mass = Var(
            self.flowsheet().time,
            ["alum", "polymer"],
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Mass flow rate of chemical solution",
        )

        self.rapid_mix_retention_time = Var(
            self.flowsheet().config.time,
            units=pyunits.seconds,
            doc="Rapid Mix Retention Time",
        )

        self.floc_retention_time = Var(
            self.flowsheet().config.time,
            units=pyunits.minutes,
            doc="Floc Retention Time",
        )

        self.rapid_mix_basin_vol = Var(
            units=pyunits.m**3, doc="Rapid Mix Basin Volume"
        )

        self.floc_basin_vol = Var(units=pyunits.m**3, doc="Floc Basin Volume")

        self.num_rapid_mixers = Var(
            units=pyunits.dimensionless, doc="Number of Rapid Mixers"
        )

        self.num_floc_mixers = Var(
            units=pyunits.dimensionless, doc="Number of Floc Mixers"
        )

        self.num_rapid_mix_processes = Var(
            units=pyunits.dimensionless, doc="Number of Rapid Mix Processes"
        )

        self.num_floc_processes = Var(
            units=pyunits.dimensionless, doc="Number of Floc Processes"
        )

        self.num_coag_processes = Var(
            units=pyunits.dimensionless, doc="Number of Coagulation Processes"
        )

        self.num_floc_injection_processes = Var(
            units=pyunits.dimensionless, doc="Number of Floc Injection Processes"
        )

        self.velocity_gradient_rapid_mix = Var(
            self.flowsheet().config.time,
            units=pyunits.s**-1,
            doc="Rapid Mix Velocity Gradient",
        )

        self.velocity_gradient_floc = Var(
            self.flowsheet().config.time,
            units=pyunits.s**-1,
            doc="Floc Velocity Gradient",
        )

        self.power_rapid_mix = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            doc="Rapid Mix Power Consumption",
        )

        self.power_floc = Var(
            self.flowsheet().config.time, units=pyunits.kW, doc="Floc Power Consumption"
        )

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            doc="Total Power Consumption",
        )

        self._fixed_perf_vars.append(self.alum_dose)
        self._fixed_perf_vars.append(self.polymer_dose)
        self._fixed_perf_vars.append(self.anion_to_cation_polymer_ratio)
        self._fixed_perf_vars.append(self.rapid_mix_retention_time)
        self._fixed_perf_vars.append(self.floc_retention_time)
        self._fixed_perf_vars.append(self.num_floc_injection_processes)
        self._fixed_perf_vars.append(self.num_floc_processes)
        self._fixed_perf_vars.append(self.num_rapid_mixers)
        self._fixed_perf_vars.append(self.num_coag_processes)
        self._fixed_perf_vars.append(self.num_floc_mixers)
        self._fixed_perf_vars.append(self.num_rapid_mix_processes)
        self._fixed_perf_vars.append(self.velocity_gradient_rapid_mix)
        self._fixed_perf_vars.append(self.velocity_gradient_floc)

        self._perf_var_dict["Alum Dosage (mg/L)"] = self.alum_dose
        self._perf_var_dict["Polymer Dosage (mg/L)"] = self.polymer_dose
        self._perf_var_dict["Alum Flow (kg/s)"] = self.chemical_flow_mass[0, "alum"]
        self._perf_var_dict["Polymer Flow (kg/s)"] = self.chemical_flow_mass[
            0, "polymer"
        ]
        self._perf_var_dict["Rapid Mix Basin Volume (m^3)"] = self.rapid_mix_basin_vol
        self._perf_var_dict["Floc Basin Volume (m^3)"] = self.floc_basin_vol
        self._perf_var_dict[
            "Rapid Mix Retention Time (s)"
        ] = self.rapid_mix_retention_time
        self._perf_var_dict["Floc Retention Time (min)"] = self.floc_retention_time
        self._perf_var_dict[
            "Rapid Mix Velocity Gradient (1/s)"
        ] = self.velocity_gradient_rapid_mix
        self._perf_var_dict[
            "Floc Velocity Gradient (1/s)"
        ] = self.velocity_gradient_floc
        self._perf_var_dict["Rapid Mix Power (kW)"] = self.power_rapid_mix
        self._perf_var_dict["Floc Power (kW)"] = self.power_floc
        self._perf_var_dict["Total Power Consumption (kW)"] = self.electricity

        def rule_rapid_mix_basin_vol(blk):
            return (
                blk.rapid_mix_basin_vol
                == blk.properties[0].flow_vol * blk.rapid_mix_retention_time[0]
            )

        self.rapid_mix_basin_vol_constraint = Constraint(rule=rule_rapid_mix_basin_vol)

        def rule_floc_basin_vol(blk):
            return blk.floc_basin_vol == pyunits.convert(
                blk.properties[0].flow_vol * blk.floc_retention_time[0],
                to_units=pyunits.m**3,
            )

        self.floc_basin_vol_constraint = Constraint(rule=rule_floc_basin_vol)

        def rule_chem_flow(blk, t, j):
            if j == "alum":
                chemical_dosage = blk.alum_dose[t]
            elif j == "polymer":
                chemical_dosage = blk.polymer_dose[t]
            return blk.chemical_flow_mass[t, j] == pyunits.convert(
                chemical_dosage * blk.properties[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        self.chemical_flow_constraint = Constraint(
            self.flowsheet().time, ["alum", "polymer"], rule=rule_chem_flow
        )

        def rule_anionic_polymer_dose(blk, t):
            return blk.anionic_polymer_dose[t] == blk.anion_to_cation_polymer_ratio[
                t
            ] * blk.polymer_dose[t] / (blk.anion_to_cation_polymer_ratio[t] + 1)

        self.anionic_polymer_dose_constraint = Constraint(
            self.flowsheet().config.time, rule=rule_anionic_polymer_dose
        )

        def rule_cationic_polymer_dose(blk, t):
            return blk.cationic_polymer_dose[t] == blk.polymer_dose[t] / (
                blk.anion_to_cation_polymer_ratio[t] + 1
            )

        self.cationic_polymer_dose_constraint = Constraint(
            self.flowsheet().config.time, rule=rule_cationic_polymer_dose
        )

        # TODO: WT3 doesn't include pump electricity. Consider adding after case study validation
        # pump_electricity(self, self.chemical_flow_vol)

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for rapid mix power consumption"
        )
        def rule_power_rapid_mix(b, t):
            return b.power_rapid_mix[t] == pyunits.convert(
                b.num_rapid_mixers
                * b.properties[t].visc_d
                * b.rapid_mix_basin_vol
                * b.velocity_gradient_rapid_mix[t] ** 2,
                to_units=pyunits.kW,
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for floc power consumption"
        )
        def rule_power_floc(b, t):
            return b.power_floc[t] == pyunits.convert(
                b.num_floc_mixers
                * b.properties[t].visc_d
                * b.floc_basin_vol
                * b.velocity_gradient_floc[t] ** 2,
                to_units=pyunits.kW,
            )

        @self.Constraint(self.flowsheet().time, doc="Total power consumption")
        def electricity_constraint(b, t):
            return b.electricity[t] == b.power_floc[t] + b.power_rapid_mix[t]

    @property
    def default_costing_method(self):
        return self.cost_coag_and_floc

    @staticmethod
    def cost_coag_and_floc(blk):
        """
        General method for costing coagulation/flocculation processes. Capital cost
        is based on the alum flowrate and the polymer flowrate of the incoming stream.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F, G, H = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_mix_a_parameter",
                "capital_mix_b_parameter",
                "capital_floc_a_parameter",
                "capital_floc_b_parameter",
                "capital_coag_inj_a_parameter",
                "capital_coag_inj_b_parameter",
                "capital_floc_inj_a_parameter",
                "capital_floc_inj_b_parameter",
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

        cost_rapid_mix = (
            A
            * pyo.units.convert(
                blk.unit_model.rapid_mix_basin_vol, to_units=pyo.units.gallons
            )
            + B
        ) * blk.unit_model.num_rapid_mix_processes

        cost_floc = (
            C
            * pyo.units.convert(
                blk.unit_model.floc_basin_vol, to_units=pyo.units.Mgallons
            )
            + D
        ) * blk.unit_model.num_floc_processes

        cost_coag_inj = (
            E
            * pyo.units.convert(
                blk.unit_model.chemical_flow_mass[t0, "alum"],
                to_units=(pyo.units.lb / pyo.units.hr),
            )
            + F
        ) * blk.unit_model.num_coag_processes

        cost_floc_inj = (
            G
            * pyo.units.convert(
                blk.unit_model.chemical_flow_mass[t0, "polymer"],
                to_units=(pyo.units.lb / pyo.units.day),
            )
            + H
        ) * blk.unit_model.num_floc_injection_processes

        expr = (
            pyo.units.convert(
                cost_rapid_mix,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
            + pyo.units.convert(
                cost_floc, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                cost_coag_inj, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            + pyo.units.convert(
                cost_floc_inj, to_units=blk.config.flowsheet_costing_block.base_currency
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
