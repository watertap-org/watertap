###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains a zero-order representation of a coagulation/flocculation unit.
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var, Expression
from idaes.core import declare_process_block_class

from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

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

        self.alum_dose = Var(self.flowsheet().time,
                             units=pyunits.mg / pyunits.L,
                             bounds=(0, None),
                             doc="Dosing rate of alum")

        self.polymer_dose = Var(self.flowsheet().time,
                             units=pyunits.mg / pyunits.L,
                             bounds=(0, None),
                             doc="Dosing rate of polymer")

        self.anion_to_cation_polymer_ratio = Var(
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Ratio of anionic to cationic polymer in dosage")

        self.chemical_flow_vol = Var(
            self.flowsheet().time,
            ["alum", "polymer"],
            units=pyunits.m ** 3 / pyunits.s,
            bounds=(0, None),
            doc="Volumetric flow rate of chemical solution")

        self.rapid_mix_retention_time = Var(self.flowsheet().config.time,
                                            units=pyunits.seconds,
                                            doc="Rapid Mix Retention Time")

        self.floc_retention_time = Var(self.flowsheet().config.time,
                                       units=pyunits.minutes,
                                       doc="Floc Retention Time")

        self.rapid_mix_basin_vol = Var(self.flowsheet().config.time,
                                          units=pyunits.gallons,
                                          doc="Rapid Mix Basin Volume")

        self.floc_basin_vol = Var(self.flowsheet().config.time,
                                  units=pyunits.gallons,
                                  doc="Floc Basin Volume")

        self._fixed_perf_vars.append(self.alum_dose)
        self._fixed_perf_vars.append(self.polymer_dose)
        self._fixed_perf_vars.append(self.anion_to_cation_polymer_ratio)
        self._fixed_perf_vars.append(self.rapid_mix_retention_time)
        self._fixed_perf_vars.append(self.floc_retention_time)

        self._perf_var_dict["Alum Dosage (mg/L)"] = self.alum_dose
        self._perf_var_dict["Polymer Dosage (mg/L)"] = self.polymer_dose
        self._perf_var_dict["Alum Flow (lb/h)"] = self.chemical_flow_vol[0, "alum"]
        self._perf_var_dict["Polymer Flow (lb/d)"] = self.chemical_flow_vol[0, "polymer"]
        self._perf_var_dict["Rapid Mix Basin Volume (gal)"] = self.rapid_mix_basin_vol
        self._perf_var_dict["Floc Basin Volume (gal)"] = self.floc_basin_vol
        self._perf_var_dict["Rapid Mix Retention Time (s)"] = self.rapid_mix_retention_time
        self._perf_var_dict["Floc Retention Time (min)"] = self.floc_retention_time

        def rule_rapid_mix_basin_vol(blk, t):
            return (blk.rapid_mix_basin_vol[t] ==
                    pyunits.convert(blk.properties[t].flow_vol
                                    * blk.rapid_mix_retention_time[t],
                                    to_units=pyunits.gallons))
        self.rapid_mix_basin_vol_constraint = Constraint(self.flowsheet().config.time, rule=rule_rapid_mix_basin_vol)

        def rule_floc_basin_vol(blk, t):
            return (blk.floc_basin_vol[t] ==
                    pyunits.convert(blk.properties[t].flow_vol
                                    * blk.floc_retention_time[t],
                                    to_units=pyunits.gallons))
        self.floc_basin_vol_constraint = Constraint(self.flowsheet().config.time, rule=rule_floc_basin_vol)

        def rule_chem_flow(blk, t, j):
            if j == 'alum':
                chemical_dosage = blk.alum_dose[t]
                units= pyunits.lb/pyunits.hour
            elif j == 'polymer':
                chemical_dosage = blk.polymer_dose[t]
                units = pyunits.lb / pyunits.day
            return (blk.chemical_flow_vol[t, j] ==
                    pyunits.convert(
                        chemical_dosage * blk.properties[t].flow_vol,
                        to_units=units))
        self.chemical_flow_constraint = Constraint(self.flowsheet().time,
                                                   ["alum", "polymer"],
                                                   rule=rule_chem_flow)

        def rule_anionic_polymer_dose(blk, t):
            return (blk.anion_to_cation_polymer_ratio * blk.polymer_dose[t]
                    / (blk.anion_to_cation_polymer_ratio + 1))
        self.anionic_polymer_dose = Expression(self.flowsheet().config.time,
                                               rule=rule_anionic_polymer_dose)

        def rule_cationic_polymer_dose(blk, t):
            return (blk.polymer_dose[t]
                    / (blk.anion_to_cation_polymer_ratio + 1))
        self.cationic_polymer_dose = Expression(self.flowsheet().config.time,
                                                rule=rule_cationic_polymer_dose)

        self._perf_expr_dict["Anionic Polymer Dosage (mg/L)"] = self.anionic_polymer_dose
        self._perf_expr_dict["Cationic Polymer Dosage (mg/L)"] = self.cationic_polymer_dose

        # pump_electricity(self, self.chemical_flow_vol)





