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
This module contains a zero-order representation of a sedimentation unit
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("SedimentationZO")
class SedimentationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Sedimentation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "sedimentation"

        build_sido(self)
        constant_intensity(self)

        # TODO: Does it really make sense for this to be indexed by time?
        self.basin_surface_area = Var(
            self.flowsheet().config.time,
            units=pyunits.ft**2,
            doc="Surface area of sedimentation tank",
        )

        self.settling_velocity = Var(
            self.flowsheet().config.time,
            units=pyunits.m / pyunits.s,
            doc="Particle settling velocity",
        )

        self._fixed_perf_vars.append(self.settling_velocity)

        self._perf_var_dict["Basin Surface Area (ft^2)"] = self.basin_surface_area
        self._perf_var_dict["Settling Velocity (m/s)"] = self.settling_velocity

        def rule_basin_surface_area(b, t):
            return b.basin_surface_area[t] == pyunits.convert(
                b.properties_in[t].flow_vol / b.settling_velocity[t],
                to_units=pyunits.ft**2,
            )

        self.basin_surface_area_constraint = Constraint(
            self.flowsheet().time, rule=rule_basin_surface_area
        )

        if self.config.process_subtype == "phosphorus_capture":
            self.phosphorus_solids_ratio = Var(
                self.flowsheet().config.time,
                units=pyunits.dimensionless,
                doc="Mass fraction of phosphorus in settleable solids",
            )

            self._fixed_perf_vars.append(self.phosphorus_solids_ratio)
            self._perf_var_dict[
                "Phosphorus-Solids Ratio (kg/kg)"
            ] = self.phosphorus_solids_ratio

            # This subtype is intended to be used explicitly for phosphorous capture.
            # If the user provides TSS, the amount of settled phosphate would be determined based on
            # an assumed fraction of phosphate in TSS. Alternatively, the user could provide phosphates
            # as the species, and the amount of solids + phosphate settled would be reported.
            # However, the user cannot provide both TSS and phosphates.
            if (
                "phosphates" in self.config.property_package.solute_set
                and "tss" in self.config.property_package.solute_set
            ):
                raise KeyError(
                    "tss and phosphates cannot both be defined in the solute_list. "
                    "Please choose one."
                )
            elif "phosphates" in self.config.property_package.solute_set:
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
                        == b.properties_byproduct[t].flow_mass_comp["phosphates"]
                        / b.phosphorus_solids_ratio[t]
                    )

                self._perf_var_dict[
                    "Final mass flow of settled solids (kg/s)"
                ] = self.final_solids_mass

            elif "tss" in self.config.property_package.solute_set:
                self.final_phosphate_mass = Var(
                    self.flowsheet().config.time,
                    units=pyunits.kg / pyunits.s,
                    doc="Phosphate mass flow in byproduct stream",
                )

                @self.Constraint(
                    self.flowsheet().time,
                    doc="Phosphate mass flow in byproduct stream constraint",
                )
                def phosphate_mass_flow_constraint(b, t):
                    return (
                        b.final_phosphate_mass[t]
                        == b.properties_byproduct[t].flow_mass_comp["tss"]
                        * b.phosphorus_solids_ratio[t]
                    )

                self._perf_var_dict[
                    "Final mass flow of settled phosphate (kg/s)"
                ] = self.final_phosphate_mass

            else:
                # Raise this error in case the user is intended to make use of the subtype but entered
                # the wrong component names.
                raise KeyError(
                    "One of the following should be specified in the solute_list: "
                    "tss or phosphates"
                )
