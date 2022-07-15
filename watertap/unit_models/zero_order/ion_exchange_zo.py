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
This module contains a zero-order representation of an ion exchange unit
operation.
"""

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
