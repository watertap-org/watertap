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
This module contains a zero-order representation of a peracetic acid (PAA) water
disinfection unit.
"""

from pyomo.environ import Var, Suffix, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData


# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("PeraceticAcidDisinfectionZO")
class PeraceticAcidDisinfectionData(ZeroOrderBaseData):
    """
    Zero-Order model for a peracetic acid water disinfection unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "peracetic_acid_disinfection"

        build_sido_reactive(self)
        constant_intensity(self)

        # Create hydraulic retention time variable
        self.HRT = Var(
            units=pyunits.hr,
            bounds=(0, None),
            doc="Hydraulic retention time of water treatment unit",
        )
        self._perf_var_dict["Hydraulic Retention Time"] = self.HRT
        self._fixed_perf_vars.append(self.HRT)

        # Create variable for mass of an E. coli cell
        self.ecoli_cell_mass = Var(
            units=pyunits.kg,
            bounds=(0, None),
            doc="Average mass of an E. coli cell",
        )
        self._perf_var_dict["E. coli Average Cell Mass"] = self.ecoli_cell_mass
        self._fixed_perf_vars.append(self.ecoli_cell_mass)

        # Create variable for weight fraction of PAA in disinfection solution
        self.disinfection_solution_wt_frac_PAA = Var(
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Weight fraction of PAA in disinfection solution",
        )
        self._perf_var_dict[
            "Weight fraction PAA in disinfection solution"
        ] = self.disinfection_solution_wt_frac_PAA
        self._fixed_perf_vars.append(self.disinfection_solution_wt_frac_PAA)

        # Create variable for disinfection solution density
        self.disinfection_solution_density = Var(
            units=pyunits.kg / pyunits.liter,
            bounds=(0, None),
            doc="Disinfection solution density",
        )
        self._perf_var_dict[
            "disinfection solution density"
        ] = self.disinfection_solution_density
        self._fixed_perf_vars.append(self.disinfection_solution_density)

        # Create variable for disinfection solution volumetric flowrate
        self.disinfection_solution_flow_vol = Var(
            self.flowsheet().time,
            units=pyunits.L / pyunits.s,
            bounds=(0, None),
            doc="Volumetric flowrate of disinfection solution",
        )
        self._perf_var_dict[
            "Disinfection solution volumetric flowrate"
        ] = self.disinfection_solution_flow_vol

        # Create constraint to calculate disinfection solution flowrate
        @self.Constraint(
            self.flowsheet().time, doc="Constraint for disinfection solution flowrate"
        )
        def disinfection_solution_flow_vol_rule(b, t):
            return pyunits.convert(
                b.inlet.flow_mass_comp[t, "peracetic_acid"],
                to_units=pyunits.kg / pyunits.s,
            ) == pyunits.convert(
                b.disinfection_solution_flow_vol[t]
                * b.disinfection_solution_density
                * b.disinfection_solution_wt_frac_PAA,
                to_units=pyunits.kg / pyunits.s,
            )

        # Create reactor volume variable
        self.reactor_volume = Var(
            units=pyunits.m**3,
            bounds=(0, None),
            doc="Volume of water treatment unit",
        )
        self._perf_var_dict["Reactor Volume"] = self.reactor_volume

        # Create constraint relating HRT, volume, and volumetric flowrate
        @self.Constraint(self.flowsheet().time, doc="Constraint for reactor volume")
        def reactor_volume_rule(b, t):
            return b.reactor_volume == (
                pyunits.convert(
                    b.HRT * b.properties_in[t].flow_vol, to_units=pyunits.m**3
                )
            )

        # Create variable for E. coli concentration at reactor inlet
        self.inlet_ecoli_conc = Var(
            self.flowsheet().time,
            units=pyunits.liter**-1,
            bounds=(0, None),
            doc="Concentration of E. coli at reactor inlet",
        )
        self._perf_var_dict["Inlet E. coli Concentration"] = self.inlet_ecoli_conc

        # Create constraint relating E. coli inlet mass flow rate and
        # concentration
        @self.Constraint(
            self.flowsheet().time, doc="Constraint for E. coli inlet concentration"
        )
        def ecoli_inlet_concentration_rule(b, t):
            return pyunits.convert(
                b.inlet.flow_mass_comp[t, "total_coliforms_fecal_ecoli"],
                to_units=pyunits.kg / pyunits.s,
            ) == pyunits.convert(
                b.inlet_ecoli_conc[t] * b.ecoli_cell_mass * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        # Create variable for E. coli concentration at reactor outlet
        self.outlet_ecoli_conc = Var(
            self.flowsheet().time,
            units=pyunits.liter**-1,
            bounds=(0, None),
            doc="Concentration of E. coli at reactor outlet",
        )
        self._perf_var_dict["Outlet E. coli Concentration"] = self.outlet_ecoli_conc

        # Create constraint relating E. coli outlet mass flow rate and
        # concentration
        @self.Constraint(
            self.flowsheet().time, doc="Constraint for E. coli outlet concentration"
        )
        def ecoli_outlet_concentration_rule(b, t):
            return pyunits.convert(
                b.treated.flow_mass_comp[t, "total_coliforms_fecal_ecoli"],
                to_units=pyunits.kg / pyunits.s,
            ) == pyunits.convert(
                b.outlet_ecoli_conc[t]
                * b.ecoli_cell_mass
                * b.properties_treated[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )
