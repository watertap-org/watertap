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
This module contains a zero-order representation of a reactive anaerobic
digestion unit operation.
"""

from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData


# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("AnaerobicDigestionReactiveZO")
class AnaerobicDigestionReactiveZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a reactive anaerobic digestion unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "anaerobic_digestion_reactive"

        build_sido_reactive(self)
        constant_intensity(self)

        if (
            self.config.process_subtype == "default"
            or self.config.process_subtype is None
        ):

            # Create water flux variable
            self.biogas_tss_ratio = Var(
                units=pyunits.m**3 / pyunits.kg,
                bounds=(0, None),
                doc="Ratio of m^3 biogas produced / kg TSS in influent",
            )
            self._perf_var_dict[
                "Ratio of m^3 biogas produced / kg TSS in influent"
            ] = self.biogas_tss_ratio
            self._fixed_perf_vars.append(self.biogas_tss_ratio)

            self.biogas_production = Var(
                self.flowsheet().time,
                units=pyunits.m**3 / pyunits.s,
                bounds=(0, None),
                doc="Biogas production",
            )
            self._perf_var_dict["Biogas production"] = self.biogas_production

            @self.Constraint(
                self.flowsheet().time, doc="Constraint for biogas production"
            )
            def biogas_prod(b, t):
                return b.biogas_production[t] == (
                    pyunits.convert(
                        b.biogas_tss_ratio * b.properties_in[t].flow_mass_comp["tss"],
                        to_units=pyunits.m**3 / pyunits.s,
                    )
                )

        if self.config.process_subtype == "GLSD_anaerobic_digester":
            self.HRT = Var(
                units=pyunits.hour,
                doc="Hydraulic retention time",
            )

            self._fixed_perf_vars.append(self.HRT)
            self._perf_var_dict["Hydraulic retention time"] = self.HRT

            self.reactor_volume = Var(
                self.flowsheet().time,
                units=pyunits.m**3,
                doc="Reactor volume",
            )

            @self.Constraint(self.flowsheet().time, doc="Constraint for reactor volume")
            def reactor_volume_cons(b, t):
                return b.reactor_volume[t] == (
                    pyunits.convert(
                        b.HRT * b.properties_in[t].flow_vol,
                        to_units=pyunits.m**3,
                    )
                )

            self._perf_var_dict["Reactor volume"] = self.reactor_volume
