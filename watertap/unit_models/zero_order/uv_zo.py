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
This module contains a zero-order representation of a UV reactor unit
operation.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.common.config import ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVZO")
class UVZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a UV unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "uv"

        build_siso(self)
        constant_intensity(self)

        self.uv_reduced_equivalent_dose = Var(
            self.flowsheet().time,
            units=pyunits.mJ / pyunits.cm**2,
            doc="Reduced equivalent dosage",
        )
        self.uv_transmittance_in = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="UV transmittance of solution at UV reactor inlet",
        )

        self.recovery_frac_mass_H2O.fix(1)
        self._fixed_perf_vars.append(self.uv_reduced_equivalent_dose)
        self._fixed_perf_vars.append(self.uv_transmittance_in)

        self._perf_var_dict[
            "UV Reduced Equivalent Dosage (mJ/cm^2)"
        ] = self.uv_reduced_equivalent_dose
        self._perf_var_dict["UV Transmittance of Feed"] = self.uv_transmittance_in
