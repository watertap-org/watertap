#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import Constraint

from idaes.core import declare_process_block_class, EnergyBalanceType
from idaes.core.base.control_volume_base import ControlVolumeBlockData as _CV
from idaes.core.base.control_volume0d import ControlVolume0DBlockData as _CV0D
from idaes.core.base.control_volume1d import ControlVolume1DBlockData as _CV1D
from idaes.core.util.exceptions import ConfigurationError

import idaes.core.util.scaling as iscale


class _IsothermalEnergyBalanceCheckerMixin:
    def _validate_isothermal_energy_balance_compatibility(self):
        if self._constructed_energy_balance_type is not EnergyBalanceType.none:
            raise ConfigurationError(
                "When using the isothermal assumption, the EnergyBalanceType must be EnergyBalanceType.none"
            )

    def add_energy_balances(self, balance_type=EnergyBalanceType.useDefault, **kwargs):
        super().add_energy_balances(balance_type=balance_type, **kwargs)
        if hasattr(self, "isothermal_assumption_eq"):
            self._validate_isothermal_energy_balance_compatibility()

    def add_isothermal_assumption(self):
        if hasattr(self, "_constructed_energy_balance_type"):
            self._validate_isothermal_energy_balance_compatibility()


# put the usual docstring on add_energy_balances
_IsothermalEnergyBalanceCheckerMixin.add_energy_balances.__doc__ = (
    _CV.add_energy_balances.__doc__
)


@declare_process_block_class("ControlVolume0DBlock")
class ControlVolume0DBlockData(_IsothermalEnergyBalanceCheckerMixin, _CV0D):
    def add_isothermal_assumption(self):
        super().add_isothermal_assumption()

        @self.Constraint(
            self.flowsheet().time,
            doc="Isothermal assumption for 0D control volume",
        )
        def isothermal_assumption_eq(b, t):
            return b.properties_out[t].temperature == b.properties_in[t].temperature

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        if hasattr(self, "isothermal_assumption_eq"):
            for t, c in self.isothermal_assumption_eq.items():
                sf = iscale.get_scaling_factor(
                    self.properties_in[t].temperature, default=1e-2, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)


@declare_process_block_class("ControlVolume1DBlock")
class ControlVolume1DBlockData(_IsothermalEnergyBalanceCheckerMixin, _CV1D):
    def add_isothermal_assumption(self):
        super().add_isothermal_assumption()

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            doc="Isothermal assumption for 1D control volume",
        )
        def isothermal_assumption_eq(b, t, x):
            if x == b.length_domain.first():
                return Constraint.Skip
            return (
                b.properties[t, b.length_domain.first()].temperature
                == b.properties[t, x].temperature
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        if hasattr(self, "isothermal_assumption_eq"):
            for (t, x), c in self.isothermal_assumption_eq.items():
                sf = iscale.get_scaling_factor(
                    self.properties[t, x].temperature, default=1e-2, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)
