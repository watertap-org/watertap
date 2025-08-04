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

"""
This module contains the default values of all the required
parameters.
"""

from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Optional

import numpy as np
from pyomo.environ import exp
from scipy import optimize


@dataclass
class UnitParams:
    """Abstract dataclass for parameters of all units"""

    energy_intensity: Optional[float] = None
    allow_shutdown: bool = False
    leakage_fraction: Optional[float] = None
    minimum_flowrate: Optional[float] = None
    nominal_flowrate: Optional[float] = None
    maximum_flowrate: Optional[float] = None
    nominal_recovery: Optional[float] = None
    minimum_recovery: Optional[float] = None
    maximum_recovery: Optional[float] = None
    minimum_uptime: Optional[int] = None
    minimum_downtime: Optional[int] = None
    startup_delay: Optional[int] = None

    @property
    def get_leakage_fraction(self):
        """Returns the leakage fraction value"""
        if self.leakage_fraction is not None:
            return self.leakage_fraction

        if self.recovery is not None:
            return 1 - self.recovery

        raise ValueError("leakage_fraction is not specified")

    @property
    def get_recovery(self):
        """Returns the recovery value"""
        if self.nominal_recovery is not None:
            return self.nominal_recovery

        if self.leakage_fraction is not None:
            return 1 - self.leakage_fraction

        raise ValueError("recovery is not specified")

    def update(self, value_map: dict):
        """Updates the values of the specified attributes"""
        # Raise an error if an unrecognized attribute is provided
        new_values = value_map.copy()
        for key in value_map:
            if key not in self.__dict__:
                if "_" + key in self.__dict__:
                    # This is a property
                    setattr(self, key, value_map[key])
                    # Remove this element
                    new_values.pop(key)

                else:
                    raise KeyError(f"Unrecognized attribute {key}")

        self.__dict__.update(new_values)


@dataclass
class IntakeParams(UnitParams):
    """Parameters for the intake unit"""

    energy_intensity: float = 0.157121734
    leakage_fraction: float = 0
    minimum_flowrate: float = 1063.5
    nominal_flowrate: float = 1063.5
    maximum_flowrate: float = 1063.5


@dataclass
class PretreatmentParams(UnitParams):
    """Parameters for the pretreatment unit"""

    allow_shutdown: bool = False
    energy_intensity: float = 0.01
    leakage_fraction: float = 0
    minimum_downtime: int = 0
    startup_delay: int = 0


@dataclass
class ROParams(UnitParams):
    """Parameters for the RO unit"""

    num_ro_skids: int = 3
    minimum_operating_skids: int = 2
    allow_shutdown: bool = True
    nominal_flowrate: float = 337.670
    minimum_recovery: float = 0.4
    nominal_recovery: float = 0.465
    maximum_recovery: float = 0.55
    minimum_uptime: int = 1
    minimum_downtime: int = 4
    startup_delay: int = 8
    allow_variable_recovery: bool = False

    def __post_init__(self):
        self._surrogate_type = "exponential_quadratic"
        self.surrogate_a = 6.180228375549232
        self.surrogate_b = 2.3684824891480476
        self.surrogate_c = 6.474944185881354
        self.surrogate_d = 1.9065595663669615e-05

    @property
    def surrogate_type(self):
        """Returns the surrogate type for RO energy intensity"""
        return self._surrogate_type

    @surrogate_type.setter
    def surrogate_type(self, value: str):
        if value in ["exponential_quadratic", "quadratic_surrogate"]:
            self._surrogate_type = value

        else:
            raise ValueError("Unrecognized surrogate type")

    @property
    def surrogate_coeffs(self):
        """Returs the coefficients of the surrogate model as a dictionary"""
        return {
            "a": self.surrogate_a,
            "b": self.surrogate_b,
            "c": self.surrogate_c,
            "d": self.surrogate_d,
        }

    def get_energy_intensity(self, recovery):
        """Returns the energy intensity for a given recovery"""
        coeffs = self.surrogate_coeffs
        if self.surrogate_type == "exponential_quadratic":
            return (
                coeffs["a"] * exp(-coeffs["b"] * recovery)
                + coeffs["c"] * recovery**2
                + coeffs["d"]
            )

        if self.surrogate_type == "quadratic_surrogate":
            return coeffs["a"] * recovery**2 + coeffs["b"] * recovery + coeffs["c"]

        return None

    def get_optimum_energy_intensity(self, recovery_lb, recovery_ub):
        """
        Returns the optimum energy intensity if it exists inside the
        interval. Returns None is the optimum is at the bounds.
        """
        # Optimum exists inside the interval, and it is unique.
        coeffs = self.surrogate_coeffs
        if self.surrogate_type == "exponenetial_quadratic":

            def _first_der(rec):
                return (
                    -coeffs["a"] * coeffs["b"] * exp(-coeffs["b"] * rec)
                    + 2 * coeffs["c"] * rec
                )

            # First derivative is a monotonically increasing function.
            # Therefore, if the first derivative has the same sign at both
            # ends of the interval, then there is no point in the interval
            # at which the derivative vanishes. So, the optimum is at its bounds
            if _first_der(recovery_lb) * _first_der(recovery_ub) > 0:
                # Optimum does not exist, so return
                return None

            root = optimize.bisect(_first_der, recovery_lb, recovery_ub, maxiter=1000)

        else:
            # This is "quadratic_surrogate":
            root = -coeffs["b"] / (2 * coeffs["a"])

        return self.get_energy_intensity(root)

    def get_energy_intensity_bounds(self, recovery_lb=None, recovery_ub=None):
        """
        Returns the bounds on energy intensity based on the bounds of
        recovery
        """
        if recovery_lb is None:
            recovery_lb = self.minimum_recovery
        if recovery_ub is None:
            recovery_ub = self.maximum_recovery

        ei_values = [
            self.get_energy_intensity(recovery=self.minimum_recovery),
            self.get_energy_intensity(recovery=self.maximum_recovery),
            self.get_optimum_energy_intensity(recovery_lb, recovery_ub),
        ]

        ei_values = list(filter(None, ei_values))  # remove None, if it exists
        return min(ei_values), max(ei_values)


@dataclass
class PosttreatmentParams(UnitParams):
    """Parameters for the posttreatment unit"""

    energy_intensity: float = 0.41
    leakage_fraction: float = 0


@dataclass
class BrineDischargeParams(UnitParams):
    """Parameters for the brine discharge unit"""

    energy_intensity: float = 0.1
    leakage_fraction: float = 0


@dataclass
class Battery:
    """Parameters for the battery"""

    energy_capacity: float = 0
    power_capacity: float = 50
    efficiency: float = 0.86
    initial_soc: float = 0.5
    minimum_soc: float = 0.2
    maximum_soc: float = 0.95


@dataclass
class FlexDesalParams:
    """Parameters for fleible desalination"""

    start_date: str = "2022-07-05 00:00:00"
    end_date: str = "2022-07-06 00:00:00"
    timestep_hours: float = 0.25

    product_water_price: float = 0
    fixed_monthly_cost: float = 766000
    customer_rate: float = 100
    constrain_to_baseline_production: bool = False
    curtailment_fraction: float = 0.0
    annual_production_AF: float = 3125  # in acre-ft / year
    production_constraint_to_objective: bool = False
    production_constraint_penalty: float = 0.6
    emissions_cost: float = 0  # Cost of emissions in $/kg

    include_demand_response: bool = False
    include_battery: bool = False
    include_onsite_solar: bool = False
    onsite_capacity: float = 0

    def __post_init__(self):
        self.intake = IntakeParams()
        self.pretreatment = PretreatmentParams()
        self.ro = ROParams()
        self.posttreatment = PosttreatmentParams()
        self.brinedischarge = BrineDischargeParams()
        self.battery = Battery()

        # datetime array
        t = np.arange(
            datetime.fromisoformat(self.start_date),
            datetime.fromisoformat(self.end_date),
            timedelta(hours=self.timestep_hours),
        ).astype(datetime)

        # length of time step in seconds
        dt_seconds = self.timestep_hours * 3600
        total_num_seconds = (t[-1] - t[0]).total_seconds() + dt_seconds

        self.num_hours = total_num_seconds / 3600
        self.num_days = self.num_hours / 24
        self.num_months = self.num_days / 31
