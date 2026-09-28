#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
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

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Optional

import numpy as np


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
    feed_cost: float = None  # in $/m3


@dataclass
class UFParams(UnitParams):
    """Parameters for the UF unit"""

    num_uf_pumps: int = 4
    minimum_operating_pumps: int = 1
    allow_shutdown: bool = True
    minimum_flowrate: float = 344
    nominal_flowrate: float = 900
    maximum_flowrate: float = 989
    nominal_recovery: float = 1
    minimum_uptime: int = 2
    minimum_downtime: int = 2
    startup_delay: int = 1
    allow_variable_recovery: bool = False
    chemical_cost: float = None  # in $/m3

    def __post_init__(self):
        # self._surrogate = # load the surrogate model here.
        self.surrogate_type: str = "quadratic_energy_intensity"
        self.surrogate_a = 1
        self.surrogate_b = 1
        self.surrogate_c = 1

    @property
    def surrogate_coeffs(self):
        """Returs the coefficients of the surrogate model as a dictionary"""
        return {
            "a": self.surrogate_a,
            "b": self.surrogate_b,
            "c": self.surrogate_c,
        }


@dataclass
class ROParams(UnitParams):
    """Parameters for the RO unit"""

    num_ro_skids: int = 4
    minimum_operating_skids: int = 2
    allow_shutdown: bool = True
    minimum_flowrate: float = 0
    nominal_flowrate: float = 337.670
    maximum_flowrate: float = 400
    minimum_recovery: float = 0.88
    nominal_recovery: float = 0.92
    maximum_recovery: float = 0.925
    minimum_uptime: int = 2
    minimum_downtime: int = 2
    max_num_skids_shutdown_per_timestep: int = 2
    startup_delay: int = 1
    allow_variable_recovery: bool = False
    replacement_types: list[str] = field(default_factory=list)
    replacement_costs: list[float] = field(default_factory=list)
    replacement_lifetimes: list[float] = field(default_factory=list)
    replacement_max_flex_penalty: list[float] = field(default_factory=list)

    def __post_init__(self):
        # self._surrogate = # load the surrogate model here.
        self.surrogate_type: str = "constant_energy_intensity"
        self.surrogate_file: Optional[str] = None
        self.surrogate_a = 1
        self.surrogate_b = 1
        self.surrogate_c = 1

    @property
    def surrogate_coeffs(self):
        """Returs the coefficients of the surrogate model as a dictionary"""
        return {
            "a": self.surrogate_a,
            "b": self.surrogate_b,
            "c": self.surrogate_c,
        }


@dataclass
class PosttreatmentParams(UnitParams):
    """Parameters for the posttreatment unit"""

    energy_intensity: float = 0.41
    leakage_fraction: float = 0
    chemical_cost: float = None  # in $/m3


@dataclass
class BrineDischargeParams(UnitParams):
    """Parameters for the brine discharge unit"""

    energy_intensity: float = 0.1
    leakage_fraction: float = 0
    brine_cost: float = None  # in $/m3


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
    """Parameters for flexible desalination"""

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
    # Other parameters not used in tutorial, but have related functions in flex_recovery_and_flow_flowsheet.py
    nonworking_hours: list[int] = field(default_factory=list)
    CAPEX_yr: float = None
    max_daily_shutdowns: Optional[int] = None

    def __post_init__(self):
        self.intake = IntakeParams()
        self.uf = UFParams()
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
