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
General costing package for zero-order processes.
"""
import os
import yaml

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import register_idaes_currency_units
from watertap.costing.watertap_costing_package import WaterTAPCostingDetailedData

# NOTE: some of these are defined in WaterTAPCostingBlockData
global_params = [
    "plant_lifetime",
    "utilization_factor",
    "land_cost_percent_FCI",
    "working_capital_percent_FCI",
    "salaries_percent_FCI",
    "benefit_percent_of_salary",
    "maintenance_costs_percent_FCI",
    "laboratory_fees_percent_FCI",
    "insurance_and_taxes_percent_FCI",
    "wacc",
    "TPEC",
    "TIC",
]


@declare_process_block_class("ZeroOrderCosting")
class ZeroOrderCostingData(WaterTAPCostingDetailedData):
    """
    General costing package for zero-order processes.
    """

    CONFIG = WaterTAPCostingDetailedData.CONFIG()
    CONFIG.declare(
        "case_study_definition",
        ConfigValue(
            default=None,
            doc="Path to YAML file defining global parameters for case study. If "
            "not provided, default values from the WaterTap database are used.",
        ),
    )

    def register_currency_definitions(self):
        # Register currency and conversion rates
        if "currency_definitions" in self._cs_def:
            pyo.units.load_definitions_from_strings(
                self._cs_def["currency_definitions"]
            )
        else:
            register_idaes_currency_units()

    def set_base_currency_base_period(self):
        # Set the base year for all costs
        self.base_currency = getattr(pyo.units, self._cs_def["base_currency"])
        # Set a base period for all operating costs
        self.base_period = getattr(pyo.units, self._cs_def["base_period"])

    def build_global_params(self):
        """
        To minimize overhead, only create global parameters for now.
        Unit-specific parameters will be added as sub-Blocks on a case-by-case
        basis as a unit of that type is costed.
        """
        # Load case study definition from file
        self._cs_def = _load_case_study_definition(self)

        super().build_global_params()

        # Define expected flows
        for f, v in self._cs_def["defined_flows"].items():
            value = v["value"]
            units = getattr(pyo.units, v["units"])
            if self.component(f + "_cost") is not None:
                self.component(f + "_cost").fix(value * units)
            else:
                self.defined_flows[f] = value * units

        # Fix all Vars from database
        for v in global_params:
            try:
                value = self._cs_def["global_parameters"][v]["value"]
                units = self._cs_def["global_parameters"][v]["units"]
                getattr(self, v).fix(value * getattr(pyo.units, units))
            except KeyError:
                raise KeyError(
                    f"Invalid case study definition file - no entry found "
                    f"for {v}, or entry lacks value and units."
                )
        # this variable floats in this package, so initialize it appropriately here
        calculate_variable_from_constraint(
            self.capital_recovery_factor, self.factor_capital_annualization_constraint
        )


def _load_case_study_definition(self):
    """
    Load data from case study definition file into a Python dict.
    If users did not provide a definition file as a config argument, the
    default definition from the WaterTap techno-economic database is used.
    """
    source_file = self.config.case_study_definition
    if source_file is None:
        source_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "techno_economic",
            "default_case_study.yaml",
        )

    try:
        with open(source_file, "r") as f:
            lines = f.read()
            f.close()
    except OSError:
        raise OSError(
            "Could not find specified case study definition file. "
            "Please check the path provided."
        )

    return yaml.load(lines, yaml.Loader)
