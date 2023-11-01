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
#
#################################################################################

__author__ = "Paul Vecchiarelli, Ben Knueven"

from collections import UserDict


class FixedKeysDict(UserDict):
    def __init__(self, d):
        self.data = d

    def __setitem__(self, k, v):
        if k not in self.data:
            raise RuntimeError(f" Key {k} not in dictionary.")
        else:
            self.data[k] = v

    def __delitem__(self, k):
        raise Exception(" Deleting keys not supported for this object.")

    def _check_value(self, k, valid_values):
        if k not in self.data:
            raise RuntimeError(f" Key {k} not in dictionary.")
        else:
            if self.data[k] not in valid_values:
                raise RuntimeError(
                    f" Value {self.data[k]} not a valid value for key {k}."
                )

    def pprint(self):
        print("-------------------------")
        for key, value in self.data.items():
            print(f" {key}\n - {value}\n")


default_oli_input_dict = FixedKeysDict(
    {
        "temperature_unit": "K",  # units
        "pressure_unit": "Pa",
        "concentration_unit": "mg/L",
        "gas_fraction_unit": "mole %",
        "alkalinity_unit": "mg HCO3/L",
        "tic_unit": "mol C/L",
        # essential input parameters for water analysis
        "electroneutrality_value": "DominantIon",
        "reconciliation_value": "EquilCalcOnly",
        "AllowSolidsToForm": False,
        "CalcAlkalnity": False,
        # optional inputs
        # required for certain specifications
        "MakeupIonBaseTag": None,
        "CO2GasFraction": None,
        "pH": None,
        "PhAcidTitrant": "HCL",
        "PhBaseTitrant": "NAOH",
        "Alkalinity": None,
        "AlkalinityPhTitrant": None,
        "AlkalinityTitrationEndPointPh": None,
        "TIC": None,
    }
)

# Full list of available optional inputs available: https://devdocs.olisystems.com/optional-inputs
default_oli_optional_properties = FixedKeysDict(
    {
        "scalingIndex": False,
        "prescalingTendencies": False,
        "prescalingTendenciesRigorous": False,
        "scalingTendencies": False,
        "MBGComposition": True,
        "materialBalanceGroup": True,
    }
)

default_oli_unit_set_info = FixedKeysDict(
    {
        "enthalpy": "J",
        "mass": "kg",
        "pt": "Pa",
        "total": "mg",
        "liq1_phs_comp": "mg",
        "solid_phs_comp": "mg",
        "vapor_phs_comp": "mg",
        "liq2_phs_comp": "mg",
        "combined_phs_comp": "mg",
        "molecularConcentration": "mg/L",
    }
)
