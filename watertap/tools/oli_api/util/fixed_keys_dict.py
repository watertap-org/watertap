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
        elif isinstance(self.data[k]["value"], list):
            if v not in self.data[k]:
                raise RuntimeError(
                    f" Value {v} not a valid value for key {k}."
                )            
        else:
            self.data[k] = v

    def __delitem__(self, k):
        raise Exception(" Deleting keys not supported for this object.")

    def pprint(self):
        print("-------------------------")
        for key, value in self.data.items():
            print(f" {key}\n - {value}\n")

default_water_analysis_properties = FixedKeysDict(
    {   
        "temperature": {"group": "Properties", "name": "Temperature", "unit": "K", "value": None},
        "pressure": {"group": "Properties", "name": "Pressure", "unit": "Pa", "value": None},
        "component": {"group": None, "name": None, "unit": "mg/L", "value": None, "charge": None},
        "electroneutrality": {"group": "Electroneutrality Options", "name": "ElectroneutralityBalanceType", "value": ["DominantIon", "ProrateCations", "ProrateAnions", "Prorate", "AutoNACL", "MakeupIon"]},
        "makeup_ion": {"group": "Electroneutrality Options", "name": "MakeupIonBaseTag", "value": None},
        "reconciliation": {"group": "Calculation Options", "name": "CalcType", "value": ["EquilCalcOnly", "ReconcilePh", "ReconcilePhAndAlkalinity", "ReconcilePhAndAlkalinityAndTic", "ReconcileCo2Gas"]},
        "ph": {"group": "Properties", "name": "pH", "value": None},
        "acid_titrant": {"group": "Calculation Options", "name": "PhAcidTitrant", "value": None},
        "base_titrant": {"group": "Calculation Options", "name": "PhBaseTitrant", "value": None},
        "alkalinity": {"group": "Properties", "name": "Alkalinity", "unit": "mg HCO3/L", "value": None},
        "alk_titrant": {"group": "Calculation Options", "name": "AlkalinityPhTitrant", "value": None},
        "AlkalinityTitrationEndPointPh": {"group": "Properties", "name": "AlkalinityTitrationEndPointpH", "value": None},        
        "tic": {"group": "Properties", "name": "TIC", "unit": "mol C/L", "value": None},
        "co2_fraction":  {"group": "Properties", "name": "CO2GasFraction", "unit": "mole %", "value": None},
        "solids": {"group":"Calculation Options", "name": "AllowSolidsToForm", "value": [False, True]},
        "calc_alk": {"group":"Calculation Options", "name": "CalcAlkalnity", "value": [False, True]},
    }
)

# Full list of available optional inputs available: https://devdocs.olisystems.com/optional-inputs
default_optional_properties = FixedKeysDict(
    {
        "scalingIndex": False,
        "prescalingTendencies": False,
        "prescalingTendenciesRigorous": False,
        "scalingTendencies": False,
        "MBGComposition": True,
        "materialBalanceGroup": True,
    }
)

default_unit_set_info = FixedKeysDict(
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