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

__author__ = "Paul Vecchiarelli, Ben Knueven, Adam Atia"

from collections import UserDict
from pyomo.environ import units as pyunits
from pyomo.core.base.units_container import _PyomoUnit
from collections.abc import Iterable


class FixedKeysDict(UserDict):
    def __init__(self, d):
        self.data = d

    def __setitem__(self, k, v):
        if k not in self.data:
            raise RuntimeError(f" Key {k} not in dictionary.")
        # also check for valid value if a list of values is given in the default dictionary
        else:
            # if user setting value in pyomo units
            if hasattr(v, "is_expression_type"):
                if isinstance(v, _PyomoUnit) or v.is_expression_type():
                    # if user assigns pyomo units as value to oli_unit, save the str to oli_unit and update pyomo_unit with pyomo units
                    if "oli_unit" in k:
                        self.data[k] = str(v)
                        self.data["pyomo_unit"] = v
                    # if user assigns pyomo units to pyomo_unit, update oli_unit with str representation of units
                    if "pyomo_unit" in k:
                        self.data[k] = v
                        self.data["oli_unit"] = str(v)
            # check if data[k] is iterable first, otherwise checking if oli_unit in data[k] throws exception
            elif isinstance(self.data[k], Iterable):
                # check if user provides str and that assignment wouldn't overwrite the oli_unit key:value pair
                if isinstance(v, str) and ("oli_unit" not in self.data[k]):
                    # if user assigns str to oli_unit, update pyomo_units with PyomoUnits representation of str
                    if "oli_unit" in k:
                        self.data[k] = v
                        self.data["pyomo_unit"] = getattr(pyunits, v)
                    else:
                        pass
                elif isinstance(v, str) and ("oli_unit" in self.data[k]):
                    self.data["oli_unit"] = v
                    self.data["pyomo_unit"] = getattr(pyunits, v)
                elif not isinstance(v, str):
                    raise RuntimeError(f"Setting {v} as the value for {k} is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units.")
                else:
                    pass
            elif not isinstance(self.data[k], Iterable):
                if isinstance(v, str) and "pyomo_unit" in k:
                    self.data[k] = getattr(pyunits, v)
                    self.data["oli_unit"] = v
                elif not isinstance(v, str):
                    raise RuntimeError(f"Setting {v} as the value for {k} is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units.")
                else:
                    pass
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


input_unit_set_temp = {
        "molecularConcentration": {
            "oli_unit": "mg/L",
            "pyomo_unit": pyunits.mg / pyunits.L,
        },
        "mass": {"oli_unit": "mg", "pyomo_unit": pyunits.mg},
        "temperature": {"oli_unit": "K", "pyomo_unit": pyunits.K},
        "pressure": {"oli_unit": "Pa", "pyomo_unit": pyunits.Pa},
        "enthalpy": {"oli_unit": "J", "pyomo_unit": pyunits.J},
        "vaporAmountMoles": {"oli_unit": "mol", "pyomo_unit": pyunits.mol},
        "vaporMolFrac": {
            "oli_unit": "mol/mol",
            "pyomo_unit": pyunits.mol / pyunits.mol,
        },
        "totalVolume": {"oli_unit": "L", "pyomo_unit": pyunits.L},
        "pipeDiameter": {"oli_unit": "m", "pyomo_unit": pyunits.meter},
        "pipeFlowVelocity": {
            "oli_unit": "m/s",
            "pyomo_unit": pyunits.meter / pyunits.second,
        },
        "diskDiameter": {"oli_unit": "m", "pyomo_unit": pyunits.meter},
        "diskRotatingSpeed": {"oli_unit": "cycle/s", "pyomo_unit": 1 / pyunits.second},
        "rotorDiameter": {"oli_unit": "m", "pyomo_unit": pyunits.meter},
        "rotorRotation": {"oli_unit": "cycle/s", "pyomo_unit": 1 / pyunits.second},
        "shearStress": {"oli_unit": "Pa", "pyomo_unit": pyunits.Pa},
        "pipeDiameter": {"oli_unit": "m", "pyomo_unit": pyunits.meter},
        "pipeRoughness": {"oli_unit": "m", "pyomo_unit": pyunits.meter},
        "liquidFlowInPipe": {
            "oli_unit": "L/s",
            "pyomo_unit": pyunits.L / pyunits.second,
        },
        "gasFlowInPipe": {"oli_unit": "L/s", "pyomo_unit": pyunits.L / pyunits.second},
        "viscAbs2ndLiq": {
            "oli_unit": "Pa-s",
            "pyomo_unit": pyunits.Pa * pyunits.second,
        },
        "alkalinity": {"oli_unit": "mg HCO3/L", "pyomo_unit": pyunits.mg / pyunits.L},
        "TIC": {"oli_unit": "mol C/L", "pyomo_unit": pyunits.mol / pyunits.L},
        "CO2GasFraction": {
            "oli_unit": "mol/mol",
            "pyomo_unit": pyunits.mol / pyunits.mol,
        },
    }

input_unit_set = FixedKeysDict({k:FixedKeysDict(v) for k,v in input_unit_set_temp.items()})
default_unit_set = FixedKeysDict({k:FixedKeysDict(v) for k,v in input_unit_set_temp.items()})

optional_properties = FixedKeysDict(
    {
        "electricalConductivity": True,
        "viscosity": True,
        "selfDiffusivityAndMobility": True,
        "heatCapacity": True,
        "thermalConductivity": True,
        "surfaceTension": True,
        "interfacialTension": True,
        "volumeStdConditions": True,
        "prescalingTendenciesEstimated": False,
        "prescalingIndexEstimated": False,
        "prescalingTendenciesRigorous": True,
        "prescalingIndexRigorous": True,
        "scalingTendencies": True,
        "scalingIndex": True,
        "hardness": True,
        "ionicStrengthXBased": True,
        "ionicStrengthMBased": True,
        "totalDissolvedSolids": True,
        "vaporToInflowMoleFraction": True,
        "partialPressure": True,
        "vaporDiffusivityMatrix": True,
        "entropyStream": True,
        "entropySpecies": True,
        "entropyStreamStandardState": True,
        "entropySpeciesStandardState": True,
        "gibbsEnergyStream": True,
        "gibbsEnergySpecies": True,
        "gibbsEnergyStreamStandardState": True,
        "gibbsEnergySpeciesStandardState": True,
        "activityCoefficientsXBased": True,
        "activityCoefficientsMBased": True,
        "fugacityCoefficients": True,
        "vaporFugacity": True,
        "kValuesXBased": True,
        "kValuesMBased": True,
        "MBGComposition": True,
        "materialBalanceGroup": True,
    }
)

# TODO: consider adding these: https://devdocs.olisystems.com/user-defined-output-unit-set
# and reducing hard-coding by using default_input_unit_set references
output_unit_set = FixedKeysDict(
    {
        "enthalpy": default_unit_set["enthalpy"]["oli_unit"],
        "mass": default_unit_set["mass"]["oli_unit"],
        "pt": default_unit_set["pressure"]["oli_unit"],
        "total": default_unit_set["mass"]["oli_unit"],
        "liq1_phs_comp": default_unit_set["mass"]["oli_unit"],
        "solid_phs_comp": default_unit_set["mass"]["oli_unit"],
        "vapor_phs_comp": default_unit_set["mass"]["oli_unit"],
        "liq2_phs_comp": default_unit_set["mass"]["oli_unit"],
        "combined_phs_comp": default_unit_set["mass"]["oli_unit"],
        "molecularConcentration": default_unit_set["molecularConcentration"]["oli_unit"],
    }
)

if __name__ == "__main__":
    unit_set=input_unit_set_temp