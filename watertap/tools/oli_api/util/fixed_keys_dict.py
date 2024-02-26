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
from pyomo.environ import units as pyunits


class FixedKeysDict(UserDict):
    def __init__(self, d):
        self.data = d

    def __setitem__(self, k, v):
        if k not in self.data:
            raise RuntimeError(f" Key {k} not in dictionary.")
        # also check for valid value if a list of values is given in the default dictionary
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


input_unit_set = FixedKeysDict(
    {
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
)

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
        "enthalpy": input_unit_set["enthalpy"]["oli_unit"],
        "mass": input_unit_set["mass"]["oli_unit"],
        "pt": input_unit_set["pressure"]["oli_unit"],
        "total": input_unit_set["mass"]["oli_unit"],
        "liq1_phs_comp": input_unit_set["mass"]["oli_unit"],
        "solid_phs_comp": input_unit_set["mass"]["oli_unit"],
        "vapor_phs_comp": input_unit_set["mass"]["oli_unit"],
        "liq2_phs_comp": input_unit_set["mass"]["oli_unit"],
        "combined_phs_comp": input_unit_set["mass"]["oli_unit"],
        "molecularConcentration": input_unit_set["molecularConcentration"]["oli_unit"],
    }
)
