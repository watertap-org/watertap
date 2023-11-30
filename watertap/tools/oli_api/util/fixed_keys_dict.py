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


# TODO: allow user to fully specify units (i.e., set to None and force unit specification)
default_water_analysis_properties = FixedKeysDict(
    {
        "Temperature": {
            "group": "Properties",
            "name": "Temperature",
            "unit": "K",
            "value": None,
        },
        "Pressure": {
            "group": "Properties",
            "name": "Pressure",
            "unit": "Pa",
            "value": None,
        },
        "ElectroNeutralityBalanceType": {
            "group": "Electroneutrality Options",
            "name": "ElectroNeutralityBalanceType",
            "value": [
                "DominantIon",
                "ProrateCations",
                "ProrateAnions",
                "Prorate",
                "AutoNACL",
                "MakeupIon",
            ],
        },
        "MakeupIonBaseTag": {
            "group": "Electroneutrality Options",
            "name": "MakeupIonBaseTag",
            "value": None,
        },
        "CalcType": {
            "group": "Calculation Options",
            "name": "CalcType",
            "value": [
                "EquilCalcOnly",
                "ReconcilePh",
                "ReconcilePhAndAlkalinity",
                "ReconcilePhAndAlkalinityAndTic",
                "ReconcileCo2Gas",
            ],
        },
        "pH": {"group": "Properties", "name": "pH", "value": None},
        "PhAcidTitrant": {
            "group": "Calculation Options",
            "name": "PhAcidTitrant",
            "value": None,
        },
        "PhBaseTitrant": {
            "group": "Calculation Options",
            "name": "PhBaseTitrant",
            "value": None,
        },
        "Alkalinity": {
            "group": "Properties",
            "name": "Alkalinity",
            "unit": "mg HCO3/L",
            "value": None,
        },
        "AlkalinityPhTitrant": {
            "group": "Calculation Options",
            "name": "AlkalinityPhTitrant",
            "value": None,
        },
        "AlkalinityTitrationEndPointPh": {
            "group": "Properties",
            "name": "AlkalinityTitrationEndPointpH",
            "value": None,
        },
        "TIC": {"group": "Properties", "name": "TIC", "unit": "mol C/L", "value": None},
        "CO2GasFraction": {
            "group": "Properties",
            "name": "CO2GasFraction",
            "unit": "mole %",
            "value": None,
        },
        "AllowSolidsToForm": {
            "group": "Calculation Options",
            "name": "AllowSolidsToForm",
            "value": [False, True],
        },
        "CalcAlkalnity": {
            "group": "Calculation Options",
            "name": "CalcAlkalnity",
            "value": [False, True],
        },
    }
)

# Full list of available optional inputs available: https://devdocs.olisystems.com/optional-inputs
default_optional_properties = FixedKeysDict(
    {
        "electricalConductivity": False,
        "viscosity": False,
        "selfDiffusivityAndMobility": False,
        "heatCapacity": False,
        "thermalConductivity": False,
        "surfaceTension": False,
        "interfacialTension": False,
        "volumeStdConditions": False,
        "prescalingTendenciesEstimated": False,
        "prescalingIndexEstimated": False,
        "prescalingTendencies": False,
        "prescalingIndex": False,
        "prescalingTendenciesRigorous": False,
        "prescalingIndexRigorous": False,
        "scalingTendencies": False,
        "scalingIndex": False,
        "hardness": False,
        "ionicStrengthXBased": False,
        "ionicStrengthMBased": False,
        "totalDissolvedSolids": False,
        "vaporToInflowMoleFraction": False,
        "partialPressure": False,
        "vaporDiffusivityMatrix": False,
        "entropyStream": False,
        "entropySpecies": False,
        "entropyStreamStandardState": False,
        "entropySpeciesStandardState": False,
        "gibbsEnergyStream": False,
        "gibbsEnergySpecies": False,
        "gibbsEnergyStreamStandardState": False,
        "gibbsEnergySpeciesStandardState": False,
        "activityCoefficientsXBased": False,
        "activityCoefficientsMBased": False,
        "fugacityCoefficients": False,
        "vaporFugacity": False,
        "kValuesXBased": False,
        "kValuesMBased": False,
        "MBGComposition": False,
        "materialBalanceGroup": False,
    }
)

default_input_unit_set = FixedKeysDict(
    {
        "molecularConcentration": {"oli_unit": "mg/L", "pyomo_unit": pyunits.mg/pyunits.L},
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
