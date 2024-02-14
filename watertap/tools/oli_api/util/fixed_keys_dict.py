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
water_analysis_properties = FixedKeysDict(
    {
        "Temperature": {
            "group": "Properties",
            "name": "Temperature",
            "unit": input_unit_set["temperature"]["oli_unit"],
            "value": None,
        },
        "Pressure": {
            "group": "Properties",
            "name": "Pressure",
            "unit": input_unit_set["pressure"]["oli_unit"],
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
        "pH": {
            "group": "Properties",
            "name": "pH",
            "value": None,
        },
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
            "unit": input_unit_set["alkalinity"]["oli_unit"],
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
        "TIC": {
            "group": "Properties",
            "name": "TIC",
            "unit": input_unit_set["TIC"]["oli_unit"],
            "value": None,
        },
        "CO2GasFraction": {
            "group": "Properties",
            "name": "CO2GasFraction",
            "unit": input_unit_set["CO2GasFraction"]["oli_unit"],
            "value": None,
        },
        "AllowSolidsToForm": {
            "group": "Calculation Options",
            "name": "AllowSolidsToForm",
            "value": [True, False],
        },
        "CalcAlkalnity": {
            "group": "Calculation Options",
            "name": "CalcAlkalnity",
            "value": [False, True],
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
# This dictionary describes the stream outputs for OLI flash calculations
stream_output_options = FixedKeysDict(
    {
        "result": [
            "MBGComposition",
            "activityCoefficientsMBased",
            "activityCoefficientsXBased",
            "entropy",
            "entropyStandardStateXBased",
            "gibbsFreeEnergy",
            "gibbsFreeEnergyStandardStateXBased",
            "mobilities",
            "molecularConcentration",
            "selfDiffusivities",
            "totalMBGMoles",
            "totalMolecularMoles",
            "totalTrueMoles",
            "trueConcentration",
        ],
        "properties": [
            "absoluteViscosity",
            "density",
            "enthalpy",
            "entropy",
            "entropyStandardState",
            "gibbsFreeEnergy",
            "gibbsFreeEnergyStandardState",
            "hardness",
            "heatCapacity",
            "idealStandardLiquidVolume",
            "interfacialTension",
            "ionicStrength",
            "ionicStrengthMBased",
            "ionicStrengthXBased",
            "mass",
            "mixHeatCapacity",
            "molarElectricalConductivity",
            "orp",
            "osmoticPressure",
            "ph",
            "pressure",
            "relativeViscosity",
            "specificElectricalConductivity",
            "surfaceTension",
            "temperature",
            "thermalConductivity",
            "totalDissolvedSolids",
            "volume",
            "volumeStdConditions",
        ],
        "waterAnalysisOutput": [
            "addedIonsToBalance",
        ],
        "additionalProperties": [
            "kValuesMBased",
            "kValuesXBased",
            "prescalingIndex",
            "prescalingTendencies",
            "scalingIndex",
            "scalingTendencies",
            "vaporToInflowMoleFraction",
        ],
    }
)
