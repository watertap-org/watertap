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

import pytest
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_equivalent
from watertap.tools.oli_api.util.fixed_keys_dict import output_unit_set, input_unit_set


@pytest.mark.unit
def test_fixed_keys_dict():
    with pytest.raises(RuntimeError):
        output_unit_set["invalid_key"] = "value"

    with pytest.raises(Exception):
        del output_unit_set["any_key"]

    with pytest.raises(RuntimeError):
        output_unit_set._check_value("mass", ["not_kg"])


input_unit_set_temp = {
    "inflows": {
        "oli_unit": "mg/L",
        "pyomo_unit": pyunits.mg / pyunits.L,
    },
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

output_unit_set_temp = {
    "enthalpy": input_unit_set_temp["enthalpy"],
    "mass": input_unit_set_temp["mass"],
    "pt": input_unit_set_temp["pressure"],
    "total": input_unit_set_temp["mass"],
    "liq1_phs_comp": input_unit_set_temp["mass"],
    "solid_phs_comp": input_unit_set_temp["mass"],
    "vapor_phs_comp": input_unit_set_temp["mass"],
    "liq2_phs_comp": input_unit_set_temp["mass"],
    "combined_phs_comp": input_unit_set_temp["mass"],
    "molecularConcentration": input_unit_set_temp["molecularConcentration"],
}

@pytest.mark.unit
def test_input_output_fixed_key_dicts():
    input_unit_set_temp.keys() == input_unit_set.keys()
    for k in input_unit_set_temp.keys():
        assert input_unit_set[k]["oli_unit"]==input_unit_set_temp[k]["oli_unit"]    
        assert_units_equivalent(input_unit_set[k]["pyomo_unit"], input_unit_set_temp[k]["pyomo_unit"])
    
    output_unit_set_temp.keys() == output_unit_set.keys()
    for k in output_unit_set_temp.keys():
        assert output_unit_set[k]["oli_unit"]==output_unit_set_temp[k]["oli_unit"]    
        assert_units_equivalent(output_unit_set[k]["pyomo_unit"], output_unit_set_temp[k]["pyomo_unit"])


@pytest.mark.unit
def test_input_unit_set():
    unit_set = input_unit_set
    # check defaults for one of the properties
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg/L"
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mg/L".lower()
    assert_units_equivalent(
        (unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg / pyunits.L
    )
    assert hasattr(
        (unit_set["molecularConcentration"]["pyomo_unit"]), "is_expression_type"
    )
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mg/L".lower()

    # reset oli_unit with a different unit, provided as string, and check that pyomo_unit follows along
    unit_set["molecularConcentration"]["oli_unit"] = "mol/L"
    assert unit_set["molecularConcentration"]["oli_unit"] == "mol/L"
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mol/L".lower()
    assert_units_equivalent(
        (unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mol / pyunits.L
    )

    # reset pyomo_unit with a different unit, provided as pint units, and check that oli_unit follows along
    unit_set["molecularConcentration"]["pyomo_unit"] = pyunits.mg
    assert_units_equivalent(
        (unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg
    )
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg"

    # reset pyomo_unit with a different unit, provided as string units, and check that oli_unit follows along
    unit_set["molecularConcentration"]["pyomo_unit"] = "mg/L"
    assert_units_equivalent(
        (unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg / pyunits.L
    )
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg/L"

    # reset oli_unit with a different unit, provided as pint units, and check that pyomo_unit follows along
    unit_set["molecularConcentration"]["oli_unit"] = pyunits.mol / pyunits.L
    assert unit_set["molecularConcentration"]["oli_unit"] == "mol/L".lower()
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mol/L".lower()
    assert_units_equivalent(
        (unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mol / pyunits.L
    )

    with pytest.raises(
        RuntimeError,
        match="Setting 1 as the value for oli_unit is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units.",
    ):
        unit_set["molecularConcentration"]["oli_unit"] = 1

    with pytest.raises(
        RuntimeError,
        match="Setting 1 as the value for pyomo_unit is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units.",
    ):
        unit_set["molecularConcentration"]["pyomo_unit"] = 1

    with pytest.raises(
        RuntimeError,
        match="Setting 1 as the value for molecularConcentration is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units.",
    ):
        unit_set["molecularConcentration"] = 1
