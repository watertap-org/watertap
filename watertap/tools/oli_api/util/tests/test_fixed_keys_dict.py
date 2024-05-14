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
from watertap.tools.oli_api.flash import Flash

@pytest.mark.unit
def test_fixed_keys_dict():
    with pytest.raises(RuntimeError):
        output_unit_set["invalid_key"] = "value"

    with pytest.raises(Exception):
        del output_unit_set["any_key"]

    with pytest.raises(RuntimeError):
        output_unit_set._check_value("mass", ["not_kg"])

    output_unit_set.pprint()


@pytest.mark.unit
def test_input_unit_set():
    unit_set = input_unit_set
    # check defaults for one of the properties
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg/L"
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mg/L".lower()
    assert_units_equivalent((unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg/pyunits.L)
    assert hasattr((unit_set["molecularConcentration"]["pyomo_unit"]), "is_expression_type")
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mg/L".lower()

    # reset oli_unit with a different unit, provided as string, and check that pyomo_unit follows along
    unit_set["molecularConcentration"]["oli_unit"] = "mol/L"
    assert unit_set["molecularConcentration"]["oli_unit"] == "mol/L"
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mol/L".lower()
    assert_units_equivalent((unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mol/pyunits.L)

    # reset pyomo_unit with a different unit, provided as pint units, and check that oli_unit follows along
    unit_set["molecularConcentration"]["pyomo_unit"] = pyunits.mg
    assert_units_equivalent((unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg)
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg"

    # reset pyomo_unit with a different unit, provided as string units, and check that oli_unit follows along
    unit_set["molecularConcentration"]["pyomo_unit"] = "mg/L"
    assert_units_equivalent((unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mg/pyunits.L)
    assert unit_set["molecularConcentration"]["oli_unit"] == "mg/L"

    # reset oli_unit with a different unit, provided as pint units, and check that pyomo_unit follows along
    unit_set["molecularConcentration"]["oli_unit"] = pyunits.mol/pyunits.L
    assert unit_set["molecularConcentration"]["oli_unit"] == "mol/L".lower()
    assert str(unit_set["molecularConcentration"]["pyomo_unit"]) == "mol/L".lower()
    assert_units_equivalent((unit_set["molecularConcentration"]["pyomo_unit"]), pyunits.mol/pyunits.L)

    with pytest.raises(RuntimeError, match="Setting 1 as the value for oli_unit is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units."):
        unit_set["molecularConcentration"]["oli_unit"] = 1

    with pytest.raises(RuntimeError, match="Setting 1 as the value for pyomo_unit is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units."):
        unit_set["molecularConcentration"]["pyomo_unit"] = 1

    with pytest.raises(RuntimeError, match="Setting 1 as the value for molecularConcentration is not permitted as a value for oli and pyomo units. Please enter units as a string type or pint units."):
        unit_set["molecularConcentration"] = 1