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

import pytest

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    watertap_to_oli,
    get_oli_names,
    oli_reverse_lookup,
    names_db,
)


@pytest.mark.parametrize(
    "input_name,expected_output",
    [
        ("NaOH", "NAOH"),
        ("B[OH]3", "BOH3"),
        ("B[OH]4_-", "BOH4ION"),
        ("K_+", "KION"),
        ("Cl_-", "CLION"),
        ("Mg_2+", "MGION"),
        ("HCO3_2-", "HCO3ION"),
    ],
)
def test_watertap_to_oli(input_name: str, expected_output: str):
    assert watertap_to_oli(input_name).oli_name == expected_output


@pytest.mark.unit
def test_case_exception():
    with pytest.raises(IOError) as excinfo:
        val = "abc_2+"
        watertap_to_oli(val)
    assert (
        str(excinfo.value)
        == f" At least 1 uppercase letter is required to specify a molecule, not '{val}'."
    )


def test_charge_exception():
    with pytest.raises(IOError) as excinfo:
        val = "Na_2#"
        watertap_to_oli(val)
    assert str(excinfo.value) == " Only + and - are valid charge indicators."


@pytest.mark.unit
def test_reverse_lookup():
    assert oli_reverse_lookup("NAION", names_db) == "Na_+"


@pytest.mark.unit
def test_reverse_lookup_exception():
    with pytest.raises(IOError) as excinfo:
        oli_reverse_lookup("Na_-", names_db)
    assert (
        str(excinfo.value)
        == " Component Na_- not found in names_db. Update this dictionary to hard code additional OLI names."
    )


@pytest.mark.unit
def test_get_oli_names():
    source_dict = {"Ca_2+": 500, "Cl_-": 1000, "Mg_2+": 500, "SO4_2-": 50}
    source_dict = get_oli_names(source_dict)


@pytest.mark.unit
def verify_names_db_contents():
    initial_names_db = {
        "NAION": "Na_+",
        "CLION": "Cl_-",
        "SO4ION": "SO4_2-",
        "MGION": "Mg_2+",
        "CAION": "Ca_2+",
        "KION": "K_+",
        "HCO3ION": "HCO3_-",
        "NA2CO3": "Na2CO3",
        "CO2": "CO2",
        "H2O": "H2O",
    }

    for k, v in initial_names_db.items():
        assert v == names_db[k]
