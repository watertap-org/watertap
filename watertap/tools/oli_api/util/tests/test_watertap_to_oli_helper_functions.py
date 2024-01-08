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
    get_charge,
    get_molar_mass,
    get_molar_mass_quantity,
)
from pyomo.environ import units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

__author__ = "Paul Vecchiarelli, Adam Atia"


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


@pytest.mark.unit
def test_charge_exceptions():
    with pytest.raises(IOError) as excinfo:
        val = "Na_2#"
        watertap_to_oli(val)
    assert (
        str(excinfo.value)
        == "Only + and - are valid charge indicators and neither was provided in 'Na_2#'."
    )

    with pytest.raises(
        IOError,
        match="Charge sign could not be determined from the string 'target_ion'",
    ):
        get_charge("target_ion")
    with pytest.raises(
        IOError, match="Charge could not be determined from the string 'my_target_ion'"
    ):
        get_charge("my_target_ion")


@pytest.mark.unit
def test_get_charge():
    z = {
        "NaCl": 0,
        "Na_+": 1,
        "Cl_-": -1,
        "Ca_2+": 2,
        "SO4_2-": -2,
    }
    for solute_name, charge_value in z.items():
        assert get_charge(solute_name) == charge_value


@pytest.mark.unit
def test_get_mw():
    z = {
        "NaCl": 58.44,
        "Na_+": 22.99,
        "Cl_-": 35.45,
        "Ca_2+": 40.08,
        "SO4_2-": 96.066,
    }
    for solute_name, mw_value in z.items():
        assert get_molar_mass(solute_name) == mw_value


@pytest.mark.unit
def test_get_mw_exception():
    with pytest.raises(
        IOError, match="Molecular weight data could not be found for foo."
    ):
        get_molar_mass("foo")


@pytest.mark.unit
def test_get_mw_quantity():
    z = {
        "NaCl": 58.44,
        "Na_+": 22.99,
        "Cl_-": 35.45,
        "Ca_2+": 40.08,
        "SO4_2-": 96.066,
    }
    for solute_name, mw_value in z.items():
        result = get_molar_mass_quantity(solute_name)
        assert value(result) == mw_value / 1000
        assert_units_equivalent(result, pyunits.kg / pyunits.mol)

        result2 = get_molar_mass_quantity(solute_name, units=pyunits.g / pyunits.mol)
        assert value(result2) == mw_value
        assert_units_equivalent(result2, pyunits.g / pyunits.mol)


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
