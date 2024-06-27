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
import pandas as pd
import pytest

from pyomo.environ import units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from watertap.core.util.chemistry import (
    get_charge,
    get_molar_mass,
    get_molar_mass_quantity,
    get_periodic_table,
)


@pytest.mark.unit
def test_charge_exceptions():
    with pytest.raises(IOError) as excinfo:
        val = "Na_2#"
        get_charge(val)
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


@pytest.fixture
def periodic_table() -> pd.DataFrame:
    return get_periodic_table()


@pytest.mark.unit
def test_periodic_table_headers(periodic_table):
    test_headers = ["AtomicMass", "Symbol"]
    assert all(header in periodic_table.columns for header in test_headers)
    size = {"cols": 28, "rows": 118}
    assert len(periodic_table.columns) == size["cols"]
    assert len(periodic_table.values) == size["rows"]
