###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Tests for WaterTAP database wrapper
"""
import pytest
import os

from watertap.core.wt_database import Database


@pytest.mark.unit
def test_default_path():
    db = Database()

    assert os.path.normpath(db._dbpath) == os.path.normpath(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "data",
            "techno_economic",
        )
    )


@pytest.mark.unit
def test_invalid_path():
    with pytest.raises(
        OSError,
        match="Could not find requested path foo. Please "
        "check that this path exists.",
    ):
        Database(dbpath="foo")


@pytest.mark.unit
def test_custom_path():
    # Pick a path we know will exist, even if it isn't a data folder
    db = Database(
        dbpath=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"
        )
    )

    assert db._dbpath == os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"
    )


class TestDatabase:
    @pytest.fixture(scope="class")
    def db(self):
        return Database()

    @pytest.mark.unit
    def test_component_list(self, db):
        # First, check that _component_list exists and is None
        assert db._component_list is None

        # Then check that component list retrieved is a dict
        assert isinstance(db.component_list, dict)

    @pytest.mark.unit
    def test_get_technology(self, db):
        assert db._cached_files == {}

        data = db._get_technology("nanofiltration")

        assert "nanofiltration" in db._cached_files
        assert db._cached_files["nanofiltration"] is data

        # Check for a few expected keys to check what we got back looks right
        assert "default" in data
        assert "removal_frac_mass_comp" in data["default"]
        assert "energy_electric_flow_vol_inlet" in data["default"]
        assert "recovery_frac_mass_H2O" in data["default"]

    @pytest.mark.unit
    def test_get_technology_invalid(self, db):
        with pytest.raises(KeyError, match="Could not find entry for foo in database."):
            db._get_technology("foo")

        assert len(db._cached_files) == 1
        assert "foo" not in db._cached_files

    @pytest.mark.unit
    def test_get_unit_operation_parameters_default(self, db):
        data = db.get_unit_operation_parameters("nanofiltration")

        assert data == db._cached_files["nanofiltration"]["default"]

        # Check for a few expected keys to check what we got back looks right
        assert "removal_frac_mass_comp" in data
        assert "energy_electric_flow_vol_inlet" in data
        assert "recovery_frac_mass_H2O" in data

    @pytest.mark.unit
    def test_get_unit_operation_parameters_invalid_subtype(self, db):
        with pytest.raises(
            KeyError,
            match="Received unrecognised subtype foo for " "technology nanofiltration.",
        ):
            db.get_unit_operation_parameters("nanofiltration", subtype="foo")

    @pytest.mark.unit
    def test_get_unit_operation_parameters_single_subtype(self, db):
        # First, insert some data for a subtype into nanofiltration entry
        db._cached_files["nanofiltration"]["subtype1"] = {
            "recovery_frac_mass_H2O": "overloaded",
            "new_param": True,
        }

        # Load data for subtype
        data = db.get_unit_operation_parameters("nanofiltration", subtype="subtype1")

        # Check data
        for k, v in data.items():
            if k == "recovery_frac_mass_H2O":
                assert v == "overloaded"
            elif k == "new_param":
                assert v is True
            else:
                # All other entries should match defaults
                assert v == db._cached_files["nanofiltration"]["default"][k]

    @pytest.mark.unit
    def test_get_unit_operation_parameters_invalid_list_of_subtype(self, db):
        with pytest.raises(
            KeyError,
            match="Received unrecognised subtype foo for " "technology nanofiltration.",
        ):
            db.get_unit_operation_parameters(
                "nanofiltration", subtype=["subtype1", "foo"]
            )

    @pytest.mark.unit
    def test_get_unit_operation_parameters_multi_subtype(self, db):
        # First, insert some data for a 2nd subtype into nanofiltration entry
        db._cached_files["nanofiltration"]["subtype2"] = {
            "recovery_frac_mass_H2O": "overloaded_again",
            "new_param_2": False,
        }

        # Load data for subtype
        data = db.get_unit_operation_parameters("nanofiltration", subtype="subtype2")

        # Check data
        for k, v in data.items():
            if k == "recovery_frac_mass_H2O":
                assert v == "overloaded_again"
            elif k == "new_param":
                assert v is True
            elif k == "new_param_2":
                assert v is False
            else:
                # All other entries should match defaults
                assert v == db._cached_files["nanofiltration"]["default"][k]

    @pytest.mark.unit
    def test_get_unit_operation_parameters_subtype_argument(self, db):
        with pytest.raises(
            TypeError,
            match="Unexpected type for subtype 12: must be " "string or list like.",
        ):
            db.get_unit_operation_parameters("nanofiltration", subtype=12)

    @pytest.mark.unit
    def test_get_solute_set_default(self, db):
        comp_set = db.get_solute_set()

        assert comp_set == [
            "boron",
            "bromide",
            "calcium",
            "chloride",
            "magnesium",
            "potassium",
            "sodium",
            "strontium",
            "sulfate",
            "tds",
            "tss",
        ]

    @pytest.mark.unit
    def test_get_solute_set_specified(self, db):
        comp_set = db.get_solute_set("seawater")

        assert comp_set == [
            "boron",
            "bromide",
            "calcium",
            "chloride",
            "magnesium",
            "potassium",
            "sodium",
            "strontium",
            "sulfate",
            "tds",
            "tss",
        ]

    @pytest.mark.unit
    def test_get_solute_set_no_default(self, db):
        # First, delete default entry from database
        del db._cached_files["water_sources"]["default"]

        with pytest.raises(
            KeyError,
            match="Database has not defined a default water "
            "source and none was provided.",
        ):
            db.get_solute_set()

    @pytest.mark.unit
    def test_flush_cache(self, db):
        db.flush_cache()

        assert db._cached_files == {}
