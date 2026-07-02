#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest
import csv
import os
from watertap.flowsheets.METAB.input_space_generation import create_samples

input_var_info = {
    "inf_fr": (5, 10),
    "temp": (22, 35),
    "hrt": (1, 12),
}

local_path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def temp_csv(tmp_path):
    return str(tmp_path / "test_input_data.csv")


def test_no_method_message(capsys, temp_csv):
    create_samples(method=None, input_var_info=input_var_info, csv_file=temp_csv)
    captured = capsys.readouterr()
    assert "Please pick a sampling method" in captured.out


def test_lhs_method(temp_csv):
    create_samples(
        method="LHS",
        input_var_info=input_var_info,
        sample_numbers=20,
        csv_file=temp_csv,
    )
    assert os.path.exists(temp_csv)

    with open(temp_csv, "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        rows = list(reader)

    assert headers == list(input_var_info.keys())
    assert len(rows) == 20


def test_lhs_samples_within_bounds(temp_csv):
    create_samples(
        method="LHS",
        input_var_info=input_var_info,
        sample_numbers=20,
        csv_file=temp_csv,
    )
    with open(temp_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            assert (
                input_var_info["inf_fr"][0]
                <= float(row["inf_fr"])
                <= input_var_info["inf_fr"][1]
            )
            assert (
                input_var_info["temp"][0]
                <= float(row["temp"])
                <= input_var_info["temp"][1]
            )
            assert (
                input_var_info["hrt"][0]
                <= float(row["hrt"])
                <= input_var_info["hrt"][1]
            )
