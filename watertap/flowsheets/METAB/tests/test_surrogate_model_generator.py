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
import os

import pandas as pd
import pytest

from watertap.flowsheets.METAB.surrogate_model_generator import (
    get_data,
    outputs_selections,
    # gen_surrogate_model,
)

local_path = os.path.dirname(os.path.abspath(__file__))

input_columns = ["inf_fr", "temp", "hrt"]

output_columns = [
    "S_su",
    "S_aa",
    "S_fa",
    "S_va",
    "S_bu",
    "S_pro",
    "S_ac",
    "S_h2",
    "S_ch4",
    "S_IC",
    "S_IN",
    "S_I",
    "X_c",
    "X_ch",
    "X_pr",
    "X_li",
    "X_su",
    "X_aa",
    "X_fa",
    "X_c4",
    "X_pro",
    "X_ac",
    "X_h2",
    "X_I",
    "VolumetricFlowrate",
]


@pytest.fixture
def csv_files(tmp_path):
    input_df = pd.DataFrame(
        {
            "inf_fr": [5, 5, 5],
            "temp": [20, 25, 30],
            "hrt": [12, 13, 14],
        }
    )
    output_df = pd.DataFrame({col: [0.1, 0.2, 0.3] for col in output_columns})

    input_file = tmp_path / "input_data.csv"
    output_file = tmp_path / "output_data.csv"
    input_df.to_csv(input_file, index=False)
    output_df.to_csv(output_file, index=False)

    return str(input_file), str(output_file)


@pytest.fixture
def data(csv_files):
    input_data, output_data = csv_files
    feed, input, output = get_data(
        input_data_file=input_data,
        output_data_file=output_data,
    )
    output = outputs_selections(output)
    feed = pd.concat([input, output], axis=1)
    return feed, input, output


@pytest.fixture
def small_data_sample(data):
    feed, input, output = data
    return (
        feed.iloc[:10].reset_index(drop=True),
        input.iloc[:10].reset_index(drop=True),
        output.iloc[:10].reset_index(drop=True),
    )


def test_get_data(csv_files):
    input_data, output_data = csv_files
    feed, input, output = get_data(
        input_data_file=input_data,
        output_data_file=output_data,
    )
    assert isinstance(feed, pd.DataFrame)
    assert isinstance(input, pd.DataFrame)
    assert isinstance(output, pd.DataFrame)
    assert list(input.columns) == input_columns
    assert len(feed.columns) == len(input.columns) + len(output.columns)
    assert len(feed) == len(input) == len(output)
    assert "Unnamed: 0" not in output.columns


def test_get_data_missing_input_file_raises(csv_files):
    _, output_data = csv_files
    with pytest.raises(FileNotFoundError):
        get_data(
            input_data_file="/nonexistent/path/input.csv",
            output_data_file=output_data,
        )


def test_get_data_missing_output_file_raises(csv_files):
    input_data, _ = csv_files
    with pytest.raises(FileNotFoundError):
        get_data(
            input_data_file=input_data,
            output_data_file="/nonexistent/path/output.csv",
        )


# TODO (diagnosis from AI): PySMO polynomial regression (polynomial_regression_fitting) creates a solution.pickle
#  file but does not reliably close the file handle, leading to a PytestUnraisableExceptionWarning / ResourceWarning:
#  unclosed file <_io.FileIO [closed]> during Python garbage collection and pytest teardown.
# def test_gen_surrogate_model_poly_saves_json(small_data_sample, tmp_path):
#     feed, input, output = small_data_sample
#     gen_surrogate_model(
#         method="poly", feed_data=feed, input_data=input, output_data=output, path=tmp_path
#     )
#     json_path = tmp_path / "poly_surrogate.json"
#     assert os.path.exists(json_path)
#     with open(json_path) as f:
#         data = json.load(f)
#     assert isinstance(data, dict)
#     assert os.path.exists(tmp_path / "poly_parity.pdf")
#
#
# def test_gen_surrogate_model_poly_none_feed_data(small_data_sample, tmp_path):
#     _, input, output = small_data_sample
#     gen_surrogate_model(method="poly", feed_data=None, input_data=input, output_data=output, path=tmp_path)
#     assert os.path.exists(tmp_path / "poly_surrogate.json")
#
#
# def test_gen_surrogate_model_kri_saves_json(small_data_sample, tmp_path):
#     feed, input, output = small_data_sample
#     gen_surrogate_model(method="kri", feed_data=feed, input_data=input, output_data=output, path=tmp_path)
#     json_path = tmp_path / "kri_surrogate.json"
#     assert os.path.exists(json_path)
#     with open(json_path) as f:
#         data = json.load(f)
#     assert isinstance(data, dict)
#     assert os.path.exists(tmp_path / "kri_parity.pdf")
#
#
# def test_gen_surrogate_model_kri_none_feed_data(small_data_sample, tmp_path):
#     _, input, output = small_data_sample
#     gen_surrogate_model(method="kri", feed_data=None, input_data=input, output_data=output, path=tmp_path)
#     assert os.path.exists(tmp_path / "kri_surrogate.json")
#
#
# def test_gen_surrogate_model_rbf_saves_json(small_data_sample, tmp_path):
#     feed, input, output = small_data_sample
#     gen_surrogate_model(method="rbf", feed_data=feed, input_data=input, output_data=output, path=tmp_path)
#     json_path = tmp_path / "rbf_surrogate.json"
#     assert os.path.exists(json_path)
#     with open(json_path) as f:
#         data = json.load(f)
#     assert isinstance(data, dict)
#     assert os.path.exists(tmp_path / "rbf_parity.pdf")
#
#
# def test_gen_surrogate_model_rbf_none_feed_data(small_data_sample, tmp_path):
#     _, input, output = small_data_sample
#     gen_surrogate_model(method="rbf", feed_data=None, input_data=input, output_data=output, path=tmp_path)
#     assert os.path.exists(tmp_path / "rbf_surrogate.json")
#
#
# def test_gen_surrogate_model_unknown_method_raises(small_data_sample, tmp_path):
#     feed, input, output = small_data_sample
#     with pytest.raises(ValueError):
#         gen_surrogate_model(
#             method="unsupported_method",
#             feed_data=feed,
#             input_data=input,
#             output_data=output,
#             path=tmp_path,
#         )
