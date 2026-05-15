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
import pandas as pd
import os
from watertap.flowsheets.METAB.model_evaluation import (
    get_input_data,
    get_eff_fr,
    get_ch4_fr,
    get_h2_fr,
    get_r1_ex_biogas_fr,
    get_mass_flowrate,
    collect_results,
    run_model,
    export_output_data,
)

# third-party
try:
    import exposan
except ImportError:
    exposan = None

local_path = os.path.dirname(os.path.abspath(__file__))


def test_get_input_data_no_filename():
    df = get_input_data()
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["inf_fr", "temp", "hrt"]
    assert list(df["inf_fr"]) == [5, 5, 5]
    assert list(df["temp"]) == [20, 25, 30]
    assert list(df["hrt"]) == [12, 13, 14]


def test_get_input_data_w_filename(tmp_path):
    input_data_path = tmp_path / "input_data.csv"
    dummy = pd.DataFrame(
        {
            "inf_fr": [5, 5, 5, 5, 5],
            "temp": [20, 25, 30, 35, 40],
            "hrt": [1, 2, 3, 4, 5],
        }
    )
    dummy.to_csv(input_data_path, index=False)
    df = get_input_data(str(input_data_path))
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["inf_fr", "temp", "hrt"]
    assert df.shape == (5, 3)


@pytest.fixture
def mock_case():
    class MockStream:
        def __init__(self, components):
            self.components = components

        def _info(self, layout, T, P, flow, composition, N, IDs):
            lines = "Stream\nflow (kmol/hr):\n"
            for name, value in self.components.items():
                lines += f" {name}   {value}\n"
            return lines

    class MockCase:
        def __init__(self):
            self.outs = [
                None,
                MockStream(
                    {"H2O": 1.0, "CH4": 20.0, "CO2": 10.0}
                ),  # index 1 for get_h2_fr,
                MockStream(
                    {"H2O": 10.0, "CH4": 75.0, "CO2": 15.0}
                ),  # index 2 for get_ch4_fr
                MockStream(
                    {"H2O": 100.0, "CH4": 50.0, "CO2": 25.0}
                ),  # index 3 for get_eff_fr
                MockStream(
                    {"H2O": 200.0, "CH4": 80.0, "CO2": 40.0}
                ),  # index 4 for get_r1_ex_biogas_fr
            ]

    return MockCase()


def test_get_eff_fr_no_case(capsys):
    result = get_eff_fr()
    captured = capsys.readouterr()
    assert "The system is off" in captured.out
    assert result is None


def test_get_eff_fr_w_case(mock_case):
    df = get_eff_fr(case=mock_case)
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (1, 3)
    assert list(df.columns) == ["H2O", "CH4", "CO2"]
    assert df.iloc[0]["H2O"] == 100.0
    assert df.iloc[0]["CO2"] == 25.0
    assert df.iloc[0]["CH4"] == 50.0


def test_get_eff_fr_w_existing_df(mock_case):
    df = get_eff_fr(case=mock_case)
    df = get_eff_fr(case=mock_case, df=df)
    assert df.shape == (2, 3)
    assert df.iloc[0]["H2O"] == 100.0
    assert df.iloc[0]["CO2"] == 25.0
    assert df.iloc[0]["CH4"] == 50.0
    assert df.iloc[1]["H2O"] == [100.0]
    assert df.iloc[1]["CO2"] == [25.0]
    assert df.iloc[1]["CH4"] == [50.0]


def test_get_ch4_fr_no_case(capsys):
    result = get_ch4_fr()
    captured = capsys.readouterr()
    assert "The system is off" in captured.out
    assert result is None


def test_get_ch4_fr_w_case(mock_case):
    df = get_ch4_fr(case=mock_case)
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (1, 3)
    assert list(df.columns) == ["H2O", "CH4", "CO2"]
    assert df.iloc[0]["H2O"] == 10.0
    assert df.iloc[0]["CH4"] == 75.0
    assert df.iloc[0]["CO2"] == 15.0


def test_get_ch4_fr_w_existing_df(mock_case):
    df = get_ch4_fr(case=mock_case)
    df = get_ch4_fr(case=mock_case, df=df)
    assert df.shape == (2, 3)
    assert df.iloc[0]["H2O"] == 10.0
    assert df.iloc[0]["CH4"] == 75.0
    assert df.iloc[0]["CO2"] == 15.0
    assert df.iloc[1]["H2O"] == [10.0]
    assert df.iloc[1]["CH4"] == [75.0]
    assert df.iloc[1]["CO2"] == [15.0]


def test_get_h2_fr_no_case(capsys):
    result = get_h2_fr()
    captured = capsys.readouterr()
    assert "The system is off" in captured.out
    assert result is None


def test_get_h2_fr_w_case(mock_case):
    df = get_h2_fr(case=mock_case)
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (1, 3)
    assert list(df.columns) == ["H2O", "CH4", "CO2"]
    assert df.iloc[0]["H2O"] == 1.0
    assert df.iloc[0]["CH4"] == 20.0
    assert df.iloc[0]["CO2"] == 10.0


def test_get_h2_fr_w_existing_df(mock_case):
    df = get_h2_fr(case=mock_case)
    df = get_h2_fr(case=mock_case, df=df)
    assert df.shape == (2, 3)
    assert df.iloc[0]["H2O"] == 1.0
    assert df.iloc[0]["CH4"] == 20.0
    assert df.iloc[0]["CO2"] == 10.0
    assert df.iloc[1]["H2O"] == [1.0]
    assert df.iloc[1]["CH4"] == [20.0]
    assert df.iloc[1]["CO2"] == [10.0]


def test_get_r1_ex_biogas_fr_no_case(capsys):
    result = get_r1_ex_biogas_fr()
    captured = capsys.readouterr()
    assert "The system is off" in captured.out
    assert result is None


def test_get_r1_ex_biogas_fr_w_case(mock_case):
    df = get_r1_ex_biogas_fr(case=mock_case)
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (1, 3)
    assert list(df.columns) == ["H2O", "CH4", "CO2"]
    assert df.iloc[0]["H2O"] == 200.0
    assert df.iloc[0]["CH4"] == 80.0
    assert df.iloc[0]["CO2"] == 40.0


def test_get_r1_ex_biogas_fr_w_existing_df(mock_case):
    df = get_r1_ex_biogas_fr(case=mock_case)
    df = get_r1_ex_biogas_fr(case=mock_case, df=df)
    assert df.shape == (2, 3)
    assert df.iloc[0]["H2O"] == 200.0
    assert df.iloc[0]["CH4"] == 80.0
    assert df.iloc[0]["CO2"] == 40.0
    assert df.iloc[1]["H2O"] == [200.0]
    assert df.iloc[1]["CH4"] == [80.0]
    assert df.iloc[1]["CO2"] == [40.0]


@pytest.fixture
def mock_stream():
    class MockComponents:
        def __str__(self):
            return "Components(H2O,CO2,CH4)"

    class MockStream:
        def __init__(self):
            self.components = MockComponents()
            self.state = [100.0, 50.0, 25.0, 10.0]

    return MockStream()


def test_get_mass_flowrate(mock_stream):
    df = get_mass_flowrate(stream=mock_stream)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["H2O", "CO2", "CH4", "Volumetric Flowrate"]
    assert df.shape == (1, 4)
    assert df.iloc[0]["H2O"] == 100.0
    assert df.iloc[0]["CO2"] == 50.0
    assert df.iloc[0]["CH4"] == 25.0
    assert df.iloc[0]["Volumetric Flowrate"] == 10.0


def test_get_mass_flowrate_drops_zero_columns(mock_stream):
    mock_stream.state = [100.0, 0.0, 25.0, 10.0]
    df = get_mass_flowrate(stream=mock_stream)
    assert "CO2" not in df.columns


def test_get_mass_flowrate_w_existing_df(mock_stream):
    df = get_mass_flowrate(stream=mock_stream)
    df = get_mass_flowrate(stream=mock_stream, df=df)
    assert df.shape == (2, 4)
    assert df.iloc[1]["H2O"] == [100.0]
    assert df.iloc[1]["Volumetric Flowrate"] == [10.0]


@pytest.fixture
def mock_collection():
    class MockComponents:
        def __init__(self, names):
            self.names = names

        def __str__(self):
            return f"Components({','.join(self.names)})"

    class MockStream:
        def __init__(self, components, state):
            self.components = MockComponents(components)
            self.state = state

        def _info(self, layout, T, P, flow, composition, N, IDs):
            lines = "Stream\nflow (kmol/hr):\n"
            for name, value in zip(self.components.names, self.state[:-1]):
                lines += f" {name}   {value}\n"
            return lines

    class MockCase:
        def __init__(self):
            self.outs = [
                None,
                MockStream(
                    ["H2O", "H2", "CO2"], [10.0, 5.0, 2.0, 1.0]
                ),  # index 1 for get_h2_fr
                MockStream(
                    ["H2O", "CH4", "CO2"], [10.0, 75.0, 15.0, 5.0]
                ),  # index 2 for get_ch4_fr
                MockStream(
                    ["H2O", "CO2", "CH4"], [100.0, 50.0, 25.0, 10.0]
                ),  # index 3 for get_eff_fr
                MockStream(
                    ["H2O", "CH4", "CO2"], [5.0, 40.0, 10.0, 2.0]
                ),  # index 4 for get_r1_ex_biogas_fr
            ]

    return MockCase()


def test_collect_results(mock_collection):
    result = collect_results(case=mock_collection)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 1
    ch4_cols = [col for col in result.columns if col.startswith("bge2_")]
    assert list(ch4_cols) == [
        "bge2_H2O",
        "bge2_CH4",
        "bge2_CO2",
        "bge2_Volumetric Flowrate",
    ]
    h2_cols = [col for col in result.columns if col.startswith("bgr2_")]
    assert list(h2_cols) == [
        "bgr2_H2O",
        "bgr2_H2",
        "bgr2_CO2",
        "bgr2_Volumetric Flowrate",
    ]
    bgs_cols = [col for col in result.columns if col.startswith("bge1_")]
    assert list(bgs_cols) == [
        "bge1_H2O",
        "bge1_CH4",
        "bge1_CO2",
        "bge1_Volumetric Flowrate",
    ]


def test_collect_results_w_existing_df(mock_collection):
    results = collect_results(case=mock_collection)
    results = collect_results(case=mock_collection, results=results)
    assert results.shape[0] == 2


def test_collect_results_mass_false(mock_collection):
    result = collect_results(case=mock_collection, mass=False)
    assert isinstance(result, pd.DataFrame)


def test_run_model_integration(tmp_path):
    if exposan is None:
        pytest.skip("exposan not available")

    input_data = pd.DataFrame(
        {
            "inf_fr": [5, 5, 5],
            "temp": [20, 25, 30],
            "hrt": [1, 5, 10],
        }
    )
    output_data = run_model(input_data)

    assert isinstance(output_data, pd.DataFrame)
    assert output_data.shape[0] == len(input_data)


def test_export_output_data_integration(tmp_path):
    if exposan is None:
        pytest.skip("exposan not available")

    input_data = pd.DataFrame(
        {
            "inf_fr": [5, 5, 5],
            "temp": [20, 25, 30],
            "hrt": [1, 5, 10],
        }
    )
    output_data = run_model(input_data)

    output_csv = str(tmp_path / "output_data.csv")
    export_output_data(output_data, filename=output_csv)

    assert os.path.exists(output_csv)

    df_read_back = pd.read_csv(output_csv)
    assert df_read_back.shape == output_data.shape
