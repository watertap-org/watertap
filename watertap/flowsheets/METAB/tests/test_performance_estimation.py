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
import numpy as np
import os
import matplotlib

matplotlib.use("Agg")

from unittest.mock import MagicMock
from watertap.flowsheets.METAB.performance_estimation import (
    performance_estimation,
    display_performance,
    display_plot,
)

local_path = os.path.dirname(os.path.abspath(__file__))

dummy_components = ["S_su", "S_aa", "S_fa"]
dummy_inputs = ["x1", "x2"]


@pytest.fixture
def mock_surrogate():
    surrogate = MagicMock()
    surrogate.input_labels.return_value = dummy_inputs
    surrogate.output_labels.return_value = dummy_components

    def evaluate_surrogate(dataframe):
        # dataframe here contains only input columns (x1, x2),
        # so generate predictions from those rather than indexing output columns
        n = len(dataframe)
        preds = np.random.rand(n, len(dummy_components)) + 1.0
        return pd.DataFrame(preds, columns=dummy_components)

    surrogate.evaluate_surrogate.side_effect = evaluate_surrogate
    return surrogate


@pytest.fixture
def mock_dataframe():
    np.random.seed(42)
    n = 20
    data = {inp: np.random.rand(n) for inp in dummy_inputs}
    data.update({comp: np.random.rand(n) + 1.0 for comp in dummy_components})
    return pd.DataFrame(data)


def test_performance_estimation_returns_dataframe(mock_surrogate, mock_dataframe):
    result = performance_estimation(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert isinstance(result, pd.DataFrame)


def test_performance_estimation_one_row_per_output(mock_surrogate, mock_dataframe):
    result = performance_estimation(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert result.shape[0] == len(dummy_components)


def test_performance_estimation_columns(mock_surrogate, mock_dataframe):
    result = performance_estimation(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert set(result.columns) >= {"Comp", "R2", "RMSE", "MAE", "MSE", "maxAE"}


def test_performance_estimation_comp_values(mock_surrogate, mock_dataframe):
    result = performance_estimation(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert list(result["Comp"]) == dummy_components


def test_performance_estimation_missing_surrogate():
    with pytest.raises(ValueError, match="surrogate"):
        performance_estimation(surrogate=None, dataframe=pd.DataFrame())


def test_performance_estimation_missing_dataframe(mock_surrogate):
    with pytest.raises(ValueError, match="dataframe"):
        performance_estimation(surrogate=mock_surrogate, dataframe=None)


def test_display_performance_returns_dataframe(mock_surrogate, mock_dataframe):
    result = display_performance(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert isinstance(result, pd.DataFrame)


def test_display_performance_columns(mock_surrogate, mock_dataframe):
    result = display_performance(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert list(result.columns) == [
        "Predicted Variables",
        "R^2",
        "MAE",
        "MSE",
        "RMSE",
        "Max AE",
    ]


def test_display_performance_row_count(mock_surrogate, mock_dataframe):
    result = display_performance(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert result.shape[0] == len(dummy_components)


def test_display_performance_index_starts_at_one(mock_surrogate, mock_dataframe):
    result = display_performance(surrogate=mock_surrogate, dataframe=mock_dataframe)
    assert list(result.index) == list(range(1, len(dummy_components) + 1))


def test_display_plot_returns_figure_list(mock_surrogate, mock_dataframe, tmp_path):
    result = display_plot(
        surrogate=mock_surrogate,
        dataframe=mock_dataframe,
        path=str(tmp_path),
        show=False,
    )
    assert isinstance(result, list)
    assert len(result) == len(dummy_components)


def test_display_plot_saves_pdf(mock_surrogate, mock_dataframe, tmp_path):
    display_plot(
        surrogate=mock_surrogate,
        dataframe=mock_dataframe,
        method="poly",
        path=str(tmp_path),
        show=False,
    )
    expected_pdf = tmp_path / "poly_parity.pdf"
    assert expected_pdf.exists()


def test_display_plot_default_path(
    mock_surrogate, mock_dataframe, tmp_path, monkeypatch
):
    # Redirect local_path so the default path stays inside tmp_path
    import watertap.flowsheets.METAB.performance_estimation as pe_module

    monkeypatch.setattr(pe_module, "local_path", str(tmp_path))
    os.makedirs(tmp_path / "results", exist_ok=True)

    result = display_plot(
        surrogate=mock_surrogate,
        dataframe=mock_dataframe,
        show=False,
    )
    assert isinstance(result, list)
