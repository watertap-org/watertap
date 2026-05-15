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
import json

from watertap.flowsheets.METAB.performance_estimation import (
    performance_estimation,
    display_performance,
    display_plot,
)

# third-party
try:
    import IPython
except ImportError:
    IPython = None

local_path = os.path.dirname(os.path.abspath(__file__))

dummy_components = ["S_su", "S_aa", "S_fa"]


@pytest.fixture
def surrogate_path(tmp_path):
    poly_data = {
        "model_encoding": {
            comp: {
                "attr": {
                    "errors": {
                        "MAE": 0.01,
                        "MSE": 0.0001,
                        "R2": 0.98,
                        "Adjusted R2": 0.97,
                    }
                }
            }
            for comp in dummy_components
        }
    }
    kri_data = {
        "model_encoding": {
            comp: {"attr": {"training_R2": 0.96, "training_rmse": 0.04}}
            for comp in dummy_components
        }
    }
    rbf_data = {
        "model_encoding": {
            comp: {"attr": {"R2": 0.95, "rmse": 0.05}} for comp in dummy_components
        }
    }

    for method, data in [("poly", poly_data), ("kri", kri_data), ("rbf", rbf_data)]:
        with open(tmp_path / f"{method}_surrogate.json", "w") as f:
            json.dump(data, f)

    return str(tmp_path) + os.sep


def test_performance_estimation_poly(surrogate_path):
    result = performance_estimation(method="poly", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == len(dummy_components)

    assert list(result.columns) == ["MAE", "MSE", "R2", "Adjusted R2", "Comp"]

    assert list(result["Comp"]) == dummy_components


def test_performance_estimation_kri(surrogate_path):
    result = performance_estimation(method="kri", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == len(dummy_components)

    assert list(result.columns) == ["R2", "RMSE", "Comp"]

    assert list(result["Comp"]) == dummy_components


def test_performance_estimation_rbf(surrogate_path):
    result = performance_estimation(method="rbf", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == len(dummy_components)

    assert list(result.columns) == ["R2", "RMSE", "Comp"]

    assert list(result["Comp"]) == dummy_components


def test_performance_estimation_file_not_found():
    with pytest.raises(FileNotFoundError):
        performance_estimation(method="poly", path="./file_not_found/")


def test_display_performance_poly(surrogate_path):
    result = display_performance(method="poly", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == [
        "Predicted Variables",
        "R^2",
        "Adjusted R^2",
        "MAE",
        "MSE",
    ]
    assert result.shape[0] == len(dummy_components)


def test_display_performance_kri(surrogate_path):
    result = display_performance(method="kri", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["R^2", "RMSE"]
    assert result.shape[0] == len(dummy_components)


def test_display_performance_rbf(surrogate_path):
    result = display_performance(method="rbf", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["R^2", "RMSE"]
    assert result.shape[0] == len(dummy_components)


def test_display_performance_invalid_method(surrogate_path):
    with pytest.raises(ValueError, match="Unsupported method"):
        display_performance(method="invalid", path=surrogate_path)


def test_display_plot_w_path():
    if IPython is None:
        pytest.skip("IPython not available")

    path = os.path.abspath(os.path.join(local_path, "..", "results"))
    result = display_plot(method="poly", path=path)

    assert hasattr(result, "src")
    assert "poly_parity.pdf" in result.src
    assert path in result.src


def test_display_plot_wo_path():
    if IPython is None:
        pytest.skip("IPython not available")

    result = display_plot(method="poly", path=None)

    assert hasattr(result, "src")
    assert "poly_parity.pdf" in result.src
