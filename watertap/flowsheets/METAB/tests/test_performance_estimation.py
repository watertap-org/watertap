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

from watertap.flowsheets.METAB.performance_estimation import (
    performance_estimation,
    display_performance,
    display_plot,
)

local_path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def surrogate_path():
    return os.path.abspath(os.path.join(local_path, "..", "results")) + os.sep


def test_performance_estimation_poly(surrogate_path):
    result = performance_estimation(method="poly", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 25

    assert list(result.columns) == ["MAE", "MSE", "R2", "Adjusted R2", "Comp"]

    expected_components = [
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
    assert list(result["Comp"]) == expected_components


def test_performance_estimation_kri(surrogate_path):
    result = performance_estimation(method="kri", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 25

    assert list(result.columns) == ["R2", "RMSE", "Comp"]

    expected_components = [
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
    assert list(result["Comp"]) == expected_components


def test_performance_estimation_rbf(surrogate_path):
    result = performance_estimation(method="rbf", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 25

    assert list(result.columns) == ["R2", "RMSE", "Comp"]

    expected_components = [
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
    assert list(result["Comp"]) == expected_components


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
    assert result.shape[0] == 25


def test_display_performance_kri(surrogate_path):
    result = display_performance(method="kri", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["R^2", "RMSE"]
    assert result.shape[0] == 25


def test_display_performance_rbf(surrogate_path):
    result = display_performance(method="rbf", path=surrogate_path)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["R^2", "RMSE"]
    assert result.shape[0] == 25

    # TODO: Need to add tests for surrogate_model_generator, and separate tests into unique files


def test_display_performance_invalid_method(surrogate_path):
    with pytest.raises(ValueError, match="Unsupported method"):
        display_performance(method="invalid", path=surrogate_path)


def test_display_plot_w_path():
    path = os.path.abspath(os.path.join(local_path, "..", "results"))
    result = display_plot(method="poly", path=path)

    assert hasattr(result, "src")
    assert "poly_parity.pdf" in result.src
    assert path in result.src


def test_display_plot_wo_path():
    result = display_plot(method="poly", path=None)

    assert hasattr(result, "src")
    assert "poly_parity.pdf" in result.src
