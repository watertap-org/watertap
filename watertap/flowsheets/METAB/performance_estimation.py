#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import os
import pandas as pd

from idaes.core.surrogate.metrics import compute_fit_metrics
from idaes.core.surrogate.plotting.sm_plotter import surrogate_parity

local_path = os.path.dirname(os.path.abspath(__file__))

__author__ = "Marcus Holly"


def performance_estimation(
    surrogate=None,
    dataframe=None,
):
    """
    Evaluates the surrogate against the provided dataframe and returns a
    DataFrame summarising fit metrics for each output variable.

    Args:
        surrogate: IDAES surrogate object
            Must implement input_labels(), output_labels(), and
            evaluate_surrogate()
        dataframe (pd.DataFrame): DataFrame containing both input and output
            columns corresponding to the surrogate's input_labels() and
            output_labels()

    Returns:
        pd.DataFrame: One row per output variable with columns:
            - ``Comp``  : output variable name
            - ``R2``    : coefficient of determination
            - ``RMSE``  : root mean squared error
            - ``MSE``   : mean squared error
            - ``MAE``   : mean absolute error
            - ``maxAE`` : maximum absolute error
            - ``SSE``   : sum of squared errors
    """

    if surrogate is None or dataframe is None:
        raise ValueError(
            "Both 'surrogate' and 'dataframe' must be provided to use compute_fit_metrics."
        )

    # Delegate all metric computation to IDAES
    metrics_by_output = compute_fit_metrics(surrogate, dataframe)

    # Reshape dict-of-dicts into a flat DataFrame, one row per output label
    rows = []
    for output_label, metrics in metrics_by_output.items():
        row = {"Comp": output_label, **metrics}
        rows.append(row)

    metrics_sum = pd.DataFrame(rows)

    return metrics_sum


def display_performance(surrogate=None, dataframe=None):
    metrics = performance_estimation(
        surrogate=surrogate,
        dataframe=dataframe,
    )
    """
    Return a formatted DataFrame of surrogate performance metrics for display.

    Args:
        surrogate: IDAES surrogate object
        dataframe: DataFrame containing input and output columns

    Returns:
        pd.DataFrame: One row per output variable, 1-based index, with columns:
            - ``Predicted Variables`` : output variable name
            - ``R^2``                 : coefficient of determination
            - ``MAE``                 : mean absolute error
            - ``MSE``                 : mean squared error
            - ``RMSE``                : root mean squared error
            - ``Max AE``              : maximum absolute error
    """

    display_metrics = pd.DataFrame(
        {
            "Predicted Variables": metrics["Comp"],
            "R^2": metrics["R2"],
            "MAE": metrics["MAE"],
            "MSE": metrics["MSE"],
            "RMSE": metrics["RMSE"],
            "Max AE": metrics["maxAE"],
        }
    )

    display_metrics.index = range(1, len(display_metrics) + 1)

    return display_metrics


def display_plot(surrogate, dataframe, method="poly", path=None, show=True):
    """
    Produces one parity plot per output variable and consolidates into a PDF

    Args:
        surrogate: IDAES surrogate object
        dataframe: DataFrame containing input and output
            columns corresponding to the surrogate's labels
        method (str, optional): Surrogate method identifier used to name the
            output PDF file (e.g. ``"poly"``, ``"kri"``, ``"rbf"``)
        path (str, optional): Directory in which to save the PDF
        show (bool, optional): Whether to display each figure interactively
            via matplotlib

    Returns:
        list[matplotlib.figure.Figure]: One figure per output variable
    """
    if path is None:
        path = os.path.join(local_path, "results")

    filename = os.path.join(path, f"{method}_parity.pdf")

    return surrogate_parity(surrogate, dataframe, filename=filename, show=show)
