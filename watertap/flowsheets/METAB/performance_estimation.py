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
import json
import pandas as pd
from IPython.display import IFrame

# third-party
try:
    import pdf2image
    from pdf2image import convert_from_path
except ImportError:
    pdf2image = None
from IPython.display import display

local_path = os.path.dirname(os.path.abspath(__file__))


def performance_estimation(
    method="poly",  # "rbf"#"kri"alamo'
    path="./results/",
):

    if method not in ("poly", "kri", "rbf"):
        raise ValueError(
            f"Unsupported method: {method}. Choose from 'poly', 'kri', or 'rbf'."
        )

    metrics_sum = pd.DataFrame()
    file = path + method + "_surrogate.json"

    with open(file, "r") as file:
        data = json.load(file)

    if method == "poly":
        for ele in data["model_encoding"]:
            metrics = data["model_encoding"][ele]["attr"]["errors"]
            metrics["Comp"] = ele
            # print(metrics)
            for key in metrics:
                metrics[key] = [metrics[key]]
            # print(metrics)
            df = pd.DataFrame.from_dict(metrics)
            metrics_sum = pd.concat([metrics_sum, df])

    elif method == "kri":
        for ele in data["model_encoding"]:
            metrics = {}
            metrics["R2"] = [data["model_encoding"][ele]["attr"]["training_R2"]]
            metrics["RMSE"] = [data["model_encoding"][ele]["attr"]["training_rmse"]]
            metrics["Comp"] = [ele]
            df = pd.DataFrame.from_dict(metrics)
            metrics_sum = pd.concat([metrics_sum, df])

    elif method == "rbf":
        for ele in data["model_encoding"]:
            metrics = {}
            metrics["R2"] = [data["model_encoding"][ele]["attr"]["R2"]]
            metrics["RMSE"] = [data["model_encoding"][ele]["attr"]["rmse"]]
            metrics["Comp"] = [ele]
            df = pd.DataFrame.from_dict(metrics)
            metrics_sum = pd.concat([metrics_sum, df])

    return metrics_sum


def display_performance(method="poly", path="./results/"):
    metrics = performance_estimation(method=method, path=path)
    display_metrics = pd.DataFrame()

    if method == "poly":
        display_metrics["Predicted Variables"] = metrics["Comp"]
        display_metrics["R^2"] = metrics["R2"]
        display_metrics["Adjusted R^2"] = metrics["Adjusted R2"]
        display_metrics["MAE"] = metrics["MAE"]
        display_metrics["MSE"] = metrics["MSE"]
    elif method in ("kri", "rbf"):
        display_metrics["R^2"] = metrics["R2"]
        display_metrics["RMSE"] = metrics["RMSE"]

    display_metrics.index = range(1, len(display_metrics) + 1)

    return display_metrics


def display_plot(method="poly", path=None):
    if path is None:
        path = os.path.join(local_path, "results")

    file_path = os.path.join(path, "{}_parity.pdf".format(method))
    return IFrame(file_path, width=700, height=500)
