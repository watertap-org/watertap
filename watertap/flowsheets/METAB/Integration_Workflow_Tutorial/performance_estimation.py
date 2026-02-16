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
import json
import pandas as pd
import PyPDF2
from pdf2image import convert_from_path
from IPython.display import display


def get_data(
    tool="pysmo",
    input_data_file="./results/input_data.csv",
    output_data_file="./results/output_data.csv",
):
    input_data = pd.read_csv(input_data_file, header=0)
    # print(input_data)
    output_data = pd.read_csv(output_data_file, header=0).iloc[:, 1:]
    # print(output_data)
    feed_data = pd.concat([input_data, output_data], axis=1)
    return feed_data, input_data, output_data


def performance_estimiation(
    method="poly",  # "rbf"#"kri"alamo'
    path="./results/",
):
    method_list = []
    file = path + method + "_surrogate.json"

    with open(file, "r") as file:
        data = json.load(file)

    if method == "poly":
        metrics_sum = pd.DataFrame()
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
        metrics_sum = pd.DataFrame()
        for ele in data["model_encoding"]:
            metrics = {}
            metrics["R2"] = [data["model_encoding"][ele]["attr"]["training_R2"]]
            metrics["RMSE"] = [data["model_encoding"][ele]["attr"]["training_rmse"]]
            metrics["Comp"] = [ele]
            df = pd.DataFrame.from_dict(metrics)
            metrics_sum = pd.concat([metrics_sum, df])

    elif method == "rbf":
        metrics_sum = pd.DataFrame()
        for ele in data["model_encoding"]:
            metrics = {}
            metrics["R2"] = [data["model_encoding"][ele]["attr"]["R2"]]
            metrics["RMSE"] = [data["model_encoding"][ele]["attr"]["rmse"]]
            metrics["Comp"] = [ele]
            df = pd.DataFrame.from_dict(metrics)
            metrics_sum = pd.concat([metrics_sum, df])

    return metrics_sum


def display_performace(method="poly", path="./results/"):
    metrics = performance_estimiation(method=method, path=path)
    if method == "poly":
        display_metrics = pd.DataFrame()
        display_metrics["Predicted Variables"] = metrics["Comp"]
        display_metrics["R^2"] = metrics["R2"]
        display_metrics["Adjusted R^2"] = metrics["Adjusted R2"]
        display_metrics["MAE"] = metrics["MAE"]
        display_metrics["MSE"] = metrics["MSE"]
        display_metrics.index = range(1, len(display_metrics) + 1)
    return display_metrics


def display_plot(method="poly", path="./results/"):

    file_path = path + method + "_parity.pdf"

    # Convert PDF pages to images
    pages = convert_from_path(file_path)

    # Get the first page as an image
    first_page_image = pages[0]

    # Display in Jupyter Notebook
    display(first_page_image)

    # NEED poppler for Windows
