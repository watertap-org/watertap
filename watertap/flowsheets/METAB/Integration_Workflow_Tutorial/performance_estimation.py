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
