#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import numpy as np
import pandas as pd
import idaes.core.surrogate.pysmo_surrogate as surrogate
from idaes.core.surrogate.pysmo import radial_basis_function
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)
from idaes.core.surrogate.metrics import compute_fit_metrics
import matplotlib.pyplot as plt

__author__ = "Hunter Barber"


def main():

    rerun = False
    min_st_regression(rerun=rerun)
    throughput_regression(rerun=rerun)


def min_st_regression(rerun=False):

    # ---------------------------------------------------------------------
    # minimum stanton number equation parameter lookup

    min_st_ninv = list(np.linspace(0.1, 0.9, 9))
    min_st_ninv.insert(0, 0.05)
    min_st_bi = (
        list(np.linspace(0.5, 10, 20))
        + list(np.linspace(15, 100, 18))
        + list(np.linspace(150, 500, 8))
    )

    ninv_list, bi_list, min_st_list = [], [], []

    min_st_param_df = pd.read_csv(
        "watertap/examples/flowsheets/pfas_treatment/min_st_parameters.csv",
        index_col=["1/n", "Bi_min", "Bi_max", "% within constant pattern"],
    )
    ninv_param_list = min_st_param_df.index.get_level_values(0).values

    for ninv in min_st_ninv:
        for bi in min_st_bi:

            ninv_param_closest = ninv_param_list[
                min(
                    range(len(ninv_param_list)),
                    key=lambda i: abs(ninv_param_list[i] - ninv),
                )
            ]

            if bi < 10:
                lookup = min_st_param_df.loc[ninv_param_closest, 0.5, 10, 0]
            else:
                lookup = min_st_param_df.loc[ninv_param_closest, 10, float("NaN"), 0]

            a0_lookup = lookup[0]
            a1_lookup = lookup[1]

            ninv_list.append(ninv)
            bi_list.append(bi)
            min_st = a0_lookup * bi + a1_lookup
            min_st_list.append(min_st)

    min_st_df = pd.DataFrame(
        {
            "ninv": ninv_list,
            "bi": bi_list,
            "min_st": min_st_list,
        }
    )

    if rerun:

        # ---------------------------------------------------------------------
        # minimum stanton number surrogate trainging

        trainer = surrogate.PysmoKrigingTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        pysmo_surr_expr = trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_kriging = surrogate.PysmoSurrogate(
            pysmo_surr_expr, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_kriging.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_kriging.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "linear"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_linear = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_linear.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_linear.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "cubic"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_cubic = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_cubic.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_cubic.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "spline"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_spline = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_spline.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_spline.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "gaussian"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_gaussian = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_gaussian.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_gaussian.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "mq"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_mq = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_mq.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_mq.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi"],
            output_labels=["min_st"],
            training_dataframe=min_st_df,
        )
        rbf_trainer.config.basis_function = "imq"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5], [1, 500]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        min_st_pysmo_surr_imq = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = min_st_pysmo_surr_imq.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_imq.json",
            overwrite=True,
        )

    else:

        # ---------------------------------------------------------------------

        min_st_pysmo_surr_kriging = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_kriging.json",
        )
        min_st_pysmo_surr_linear = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_linear.json",
        )
        min_st_pysmo_surr_cubic = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_cubic.json",
        )
        min_st_pysmo_surr_spline = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_spline.json",
        )
        min_st_pysmo_surr_gaussian = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_gaussian.json",
        )
        min_st_pysmo_surr_mq = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_mq.json",
        )
        min_st_pysmo_surr_imq = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_imq.json",
        )

        # ---------------------------------------------------------------------

        surrogate_models = [
            min_st_pysmo_surr_kriging,
            min_st_pysmo_surr_linear,
            min_st_pysmo_surr_cubic,
            min_st_pysmo_surr_spline,
            min_st_pysmo_surr_gaussian,
            min_st_pysmo_surr_mq,
            min_st_pysmo_surr_imq,
        ]
        res_names = []
        for m in range(0, len(surrogate_models)):
            err = compute_fit_metrics(surrogate_models[m], min_st_df)
            err = pd.DataFrame.from_dict(err)
            print(err)
            res_names.append(err)

        # Plot metrics
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["R2"],
                "Linear-RBF": res_names[1].loc["R2"],
                "Cubic-RBF": res_names[2].loc["R2"],
                "Spline-RBF": res_names[3].loc["R2"],
                "Gaussian-RBF": res_names[4].loc["R2"],
                "MQ-RBF": res_names[4].loc["R2"],
                "IMQ-RBF": res_names[4].loc["R2"],
            }
        )
        df.plot.bar(rot=0, ax=axes[0, 0])
        axes[0, 0].set_title("R2")
        axes[0, 0].set_ylabel("R2")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["RMSE"],
                "Linear-RBF": res_names[1].loc["RMSE"],
                "Cubic-RBF": res_names[2].loc["RMSE"],
                "Spline-RBF": res_names[3].loc["RMSE"],
                "Gaussian-RBF": res_names[4].loc["RMSE"],
                "MQ-RBF": res_names[4].loc["RMSE"],
                "IMQ-RBF": res_names[4].loc["RMSE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[0, 1], logy=True)
        axes[0, 1].set_title("RMSE")
        axes[0, 1].set_ylabel("Log(RMSE)")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["MAE"],
                "Linear-RBF": res_names[1].loc["MAE"],
                "Cubic-RBF": res_names[2].loc["MAE"],
                "Spline-RBF": res_names[3].loc["MAE"],
                "Gaussian-RBF": res_names[4].loc["MAE"],
                "MQ-RBF": res_names[4].loc["MAE"],
                "IMQ-RBF": res_names[4].loc["MAE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[1, 0], logy=True)
        axes[1, 0].set_title("MAE")
        axes[1, 0].set_ylabel("Log(MAE)")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["maxAE"],
                "Linear-RBF": res_names[1].loc["maxAE"],
                "Cubic-RBF": res_names[2].loc["maxAE"],
                "Spline-RBF": res_names[3].loc["maxAE"],
                "Gaussian-RBF": res_names[4].loc["maxAE"],
                "MQ-RBF": res_names[4].loc["maxAE"],
                "IMQ-RBF": res_names[4].loc["maxAE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[1, 1], logy=True)
        axes[1, 1].set_ylabel("Log(maxAE)")
        axes[1, 1].set_title("maxAE")


def throughput_regression(rerun=False):

    # ---------------------------------------------------------------------
    # throughput equation parameter lookup

    throughput_df = pd.read_csv(
        "watertap/examples/flowsheets/pfas_treatment/throughput_parameters.csv",
    )
    throughput_ninv_bi = throughput_df.loc[:, ["1/n", "Bi"]]
    throughput_conc_ratio = list(np.linspace(0.05, 0.95, 19))
    throughput_conc_ratio.insert(0, 0.01)
    throughput_conc_ratio.insert(-1, 0.99)

    ninv_list, bi_list, conc_ratio_list, throughput_list = [], [], [], []

    throughput_param_df = pd.read_csv(
        "watertap/examples/flowsheets/pfas_treatment/throughput_parameters.csv",
        index_col=["1/n", "Bi"],
    )
    ninv_param_list = throughput_param_df.index.get_level_values(0).values

    for ind, row in throughput_ninv_bi.iterrows():
        ninv = row["1/n"]
        bi = row["Bi"]
        for conc_ratio in throughput_conc_ratio:

            ninv_param_closest = ninv_param_list[
                min(
                    range(len(ninv_param_list)),
                    key=lambda i: abs(ninv_param_list[i] - ninv),
                )
            ]
            N_Bi_param_list = throughput_param_df.loc[
                ninv_param_closest, :
            ].index.values
            N_Bi_param_closest = N_Bi_param_list[
                min(
                    range(len(N_Bi_param_list)),
                    key=lambda i: abs(N_Bi_param_list[i] - bi),
                )
            ]
            lookup = throughput_param_df.loc[ninv_param_closest, N_Bi_param_closest]
            b0_lookup = lookup[0]
            b1_lookup = lookup[1]
            b2_lookup = lookup[2]
            b3_lookup = lookup[3]
            b4_lookup = lookup[4]

            ninv_list.append(ninv)
            bi_list.append(bi)
            conc_ratio_list.append(conc_ratio)
            throughput = (
                b0_lookup
                + b1_lookup * (conc_ratio**b2_lookup)
                + b3_lookup / (1.01 - (conc_ratio**b4_lookup))
            )
            throughput_list.append(throughput)

    throughput_df = pd.DataFrame(
        {
            "ninv": ninv_list,
            "bi": bi_list,
            "conc_ratio": conc_ratio_list,
            "throughput": throughput_list,
        }
    )

    if rerun:

        # ---------------------------------------------------------------------
        # throughput surrogate trainging

        trainer = surrogate.PysmoKrigingTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        pysmo_surr_expr = trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_kriging = surrogate.PysmoSurrogate(
            pysmo_surr_expr, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_kriging.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_kriging.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "linear"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_linear = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_linear.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_linear.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "cubic"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_cubic = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_cubic.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_cubic.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "spline"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_spline = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_spline.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_spline.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "gaussian"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_gaussian = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_gaussian.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_gaussian.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "mq"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_mq = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_mq.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_mq.json",
            overwrite=True,
        )

        # ---------------------------------------------------------------------

        rbf_trainer = surrogate.PysmoRBFTrainer(
            input_labels=["ninv", "bi", "conc_ratio"],
            output_labels=["throughput"],
            training_dataframe=throughput_df,
        )
        rbf_trainer.config.basis_function = "imq"
        rbf_train = rbf_trainer.train_surrogate()

        input_labels = trainer._input_labels
        output_labels = trainer._output_labels
        xmin, xmax = [0, 0.5, 0.001], [1, 500, 0.999]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }

        throughput_pysmo_surr_imq = surrogate.PysmoSurrogate(
            rbf_train, input_labels, output_labels, input_bounds
        )

        # To save a model
        model = throughput_pysmo_surr_imq.save_to_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_imq.json",
            overwrite=True,
        )

    else:

        # ---------------------------------------------------------------------

        throughput_pysmo_surr_kriging = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_kriging.json",
        )
        throughput_pysmo_surr_linear = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_linear.json",
        )
        throughput_pysmo_surr_cubic = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_cubic.json",
        )
        throughput_pysmo_surr_spline = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_spline.json",
        )
        throughput_pysmo_surr_gaussian = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_gaussian.json",
        )
        throughput_pysmo_surr_mq = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_mq.json",
        )
        throughput_pysmo_surr_imq = surrogate.PysmoSurrogate.load_from_file(
            "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_imq.json",
        )

        # ---------------------------------------------------------------------

        surrogate_models = [
            throughput_pysmo_surr_kriging,
            throughput_pysmo_surr_linear,
            throughput_pysmo_surr_cubic,
            throughput_pysmo_surr_spline,
            throughput_pysmo_surr_gaussian,
            throughput_pysmo_surr_mq,
            throughput_pysmo_surr_imq,
        ]
        res_names = []
        for m in range(0, len(surrogate_models)):
            err = compute_fit_metrics(surrogate_models[m], throughput_df)
            err = pd.DataFrame.from_dict(err)
            print(err)
            res_names.append(err)

        # Plot metrics
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["R2"],
                "Linear-RBF": res_names[1].loc["R2"],
                "Cubic-RBF": res_names[2].loc["R2"],
                "Spline-RBF": res_names[3].loc["R2"],
                "Gaussian-RBF": res_names[4].loc["R2"],
                "MQ-RBF": res_names[4].loc["R2"],
                "IMQ-RBF": res_names[4].loc["R2"],
            }
        )
        df.plot.bar(rot=0, ax=axes[0, 0])
        axes[0, 0].set_title("R2")
        axes[0, 0].set_ylabel("R2")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["RMSE"],
                "Linear-RBF": res_names[1].loc["RMSE"],
                "Cubic-RBF": res_names[2].loc["RMSE"],
                "Spline-RBF": res_names[3].loc["RMSE"],
                "Gaussian-RBF": res_names[4].loc["RMSE"],
                "MQ-RBF": res_names[4].loc["RMSE"],
                "IMQ-RBF": res_names[4].loc["RMSE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[0, 1], logy=True)
        axes[0, 1].set_title("RMSE")
        axes[0, 1].set_ylabel("Log(RMSE)")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["MAE"],
                "Linear-RBF": res_names[1].loc["MAE"],
                "Cubic-RBF": res_names[2].loc["MAE"],
                "Spline-RBF": res_names[3].loc["MAE"],
                "Gaussian-RBF": res_names[4].loc["MAE"],
                "MQ-RBF": res_names[4].loc["MAE"],
                "IMQ-RBF": res_names[4].loc["MAE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[1, 0], logy=True)
        axes[1, 0].set_title("MAE")
        axes[1, 0].set_ylabel("Log(MAE)")

        df = pd.DataFrame(
            {
                "Kriging": res_names[0].loc["maxAE"],
                "Linear-RBF": res_names[1].loc["maxAE"],
                "Cubic-RBF": res_names[2].loc["maxAE"],
                "Spline-RBF": res_names[3].loc["maxAE"],
                "Gaussian-RBF": res_names[4].loc["maxAE"],
                "MQ-RBF": res_names[4].loc["maxAE"],
                "IMQ-RBF": res_names[4].loc["maxAE"],
            }
        )
        df.plot.bar(rot=0, ax=axes[1, 1], logy=True)
        axes[1, 1].set_ylabel("Log(maxAE)")
        axes[1, 1].set_title("maxAE")

        plt.show()

        """
        show = True
        surrogate_scatter2D(
            throughput_pysmo_surr_linear,
            throughput_df,
            filename=None,
            show=show,
        )
        surrogate_scatter2D(
            throughput_pysmo_surr_linear,
            throughput_df,
            filename=None,
            show=show,
        )
        surrogate_parity(
            throughput_pysmo_surr_linear,
            throughput_df,
            filename=None,
            show=show,
        )
        surrogate_residual(
            throughput_pysmo_surr_linear,
            throughput_df,
            filename=None,
            show=show,
        )
        """


if __name__ == "__main__":
    main()
