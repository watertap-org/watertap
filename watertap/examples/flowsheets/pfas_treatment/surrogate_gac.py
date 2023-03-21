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
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)

__author__ = "Hunter Barber"


def main():

    min_st_regression()
    throughput_regression()


def min_st_regression():
    # ---------------------------------------------------------------------
    # minimum stanton number equation parameter lookup

    min_st_ninv = list(np.linspace(0.1, 0.9, 9))
    min_st_ninv.insert(0, 0.05)
    min_st_bi = list(np.linspace(0.5, 10, 20)) + list(np.linspace(10, 100, 19))

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

    # ---------------------------------------------------------------------
    # minimum stanton number surrogate training

    trainer = surrogate.PysmoKrigingTrainer(
        input_labels=["ninv", "bi"],
        output_labels=["min_st"],
        training_dataframe=min_st_df,
    )
    pysmo_surr_expr = trainer.train_surrogate()

    input_labels = trainer._input_labels
    output_labels = trainer._output_labels
    xmin, xmax = [0, 0.5], [1, 100]
    input_bounds = {
        input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
    }

    show = False
    pysmo_surr = surrogate.PysmoSurrogate(
        pysmo_surr_expr, input_labels, output_labels, input_bounds
    )
    surrogate_scatter2D(
        pysmo_surr,
        min_st_df,
        filename="watertap/examples/flowsheets/pfas_treatment/min_st_scatter2D.pdf",
        show=show,
    )
    surrogate_scatter2D(
        pysmo_surr,
        min_st_df,
        filename="watertap/examples/flowsheets/pfas_treatment/min_st_scatter3D.pdf",
        show=show,
    )
    surrogate_parity(
        pysmo_surr,
        min_st_df,
        filename="watertap/examples/flowsheets/pfas_treatment/min_st_parity.pdf",
        show=show,
    )
    surrogate_residual(
        pysmo_surr,
        min_st_df,
        filename="watertap/examples/flowsheets/pfas_treatment/min_st_residual.pdf",
        show=show,
    )

    # To save a model
    model = pysmo_surr.save_to_file(
        "watertap/examples/flowsheets/pfas_treatment/min_st_surrogate.json",
        overwrite=True,
    )


def throughput_regression():
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

    # ---------------------------------------------------------------------
    # throughput equation surrogate training

    trainer = surrogate.PysmoKrigingTrainer(
        input_labels=["ninv", "bi", "conc_ratio"],
        output_labels=["throughput"],
        training_dataframe=throughput_df,
    )
    pysmo_surr_expr = trainer.train_surrogate()

    input_labels = trainer._input_labels
    output_labels = trainer._output_labels
    xmin, xmax = [0, 0.5, 0.01], [1, 100, 0.99]
    input_bounds = {
        input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
    }

    show = False
    pysmo_surr = surrogate.PysmoSurrogate(
        pysmo_surr_expr, input_labels, output_labels, input_bounds
    )
    surrogate_scatter2D(
        pysmo_surr,
        throughput_df,
        filename="watertap/examples/flowsheets/pfas_treatment/throughput_scatter2D.pdf",
        show=show,
    )
    surrogate_scatter2D(
        pysmo_surr,
        throughput_df,
        filename="watertap/examples/flowsheets/pfas_treatment/throughput_scatter3D.pdf",
        show=show,
    )
    surrogate_parity(
        pysmo_surr,
        throughput_df,
        filename="watertap/examples/flowsheets/pfas_treatment/throughput_parity.pdf",
        show=show,
    )
    surrogate_residual(
        pysmo_surr,
        throughput_df,
        filename="watertap/examples/flowsheets/pfas_treatment/throughput_residual.pdf",
        show=show,
    )

    # To save a model
    model = pysmo_surr.save_to_file(
        "watertap/examples/flowsheets/pfas_treatment/throughput_surrogate.json",
        overwrite=True,
    )


if __name__ == "__main__":
    main()
