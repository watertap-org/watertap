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
import numpy as np
import pandas as pd

# Import Pyomo libraries
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    Var,
    Constraint,
    Set,
    Objective,
    maximize,
)
from pyomo.common.timing import TicTocTimer

# Import IDAES libraries
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer, PysmoSurrogate
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock

from idaes.core.surrogate.pysmo_surrogate import (
    PysmoPolyTrainer,
    PysmoRBFTrainer,
    PysmoKrigingTrainer,
    PysmoSurrogate,
)
from idaes.core.surrogate.alamopy import AlamoTrainer, AlamoSurrogate


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


def outputs_selections(output_data):
    outputs_list = [
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

    # output_data = output_data[(output_data >= 0).all(axis=1)]
    print(output_data)
    output_data.columns = output_data.columns.str.replace(" ", "")
    output_data = output_data[outputs_list]

    print(output_data)

    return output_data


def gen_surrogate_model(
    tool="idaes", method="poly", feed_data=None, input_data=None, output_data=None
):
    if feed_data == None:
        feed_data = pd.concat([input_data, output_data], axis=1)

    input_labels = list(input_data.columns)
    output_labels = list(output_data.columns)
    xmin, xmax = input_data.min().tolist(), input_data.max().tolist()
    input_bounds = {
        input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
    }
    if method == "poly":
        # Create PySMO trainer object
        trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=feed_data,
        )

        # Set PySMO options
        trainer.config.maximum_polynomial_order = 6
        trainer.config.multinomials = True
        trainer.config.training_split = 0.8
        trainer.config.number_of_crossvalidations = 3
        # Train surrogate (calls PySMO through IDAES Python wrapper)
        poly_train = trainer.train_surrogate()

        poly_surr = PysmoSurrogate(
            poly_train, input_labels, output_labels, input_bounds
        )

        poly_surr.save_to_file(
            "./results/{}_surrogate.json".format(method), overwrite=True
        )

        surrogate_parity(
            poly_surr, feed_data, filename="./results/{}_parity.pdf".format(method)
        )

    elif method == "kri":
        trainer = PysmoKrigingTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=feed_data,
        )
        # Set PySMO options
        trainer.config.numerical_gradients = True
        trainer.config.regularization = True
        # Train surrogate (calls PySMO through IDAES Python wrapper)
        krig_train = trainer.train_surrogate()

        krig_surr = PysmoSurrogate(
            krig_train, input_labels, output_labels, input_bounds
        )
        krig_surr.save_to_file(
            "./results/{}_surrogate.json".format(method), overwrite=True
        )
        surrogate_parity(
            krig_surr, feed_data, filename="./results/{}_parity.pdf".format(method)
        )

    elif method == "rbf":
        trainer = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=feed_data,
        )

        trainer.config.basis_function = "cubic"
        rbf_train = trainer.train_surrogate()
        rbf_surr = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds)
        rbf_surr.save_to_file(
            "./results/{}_surrogate.json".format(method), overwrite=True
        )
        surrogate_parity(
            rbf_surr, feed_data, filename="./results/{}_parity.pdf".format(method)
        )

    elif method == "alamo":
        trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=feed_data,
        )
        trainer.config.constant = 1
        trainer.config.linfcns = 1
        trainer.config.expfcns = 1
        trainer.config.logfcns = 1
        trainer.config.sinfcns = 1
        trainer.config.cosfcns = 1
        trainer.config.monomialpower = [2, 3, 4]
        trainer.config.multi2power = [1, 2, 3]
        trainer.config.multi3power = [1, 2, 3]
        success, alm_surr, msg = trainer.train_surrogate()
        surrogate_expressions = trainer._results["Model"]
        alm_surr = AlamoSurrogate(
            surrogate_expressions, input_labels, output_labels, input_bounds
        )
        alm_surr.save_to_file(
            "./results/{}_surrogate.json".format(method), overwrite=True
        )
        surrogate_parity(
            alm_surr, feed_data, filename="./results/{}_parity.pdf".format(method)
        )


if __name__ == "__main__":

    feed_data, input_data, output_data = get_data()
    output_data = outputs_selections(output_data)
    gen_surrogate_model(
        tool="idaes",
        method="rbf",  # kri, poly,rbf,alamo
        feed_data=feed_data,
        input_data=input_data,
        output_data=output_data,
    )
