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

import pyomo.environ as pyo
from idaes.core.surrogate import AlamoTrainer, AlamoSurrogate
from idaes.core.util.model_statistics import (
    variables_near_bounds_generator,
    total_constraints_set,
)
from idaes.core.util.scaling import (
    badly_scaled_var_generator,
)
import pandas as pd
from pyomo.environ import (
    Var,
    Param,
    Set,
    Constraint,
    Piecewise,
    Suffix,
    NonNegativeReals,
    PositiveIntegers,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In
from enum import Enum, auto

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import idaes.core.surrogate.pysmo_surrogate as surrogate
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)
import os

from watertap.core import ControlVolume0DBlock, InitializationMixin

__author__ = "Hunter Barber"

_log = idaeslog.getLogger(__name__)


freund_ninv = 0.564
N_Bi = 23.4

# minimum stanton number equation parameter lookup
min_st_param_df = pd.read_csv(
    "watertap/examples/flowsheets/pfas_treatment/min_st_parameters.csv",
    index_col=["1/n", "Bi_min", "Bi_max", "% within constant pattern"],
)
print(min_st_param_df)
freund_ninv_param_list = min_st_param_df.index.get_level_values(0).values
freund_ninv_param_closest = freund_ninv_param_list[
    min(
        range(len(freund_ninv_param_list)),
        key=lambda i: abs(freund_ninv_param_list[i] - freund_ninv),
    )
]
if N_Bi < 10:
    lookup = min_st_param_df.loc[freund_ninv_param_closest, 0.5, 10, 0]
else:
    lookup = min_st_param_df.loc[freund_ninv_param_closest, 10, float("NaN"), 0]
a0_lookup = lookup[0]
a1_lookup = lookup[1]

print("a0_lookup", a0_lookup)
print("a1_lookup", a1_lookup)

# throughput equation parameter lookup
throughput_param_df = pd.read_csv(
    "watertap/examples/flowsheets/pfas_treatment/throughput_parameters.csv",
    index_col=["1/n", "Bi"],
)
print(throughput_param_df)
freund_ninv_param_list = throughput_param_df.index.get_level_values(0).values
freund_ninv_param_closest = freund_ninv_param_list[
    min(
        range(len(freund_ninv_param_list)),
        key=lambda i: abs(freund_ninv_param_list[i] - freund_ninv),
    )
]
N_Bi_param_list = throughput_param_df.loc[freund_ninv_param_closest, :].index.values
N_Bi_param_closest = N_Bi_param_list[
    min(
        range(len(N_Bi_param_list)),
        key=lambda i: abs(N_Bi_param_list[i] - N_Bi),
    )
]
lookup = throughput_param_df.loc[freund_ninv_param_closest, N_Bi_param_closest]
b0_lookup = lookup[0]
b1_lookup = lookup[1]
b2_lookup = lookup[2]
b3_lookup = lookup[3]
b4_lookup = lookup[4]

print("b0_lookup", b0_lookup)
print("b1_lookup", b1_lookup)
print("b2_lookup", b2_lookup)
print("b3_lookup", b3_lookup)
print("b4_lookup", b4_lookup)

# trainer min stanton
throughput_param_df = pd.read_csv(
    "watertap/examples/flowsheets/pfas_treatment/throughput_parameters.csv",
)
# after reading or generating a DataFrame object called `data_training`
trainer = surrogate.PysmoKrigingTrainer(
    input_labels=["1/n", "Bi"],
    output_labels=["B0", "B1", "B2", "B3", "B4"],
    training_dataframe=throughput_param_df,
)
pysmo_surr_expr = trainer.train_surrogate()

input_labels = trainer._input_labels
output_labels = trainer._output_labels
xmin, xmax = [0, 0.5], [1, 100]
input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}

pysmo_surr = surrogate.PysmoSurrogate(
    pysmo_surr_expr, input_labels, output_labels, input_bounds
)

surrogate_scatter2D(
    pysmo_surr, throughput_param_df, filename="scatter2D.pdf", show=True
)
surrogate_scatter2D(
    pysmo_surr, throughput_param_df, filename="scatter3D.pdf", show=True
)
surrogate_parity(pysmo_surr, throughput_param_df, filename="parity.pdf", show=True)
surrogate_residual(pysmo_surr, throughput_param_df, filename="residual.pdf", show=True)

# trainer throughput
throughput_param_df = pd.read_csv(
    "watertap/examples/flowsheets/pfas_treatment/throughput_parameters.csv",
)
# after reading or generating a DataFrame object called `data_training`
trainer = surrogate.PysmoKrigingTrainer(
    input_labels=["1/n", "Bi"],
    output_labels=["B0", "B1", "B2", "B3", "B4"],
    training_dataframe=throughput_param_df,
)
pysmo_surr_expr = trainer.train_surrogate()

input_labels = trainer._input_labels
output_labels = trainer._output_labels
xmin, xmax = [0, 0.5], [1, 100]
input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}

pysmo_surr = surrogate.PysmoSurrogate(
    pysmo_surr_expr, input_labels, output_labels, input_bounds
)

surrogate_scatter2D(
    pysmo_surr, throughput_param_df, filename="scatter2D.pdf", show=True
)
surrogate_scatter2D(
    pysmo_surr, throughput_param_df, filename="scatter3D.pdf", show=True
)
surrogate_parity(pysmo_surr, throughput_param_df, filename="parity.pdf", show=True)
surrogate_residual(pysmo_surr, throughput_param_df, filename="residual.pdf", show=True)
