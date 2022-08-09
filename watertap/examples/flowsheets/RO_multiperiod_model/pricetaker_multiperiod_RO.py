###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import os
import numpy as np

from pyomo.environ import (
    NonNegativeReals,
    ConcreteModel,
    Var,
    units as pyunits,
)

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
import watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO as mp_swro


def _get_lmp(time_steps, data_path="dagget_CA_LMP_hourly_2015.csv"):
    # read data into array from file
    lmp_data = np.genfromtxt(data_path, delimiter=",")

    # index only the desired number of timesteps
    return lmp_data[:time_steps]


if __name__ == "__main__":
    print(_get_lmp(24))
