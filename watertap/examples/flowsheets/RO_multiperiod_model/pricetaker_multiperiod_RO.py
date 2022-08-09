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
    Param,
    ConcreteModel,
    Var,
    units as pyunits,
)

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as swro
from watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO import (
    create_multiperiod_swro_model,
)


def _get_lmp(time_steps, data_path):
    """
    Get price signals from data set
    :param time_steps: Number of time steps considered in MP analysis
    :param data_path: Price [$/kWh] on the same interval as time_steps
    :return: reshaped data
    """
    # read data into array from file
    lmp_data = np.genfromtxt(data_path, delimiter=",")

    # index only the desired number of timesteps
    return lmp_data[:time_steps]


def run_pricetaker_analysis(ndays=1, data_path="dagget_CA_LMP_hourly_2015.csv"):

    # number of time steps
    n_steps = ndays * 24

    # get data
    lmp = _get_lmp(n_steps, data_path)

    # create mp model
    mp_swro = create_multiperiod_swro_model(n_time_points=n_steps)

    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()

    for count, blk in enumerate(t_blocks):
        blk_swro = blk.ro_mp
        blk.lmp_signal = Param(default=0, mutable=True)

        # set the electricity_price in each flowsheet
        blk_swro.fs.costing.electricity_price = lmp[count]


if __name__ == "__main__":
    print(_get_lmp(24))
