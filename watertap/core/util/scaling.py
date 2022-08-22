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
"""
This module contains utility functions for model statistics of WaterTAP models.
"""

import pyomo.environ as pyo
from pyomo.environ import check_optimal_termination
from pyomo.common.collections import ComponentSet

import idaes.logger
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    badly_scaled_var_generator,
    get_scaling_factor,
    unset_scaling_factor,
)

solver = get_solver()

__author__ = "Hunter Barber"


# TODO: handle user specified variables as arg (ex/ m.fs.ro.area)
def variable_sens_generator(
    blk, lb_scale=1e-1, ub_scale=1e3, tol=1e2, zero=1e-10
):
    """

    Keyword Arguments:

    Returns:
    """

    sv_hist = {}

    for scale in [lb_scale, ub_scale]:

        temp_blk = blk.clone()  # clone flowsheet to resolve at conditions supplied in scale

        # loop through all vars to wipe previous sf and reset dsf considering new scale
        for (var_name, var_index), var_obj in temp_blk.component_data_iterindex(
            ctype=pyo.Var, active=True, descend_into=True
        ):

            unset_scaling_factor(var_obj)  # remove prior sf which are reestablished on init

            # if the flow term has a default scaling factor and is fixed, as standard by flowsheets
            if (
                    blk.fs.properties.get_default_scaling(var_name, index=var_index) is not None
                    and "flow" in var_name
                    and var_obj.fixed
            ):

                # overwrite old dsf with one considering the new scale
                dsf = blk.fs.properties.get_default_scaling(var_name, index=var_index)
                temp_blk.fs.properties.set_default_scaling(var_name, dsf * scale ** -1, index=var_index)

                # overwrite the fixed variable value
                var_obj.fix(var_obj * scale)

        # resolve the model at the new scale
        calculate_scaling_factors(temp_blk)
        temp_blk.fs.unit.initialize(outlvl=idaes.logger.ERROR)
        results = solver.solve(temp_blk)

        # ensure model solves wrt new scale
        if not check_optimal_termination(results):
            yield "Failed run on", scale, "scale"
            break

        # check variable scaling wrt process flow scale
        for b in list(
            badly_scaled_var_generator(temp_blk, large=tol, small=tol ** -1)
        ):
            yield f"badly scaled variable for scaled flow of {scale} ", b[0].name, b[1]

        # store results for scale to sv_hist for current model copy
        for v in ComponentSet(
            temp_blk.component_data_objects(pyo.Var, active=True, descend_into=True)
        ):

            val = pyo.value(v, exception=False)
            if val is None:
                continue

            sf = get_scaling_factor(v, default=1)
            sv = abs(val * sf)
            if sv < zero:
                continue

            # store results in sv_hist dict
            if v.name not in sv_hist.keys():
                sv_hist[v.name] = [sv]
            else:
                sv_hist[v.name].append(sv)


    # loop through svs store in sv_hist to determine change in scaled variable
    for var_key, var_sv in sv_hist.items():

        # get order of magnitude difference
        sens = var_sv[0] / var_sv[1]

        # check whether sensitivity is within tolerance
        if sens > tol or sens < tol**-1:
            yield "high sensitivity of scaled variable", var_key, sens

