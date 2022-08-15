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


def variable_sens_generator(blk, lb_scale=1e-1, ub_scale=1e3, tol=1e2, zero=1e-10):
    """

    Keyword Arguments:

    Returns:

    """

    test_scale = [lb_scale, ub_scale]
    var_hist = {}
    for scale in test_scale:
        temp_blk = blk.clone()
        # loop through variables and scale previously fixed inlet flow
        for vt, vm in temp_blk.component_data_iterindex(
            ctype=pyo.Var, active=True, descend_into=True
        ):
            unset_scaling_factor(vm)  # remove prior sf which are reestablished on init
            if (
                temp_blk.fs.properties.get_default_scaling(vt[0], index=vt[1])
                is not None
                and "flow" in vt[0]
                and vm.fixed
            ):
                dsf = temp_blk.fs.properties.get_default_scaling(vt[0], index=vt[1])
                temp_blk.fs.properties.set_default_scaling(
                    vt[0], dsf * scale**-1, index=vt[1]
                )
                vm.fix(vm * scale)
        calculate_scaling_factors(temp_blk)
        temp_blk.fs.unit.initialize(outlvl=idaes.logger.ERROR)
        results = solver.solve(temp_blk)
        if not check_optimal_termination(results):
            print("Failed run on", scale, "scale")

        # store results for scale to var_hist
        for v in ComponentSet(
            temp_blk.component_data_objects(pyo.Var, active=True, descend_into=True)
        ):
            val = pyo.value(v, exception=False)
            if val is None:
                continue
            sf = get_scaling_factor(v, default=1)
            sv = abs(val * sf)  # scaled value
            if sv < zero:
                continue
            if v.name not in var_hist.keys():
                var_hist[v.name] = [sv]
            else:
                var_hist[v.name].append(sv)

    badly_scaled_vars = list(
        badly_scaled_var_generator(temp_blk, large=1e2, small=1e-2)
    )
    for b in badly_scaled_vars:
        yield "badly scaled variable", b[0].name, b[1]

    # loop through variables to select minmax
    for k in var_hist.keys():
        sens = max(var_hist[k]) / min(var_hist[k])
        if sens > tol or sens < tol**-1:
            yield "high sensitivity of scaled variable", k, sens
