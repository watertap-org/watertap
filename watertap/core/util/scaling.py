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
This module contains utility functions for testing the scaling of WaterTAP models.
"""

import pyomo.environ as pyo
from idaes.core import PhysicalParameterBlock

import idaes.logger
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    badly_scaled_var_generator,
    get_scaling_factor,
    unset_scaling_factor,
)

__author__ = "Hunter Barber"


def variable_sens_generator(blk, lb_scale=1e-2, ub_scale=1e2, tol=1e3, zero=1e-10):
    """
    Generator which varies the magnitude of inlet flow(s) to some lower and upper specification
        for a square and solved flowsheet. Previous scaling factors in the model are cleared and default
        scaling factors are adjusted to account for the new flow magnitude and again propagated through
        the variables on initialization. The model is re-solved and tested for: (i) are any variables
        poorly scaled at different process scales, (ii) are any scaled variables with static scaling
        factors proportional the process flow and therefore candidate for implementing dynamic scaling
        factors, (iii) do any variables touch their bounds between process scales that were not
        explicitly tested in other test frames, leading to non-optimal solves.

    Keyword Arguments:
        blk : pyomo ConcreteModel block which holds the flowsheet to be tested
        lb_scale : the lower bound scale which adjusts the magnitude for the
            fixed inlet flows of the process (default = 1e-2)
        ub_scale : the upper bound scale which adjusts the magnitude for the
            fixed inlet flows of the process (default = 1e2)
        tol : tolerance for the change in magnitude of scaled variables between
            the lower and upper bounds (default = 1e3), should always be smaller
            than the difference in magnitude of lb_scale and up_scale
        zero: magnitude that is considered to be zero, scaled variables with an absolute value
            of less than this value are okay, and not reported (default = 1e-10)

    Yields:
        string denoting what step failed, variable name as string,
            sensitivity of sv between lower and upper solutions
    """

    solver = get_solver()
    sv_hist = {}

    # get property block
    # TODO: Let multiple property packages be handled correctly
    for fs_blk in blk.block_data_objects(
        active=None, sort=False, descend_into=True, descent_order=None
    ):
        if isinstance(fs_blk, PhysicalParameterBlock):
            property_blk = fs_blk

    if property_blk is None:
        raise TypeError("No PhysicalParameterBlock exists in the Flowsheet")

    for scale in [lb_scale, ub_scale]:

        # clone flowsheet to re-solve at conditions supplied by scale
        temp_blk = blk.clone()

        # loop through all vars to wipe previous sf and reset dsf considering new scale
        for (var_name, var_index), var_obj in temp_blk.component_data_iterindex(
            ctype=pyo.Var, active=True, descend_into=True
        ):

            # unset scaling factor of each variable
            unset_scaling_factor(var_obj)
            # catch indexed variables that obtain sf propagated from an un-indexed parent
            if var_index is not None:
                unset_scaling_factor(var_obj.parent_component())

            # if the flow term has a default scaling factor, overwrite the dsf to new scale magnitude
            if (
                property_blk.get_default_scaling(var_name, index=var_index) is not None
                and "flow" in var_name
            ):

                # overwrite old dsf with one considering the new scale
                dsf = blk.fs.properties.get_default_scaling(var_name, index=var_index)
                temp_blk.fs.properties.set_default_scaling(
                    var_name, dsf / scale, index=var_index
                )

                if var_obj.fixed:
                    # overwrite the fixed variable value
                    var_obj.fix(var_obj * scale)

                # TODO: Determine how to include specified sf, e.g. ro.area

        # resolve the model at the new scale
        calculate_scaling_factors(temp_blk)
        temp_blk.fs.unit.initialize(outlvl=idaes.logger.ERROR)
        results = solver.solve(temp_blk)
        # ensure model solves wrt new scale
        if not pyo.check_optimal_termination(results):
            yield "Failed run on", scale, "scale"
            break

        # check variable scaling wrt new scale
        for var, sv in badly_scaled_var_generator(temp_blk, zero=zero):
            yield f"badly scaled variable for scaled flow of {scale} ", var.name, sv

        # store results for sv to sv_hist dict for current model copy
        for v in temp_blk.component_data_objects(
            pyo.Var, active=True, descend_into=True
        ):

            if not v.stale:
                continue

            sf = get_scaling_factor(v, default=1)
            val = pyo.value(v, exception=False)
            sv = abs(val * sf)
            if sv < zero:
                continue

            # store results in sv_hist dict
            if v.name not in sv_hist.keys():
                sv_hist[v.name] = [sv]
            else:
                sv_hist[v.name].append(sv)

        # clear temp_blk before running other scale
        temp_blk.clear()

    # loop through svs store in sv_hist to determine change in scaled variable
    for var_key, var_sv in sv_hist.items():

        # get order of magnitude difference (sensitivity)
        sens = var_sv[0] / var_sv[1]

        # check whether sensitivity is within tolerance
        if sens > tol or sens < tol**-1:
            # positive or negative can tell whether it is proportional or inversely proportional,
            # but will still be caught by tol value
            yield "high sensitivity of scaled variable", var_key, sens
