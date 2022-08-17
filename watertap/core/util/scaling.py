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
from pyomo.network import Port
from pyomo.environ import ConcreteModel

import idaes.logger
from idaes.core import MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    badly_scaled_var_generator,
    get_scaling_factor,
    unset_scaling_factor,
)
from idaes.core.util.model_statistics import fixed_variables_set

solver = get_solver()

__author__ = "Hunter Barber"


# TODO: handle user specified varibles as arg (ex/ m.fs.ro.area)
def variable_sens_generator(
    inlet_port, lb_scale=1e-1, ub_scale=1e3, tol=1e2, zero=1e-10
):
    """

    Keyword Arguments:
        inlet_port : port objects to be scaled by model

    Returns:

    """

    # TODO: handle multiple inlet ports
    fs = inlet_port.flowsheet()
    print(inlet_port)
    print(fs)
    print(fs.parent_block())

    var_hist = {}
    for scale in [lb_scale, ub_scale]:
        temp_fs = (
            fs.clone()
        )  # clone flowsheet for loop of lower bound scale and upper bound scale
        # loop through variables and scale previously fixed inlet flow
        print(temp_fs.inlet_port)
        for vt, vm in temp_fs.component_data_iterindex(
            ctype=pyo.Var, active=True, descend_into=True
        ):
            unset_scaling_factor(
                vm
            )  # remove prior sf which are reestablished on init, reapplied on
            if (
                temp_fs.properties.get_default_scaling(vt[0], index=vt[1]) is not None
                and "flow" in vt[0]
                and vm.fixed
            ):
                dsf = temp_fs.properties.get_default_scaling(vt[0], index=vt[1])
                temp_fs.properties.set_default_scaling(
                    vt[0], dsf * scale**-1, index=vt[1]
                )
                vm.fix(vm * scale)
        calculate_scaling_factors(
            temp_fs
        )  # recalculate scaling factors that were removed
        temp_fs.unit.initialize(
            outlvl=idaes.logger.ERROR
        )  # reinitialize model and supress logger
        results = solver.solve(temp_fs)  # resolve model copy
        if not check_optimal_termination(results):
            print(
                "Failed run on", scale, "scale"
            )  # make sure results are feasible for model copy

        # store results for scale to var_hist for current model copy
        for v in ComponentSet(
            temp_fs.component_data_objects(pyo.Var, active=True, descend_into=True)
        ):
            val = pyo.value(v, exception=False)
            if val is None:
                continue
            sf = get_scaling_factor(v, default=1)
            sv = abs(val * sf)  # scaled value
            if sv < zero:
                continue
            # store results in var_hist dict
            if v.name not in var_hist.keys():
                var_hist[v.name] = [sv]
            else:
                var_hist[v.name].append(sv)

        # check variable scaling wrt process flow scale
        badly_scaled_vars = list(
            badly_scaled_var_generator(temp_fs, large=1e2, small=1e-2)
        )
        for b in badly_scaled_vars:
            yield f"badly scaled variable for scaled flow of {scale} ", b[0].name, b[1]

    # loop through  store in var_hist to determine change in scaled variable
    for k in var_hist.keys():
        sens = max(var_hist[k]) / min(var_hist[k])
        if sens > tol or sens < tol**-1:
            yield "high sensitivity of scaled variable", k, sens

    """ 
    print(blk.fs.unit.process_flow.properties_in[0])
    print(blk.fs.unit.process_flow.properties_in[0].parent_block())
    print(blk.fs.unit.process_flow.properties_in[0].flowsheet())

    print(blk.fs.properties.build_state_block(
        [0], default={"defined_state": True}))

    blk.fs.state_block_ref = blk.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )
    state_vars_dict = blk.fs.state_block_ref[0].define_state_vars()
    for key in state_vars_dict.keys():
        if "flow" in key:
            print(key)

    for var in fixed_variables_set(blk):
        print(var)
        print(var.local_name)
        print(var.index())


    for port in blk.component_data_objects(Port, active=True, descend_into=True):
        print(port)
        print(port.doc)
        for var in port.vars:
            print(var)

    state = blk.fs.properties.build_state_block([0])
    print(blk.fs.unit.process_flow.properties_in[0].get_material_flow_basis())
    """
