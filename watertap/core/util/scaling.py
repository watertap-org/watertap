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
from pyomo.environ import Block
from pyomo.environ import ConcreteModel

import idaes.logger
from idaes.core import MaterialFlowBasis
from idaes.core.util.model_statistics import fixed_variables_generator
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


# TODO: handle user specified variables as arg (ex/ m.fs.ro.area)
def variable_sens_generator(
    blk, lb_scale=1e-1, ub_scale=1e3, tol=1e2, zero=1e-10
):
    """

    Keyword Arguments:
        inlet_port : port objects to be scaled by model

    Returns:

    """
    '''
    var_hist = {}
    for scale in [lb_scale, ub_scale]:

        # clone 'm' model with all existing specifications
        temp_blk = blk.clone()

        # unset all scaling factors in the cloned model
        for var in temp_blk.component_data_objects(
                ctype=pyo.Var, active=True, sort=False, descend_into=True, descent_order=None
        ):
            unset_scaling_factor(var)

        # unset all scaling factors in the cloned model
        for tup_key, dsf in temp_blk.fs.properties.default_scaling_factor.items():
            if "flow" in tup_key[0]:
                # get default scaling factor from original model
                print(temp_blk.fs.properties.default_scaling_factor)
                dsf = blk.fs.properties.get_default_scaling(tup_key[0], index=tup_key[1])
                print(tup_key, dsf)
                # scale default scaling factor in temp model for resolve at new scale
                temp_blk.fs.properties.set_default_scaling(tup_key[0],  dsf * scale ** -1, index=tup_key[1])

        print(temp_blk.fs.properties.default_scaling_factor)
    '''
    '''
    for var1 in blk.component_data_objects(ctype=None, active=None, sort=False, descend_into=True, descent_order=None):
        print(var1, var1.ctype, type(var1), var1.__class__.__name__)
    '''
    '''
    for varP in blk.component_data_objects(ctype=Port, active=True, descend_into=True):
        print(varP)
        varPdict = varP.vars
        print(varPdict)
        for varP2 in varPdict.keys():
            print(varP2, varPdict[varP2])
    '''

    # TODO: handle multiple inlet ports

    var_hist = {}

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

        # store results for scale to var_hist for current model copy
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
            # store results in var_hist dict
            if v.name not in var_hist.keys():
                var_hist[v.name] = [sv]
            else:
                var_hist[v.name].append(sv)


    # loop through  store in var_hist to determine change in scaled variable
    for var_key, var_svs in var_hist.items():
        sens = var_svs[0] / var_svs[1]
        if sens > tol or sens < tol**-1:
            yield "high sensitivity of scaled variable", var_key, sens

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
