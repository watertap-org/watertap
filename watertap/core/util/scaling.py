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
        zero: magnitude that is considered to be zero, scaled variables with a value of
            zero are okay, and not reported (default = 1e-10)

    Yields:
        string denoting what step failed, variable name as string,
            sensitivity of sv between lower and upper solutions
    """

    # TODO: Check DOF is 0

    sv_hist = {}

    for scale in [lb_scale, ub_scale]:

        temp_blk = (
            blk.clone()
        )  # clone flowsheet to re-solve at conditions supplied by scale

        # need to handle objects with no data but that have sf that are propagated to data objects
        for obj in temp_blk.component_objects(
            ctype=pyo.Var,
            active=None,
            sort=False,
            descend_into=True,
            descent_order=None,
        ):
            unset_scaling_factor(obj)
            if obj.is_indexed:
                for (index, indexed_obj) in obj.items():
                    unset_scaling_factor(obj[index])

        # loop through all vars to wipe previous sf and reset dsf considering new scale
        for (var_name, var_index), var_obj in temp_blk.component_data_iterindex(
            ctype=pyo.Var, active=True, descend_into=True
        ):

            # repetitive with above loop
            unset_scaling_factor(var_obj)

            # if the flow term has a default scaling factor and is fixed, as standard by flowsheets,
            # overwrite the dsf to new scale magnitude
            # TODO: Check multiple inlet ports are handled correctly
            if (
                blk.fs.properties.get_default_scaling(var_name, index=var_index)
                is not None
                and "flow" in var_name
                and var_obj.fixed
            ):

                print(var_obj)

                # overwrite old dsf with one considering the new scale
                dsf = blk.fs.properties.get_default_scaling(var_name, index=var_index)
                temp_blk.fs.properties.set_default_scaling(
                    var_name, dsf * scale**-1, index=var_index
                )

                # overwrite the fixed variable value
                var_obj.fix(var_obj * scale)

                # TODO: Determine how to consider user specified sf, e.g. ro.area

        # resolve the model at the new scale
        calculate_scaling_factors(temp_blk)
        temp_blk.fs.unit.initialize(outlvl=idaes.logger.ERROR)
        results = solver.solve(temp_blk)
        print(results)
        # ensure model solves wrt new scale
        if not pyo.check_optimal_termination(results):
            yield "Failed run on", scale, "scale"
            break

        # check variable scaling wrt new scale
        for b in list(badly_scaled_var_generator(temp_blk, zero=zero)):
            yield f"badly scaled variable for scaled flow of {scale} ", b[0].name, b[1]

        # store results for sv to sv_hist dict for current model copy
        # TODO: Can this be moved to eliminate the need for the sv_hist dict
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

        # clear temp_blk before running other scale
        temp_blk.clear()

    # loop through svs store in sv_hist to determine change in scaled variable
    # TODO: Can the objects of the original blk be output, as opposed to variable
    #   name strings taken from the temp_blk
    for var_key, var_sv in sv_hist.items():

        # get order of magnitude difference (sensitivity)
        sens = var_sv[0] / var_sv[1]

        # check whether sensitivity is within tolerance
        if sens > tol or sens < tol**-1:
            # positive or negative can tell whether it is proportional or inversely proportional,
            # but will still be caught by tol value
            yield "high sensitivity of scaled variable", var_key, sens


'''
def _component_data_iteritems(self, ctype, active, sort, dedup):
    """return the name, index, and component data for matching ctypes

    Generator that returns a nested 2-tuple of

        ((component name, index value), _ComponentData)

    for every component data in the block matching the specified
    ctype(s).

    Parameters
    ----------
    ctype:  None or type or iterable
        Specifies the component types (`ctypes`) to include

    active: None or bool
        Filter components by the active flag

    sort: None or bool or SortComponents
        Iterate over the components in a specified sorted order

    dedup: _DeduplicateInfo
        Deduplicator to prevent returning the same _ComponentData twice
    """
    _sort_indices = SortComponents.sort_indices(sort)
    _subcomp = PseudoMap(self, ctype, active, sort)
    for name, comp in _subcomp.items():
        # NOTE: Suffix has a dict interface (something other derived
        #   non-indexed Components may do as well), so we don't want
        #   to test the existence of iteritems as a check for
        #   component datas. We will rely on is_indexed() to catch
        #   all the indexed components.  Then we will do special
        #   processing for the scalar components to catch the case
        #   where there are "sparse scalar components"
        if comp.is_indexed():
            _items = comp.items()
            if _sort_indices:
                _items = sorted_robust(_items, key=itemgetter(0))
        elif hasattr(comp, '_data'):
            # This is a Scalar component, which may be empty (e.g.,
            # from Constraint.Skip on a scalar Constraint).  Only
            # return a ComponentData if one officially exists.
            # Sorting is not a concern as this component has either
            # 0 or 1 datas
            assert len(comp._data) <= 1
            _items = comp._data.items()
        else:
            # This is a non-IndexedComponent Component.  Return it.
            _items = ((None, comp),)

        if active is None or not isinstance(comp, ActiveIndexedComponent):
            _items = (((name, idx), compData) for idx, compData in _items)
        else:
            _items = (((name, idx), compData) for idx, compData in _items
                      if compData.active == active)

        yield from dedup.unique(comp, _items, False)


def component_data_iterindex(self,
                             ctype=None,
                             active=None,
                             sort=False,
                             descend_into=True,
                             descent_order=None):
    """
    Return a generator that returns a tuple for each
    component data object in a block.  By default, this
    generator recursively descends into sub-blocks.  The
    tuple is

        ((component name, index value), _ComponentData)

    """
    dedup = _DeduplicateInfo()
    for _block in self.block_data_objects(
            active, sort, descend_into, descent_order):
        yield from _block._component_data_iteritems(
            ctype, active, sort, dedup)
'''
