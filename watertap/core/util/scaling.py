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
"""
This module contains a utility function for the scaling of WaterTAP property model constraints.
"""
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def transform_property_constraints(self):
    for p in self.params.get_metadata().properties.list_supported_properties():
        var_str = p.name
        if p.method is not None and self.is_property_constructed(var_str):
            var = getattr(self, var_str)
            # Some property models may not use `None` as the method for state variables and wouldn't have a constraint
            if var_str in list(self.define_state_vars()):
                continue
            msg = (
                f"If there was a property constraint written for the variable, {var}, that constraint was not "
                f"scaled. The transform_property_constraints tool expects constraints to have the following naming "
                f"convention: 'eq_' + '{var_str}'. This suggests that the user may have defined a property in "
                f"metadata but failed to follow the naming convention for its constraint. If there is no property "
                f"constraint associated with the {var_str}, this warning can be ignored."
            )
            if not isinstance(var, pyo.Var):
                continue  # properties that are not vars do not have constraints
            # adding a conditional to check if a constraint exists for property; in the case when we only add an object reference, there would not be a constraint
            if hasattr(self, "eq_" + var_str):
                con = getattr(self, "eq_" + var_str)
                for ind, c in con.items():
                    sf = iscale.get_scaling_factor(var[ind], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)
            # in some cases, a variable could be fixed within the parameter block, and the logger warning shouldn't apply.
            elif not var.is_indexed() and var.is_fixed():
                pass
            # in some cases, an indexed variable could be fixed within the parameter block, and the logger warning shouldn't apply
            # to those indexed variables that are fixed.
            elif var.is_indexed():
                for indexed_var in var.values():
                    if indexed_var.is_fixed():
                        pass
                    else:
                        _log.warning(msg)
            else:
                _log.warning(msg)
