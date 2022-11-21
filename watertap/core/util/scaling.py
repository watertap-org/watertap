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

from pyomo.environ import Expression, Param
import idaes.core.util.scaling as iscale


def transform_property_constraints(self):
    for metadata_dic in self.params.get_metadata().properties.values():
        var_str = metadata_dic["name"]
        if metadata_dic["method"] is not None and self.is_property_constructed(var_str):
            var = getattr(self, var_str)
            if isinstance(var, Expression):
                continue  # properties that are expressions do not have constraints
            if isinstance(var, Param):
                continue  # properties that are parameters do not have constraints
            con = getattr(self, "eq_" + var_str)
            for ind, c in con.items():
                sf = iscale.get_scaling_factor(var[ind], default=1, warning=True)
                iscale.constraint_scaling_transform(c, sf)
