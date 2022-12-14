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

import pyomo.environ as pyo
import idaes.core.util.scaling as iscale


def transform_property_constraints(self):
    for metadata_dic in self.params.get_metadata().properties.values():
        var_str = metadata_dic["name"]
        if metadata_dic["method"] is not None and self.is_property_constructed(var_str):
            var = getattr(self, var_str)
            if not isinstance(var, pyo.Var):
                continue  # properties that are not vars do not have constraints
            con = getattr(self, "eq_" + var_str)
            for ind, c in con.items():
                sf = iscale.get_scaling_factor(var[ind], default=1, warning=True)
                iscale.constraint_scaling_transform(c, sf)
