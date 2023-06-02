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
This module contains a utility functions for scaling WaterTAP models.
"""
import numpy as np
import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.common.modeling import unique_component_name
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def set_autoscaling_factors(
    m,
    ignore_constraint_scaling=False,
    ignore_variable_scaling=False,
    max_grad=100,
    min_scale=1e-6,
    equilibriate_variables=True,
):
    # Pynumero requires an objective, but I don't, so let's see if we have one
    n_obj = 0
    for c in m.component_data_objects(pyo.Objective, active=True):
        n_obj += 1
    # Add an objective if there isn't one
    if n_obj == 0:
        dummy_objective_name = unique_component_name(m, "objective")
        setattr(m, dummy_objective_name, pyo.Objective(expr=0))
    # Create NLP and calculate the objective
    if not AmplInterface.available():
        raise RuntimeError("Pynumero not available.")
    nlp = PyomoNLP(m)
    jac = nlp.evaluate_jacobian().tocsc()
    if not ignore_variable_scaling:
        for i, v in enumerate(nlp.get_pyomo_variables()):
            current_factor = iscale.get_scaling_factor(v, default=1.9)
            if equilibriate_variables:
                abs_data = np.abs(jac.getcol(i).data)
                if len(abs_data) == 0:
                    continue
                new_factor = abs_data.max()
                iscale.set_scaling_factor(v, new_factor, data_objects=False)
            else:
                new_factor = current_factor
            # if current_factor != new_factor:
            #    print(f"updated scaling factor for {v.name} from {current_factor} to {new_factor}")
            jac[:, i] *= 1.0 / new_factor

    jac = jac.tocsr()
    # calculate constraint scale factors
    for i, c in enumerate(nlp.get_pyomo_constraints()):
        sc = iscale.get_scaling_factor(c, default=1)
        if ignore_constraint_scaling or iscale.get_scaling_factor(c) is None:
            sc = 1
            abs_row = np.abs(jac.getrow(i).data)
            mg = abs_row.max()
            if mg > max_grad:
                sc = max(min_scale, max_grad / mg)
            iscale.set_scaling_factor(c, sc, data_objects=False)

    # delete dummy objective
    if n_obj == 0:
        delattr(m, dummy_objective_name)
    return nlp


def transform_property_constraints(self):
    """
    This is a utility function for the scaling of WaterTAP property model constraints.
    """
    for p in self.params.get_metadata().properties.list_supported_properties():
        var_str = p.name
        if p.method is not None and self.is_property_constructed(var_str):
            var = getattr(self, var_str)
            if not isinstance(var, pyo.Var):
                continue  # properties that are not vars do not have constraints
            # adding a conditional to check if a constraint exists for property; in the case when we only add and object reference, there would not be a constraint
            if hasattr(self, "eq_" + var_str):
                con = getattr(self, "eq_" + var_str)
                for ind, c in con.items():
                    sf = iscale.get_scaling_factor(var[ind], default=1, warning=True)
                    iscale.constraint_scaling_transform(c, sf)
            else:
                msg = (
                    f"If there was a property constraint written for the variable, {var}, that constraint was not "
                    f"scaled. The transform_property_constraints tool expects constraints to have the following naming "
                    f"convention: 'eq_' + '{var_str}'. This suggests that the user may have defined a property in "
                    f"metadata but failed to follow the naming convention for its constraint. If there is no property "
                    f"constraint associated with the {var_str}, this warning can be ignored."
                )
                _log.warning(msg)
