import pyomo.environ as pyo
import copy
import numpy as np
import os
import idaes.logger as idaeslog
import platform
import time
from pyomo.common.collections import ComponentSet, ComponentMap


from idaes.core.util import to_json, from_json

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver


def do_async_sweep(self, model, sweep_params, outputs, local_values):
    try:
        import ray
        from watertap.tools.parameter_sweep.async_param_sweep import rayio_sweep

        return rayio_sweep.do_rayio_sweep(
            self, model, sweep_params, outputs, local_values
        )
    except ImportError:
        from watertap.tools.parameter_sweep.async_param_sweep import (
            multiprocessing_sweep,
        )

        return multiprocessing_sweep.do_mp_sweep(
            self, model, sweep_params, outputs, local_values
        )
