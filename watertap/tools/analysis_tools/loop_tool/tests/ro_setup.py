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

import watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery as ro_erd
from watertap.tools.parameter_sweep.parameter_sweep import (
    ParameterSweep,
    RecursiveParameterSweep,
)
from watertap.tools.parameter_sweep.parameter_sweep_reader import ParameterSweepReader

from watertap.tools.parameter_sweep.parameter_sweep_differential import (
    DifferentialParameterSweep,
)
from idaes.core.solvers import get_solver
import os

__author__ = "Alexander V. Dudchenko (SLAC)"


def ro_build(**kwargs):
    m = ro_erd.build(**kwargs)
    return m


def ro_init(m, solver=None, **kwargs):
    ro_erd.set_operating_conditions(
        m, water_recovery=0.5, over_pressure=0.3, solver=solver
    )
    ro_erd.initialize_system(m, solver=solver)
    ro_erd.optimize_set_up(m)


def ro_solve(m, solver=None, **kwargs):
    result = solver.solve(m)
    return result
