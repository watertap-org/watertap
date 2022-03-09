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

import sys
import os
import time

from watertap.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
from watertap.examples.flowsheets.lsrro.analysis import lsrro_paper_analysis as lsrro_case


def run_case(number_of_stages, Cin, water_recovery,  A_fixed, permeate_quality_limit, has_CP, nx):
    sweep_params = {}
    outputs = {}

    m = lsrro_case.build(number_of_stages, has_CP=has_CP)
    lsrro_case.set_operating_conditions(m, Cin=Cin)
    lsrro_case.initialize(m)
    lsrro_case.solve(m)
    lsrro_case.optimize_set_up(m, water_recovery=water_recovery, A_fixed=A_fixed, permeate_quality_limit=permeate_quality_limit)
    m, res = lsrro_case.solve(m, raise_on_failure=False, tee=False)

    sweep_params['Recovery Rate'] = LinearSample(m.fs.water_recovery, 0.3, 0.9, nx)
    # sweep_params['Feed Concentration'] = LinearSample(m.fs.feed.) #TODO: link to feed concentration, not mass flow NaCl
    outputs['LCOW'] = m.fs.costing.LCOW


    return m, res


if __name__ == "__main__":

    m, res = run_case(number_of_stages=3,
             Cin=70,
             water_recovery=0.5,
             A_fixed= 5/3.6e11,
             permeate_quality_limit=1000e-6,
             has_CP=True,
             nx=4
            )