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
Tests for ADM1 flowsheet example.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

"""

# Some more information about this module
__author__ = "Alejandro Garciadiego"

import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent
import numpy as np
import pandas as pd
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from pyomo.opt import TerminationCondition, SolverStatus

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.BSM2.ASM1_ADM1_flowsheet import (
    build_flowsheet,
)


def model():
    m, res = build_flowsheet()

    N = 2  # number of samples
    x_s = 0.2
    st_dev_x = 0.04

    m.results = res

    solver = get_solver(options={"bound_push": 1e-8})

    Feed = []
    S_I = []
    S_S = []
    X_I = []
    X_S = []
    X_BH = []
    X_BA = []
    X_P = []

    Feedt = []
    S_It = []
    S_St = []
    X_It = []
    X_St = []
    X_BHt = []
    X_BAt = []
    X_Pt = []

    Feedg = []
    ch4 = []
    co2 = []
    flow = []

    for i in range(N):
        noise_x = np.random.normal(loc=0, scale=st_dev_x)
        print(value(m.fs.FeedWater.conc_mass_comp[0, "X_S"]))
        m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(x_s + noise_x)
        print(i)
        print(value(m.fs.FeedWater.conc_mass_comp[0, "X_S"]))

        try:

            # solve the model
            status = solver.solve(m, tee=False)

            if (status.solver.status == SolverStatus.ok) and (
                status.solver.termination_condition == TerminationCondition.optimal
            ):
                print("yes")

                Feed.append(value(m.fs.FeedWater.conc_mass_comp[0, "X_S"]))
                S_I.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "S_I"]))
                S_S.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "S_S"]))
                X_I.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "X_I"]))
                X_S.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "X_S"]))
                X_BH.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "X_BH"]))
                X_BA.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "X_BA"]))
                X_P.append(value(m.fs.adm_asm.outlet.conc_mass_comp[0, "X_P"]))

                Feedt.append(value(m.fs.FeedWater.conc_mass_comp[0, "X_S"]))
                S_It.append(value(m.fs.Treated.conc_mass_comp[0, "S_I"]))
                S_St.append(value(m.fs.Treated.conc_mass_comp[0, "S_S"]))
                X_It.append(value(m.fs.Treated.conc_mass_comp[0, "X_I"]))
                X_St.append(value(m.fs.Treated.conc_mass_comp[0, "X_S"]))
                X_BHt.append(value(m.fs.Treated.conc_mass_comp[0, "X_BH"]))
                X_BAt.append(value(m.fs.Treated.conc_mass_comp[0, "X_BA"]))
                X_Pt.append(value(m.fs.Treated.conc_mass_comp[0, "X_P"]))

                Feedg.append(value(m.fs.FeedWater.conc_mass_comp[0, "X_S"]))
                ch4.append(value(m.fs.RADM.vapor_outlet.conc_mass_comp[0, "S_ch4"]))
                co2.append(value(m.fs.RADM.vapor_outlet.conc_mass_comp[0, "S_co2"]))
                flow.append(value(m.fs.RADM.vapor_outlet.flow_vol[0]))

        except ValueError:
            Feed[i] = "NaN"

    df = pd.DataFrame(
        {
            "Feed": Feed,
            "S_I": S_I,
            "S_S": S_S,
            "X_I": X_I,
            "X_S": X_S,
            "X_BH": X_BH,
            "X_BA": X_BA,
            "X_P": X_P,
        }
    )
    df.to_csv("results.csv")
    dft = pd.DataFrame(
        {
            "Feedt": Feed,
            "S_It": S_It,
            "S_St": S_St,
            "X_It": X_It,
            "X_St": X_St,
            "X_BHt": X_BHt,
            "X_BA": X_BAt,
            "X_P": X_Pt,
        }
    )
    dft.to_csv("resultst.csv")
    dfg = pd.DataFrame({"Feedg": Feed, "ch4": ch4, "co2": co2, "flow": flow})
    dfg.to_csv("resultstg.csv")

    return m


x = model()
