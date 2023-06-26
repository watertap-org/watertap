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

from pyomo.environ import assert_optimal_termination, value, units
from pyomo.util.check_units import assert_units_consistent
import numpy as np
import pandas as pd
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from pyomo.opt import TerminationCondition, SolverStatus

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.electroNP.electroNP_flowsheet import (
    build_flowsheet,
)


def model():
    m, res = build_flowsheet()

    N = 100  # number of samples
    x_s = 3.4655
    st_dev_x = 0.05

    m.results = res

    solver = get_solver(options={"bound_push": 1e-8})

    Feed = []
    S_A = []
    S_F = []
    S_I = []
    S_N2 = []
    S_NH4 = []
    S_NO3 = []
    S_O2 = []
    S_PO4 = []
    X_AUT = []
    X_H = []
    X_I = []
    X_PAO = []
    X_PHA = []
    X_PP = []
    X_S = []

    c_cost = []
    o_cost = []
    lcow = []
    e_cost = []

    for i in range(N):
        noise_x = np.random.normal(loc=0, scale=st_dev_x)
        # print(value(m.fs.FeedWater.conc_mass_comp[0, "S_su"]))
        m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"].fix(x_s + noise_x)
        print(i)
        print(value(m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"]))

        try:

            # solve the model
            status = solver.solve(m, tee=False)

            if (status.solver.status == SolverStatus.ok) and (
                status.solver.termination_condition == TerminationCondition.optimal
            ):
                print("yes")

                Feed.append(value(m.fs.AD.inlet.conc_mass_comp[0, "X_PAO"]))
                S_A.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_A"]))
                S_F.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_F"]))
                S_I.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_I"]))
                S_N2.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_N2"]))
                S_NH4.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_NH4"]))
                S_NO3.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_NO3"]))
                S_O2.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_O2"]))
                S_PO4.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "S_PO4"]))
                X_AUT.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_AUT"]))
                X_H.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_H"]))
                X_I.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_I"]))
                X_PAO.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_PAO"]))
                X_PHA.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_PHA"]))
                X_PP.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_PP"]))
                X_S.append(value(m.fs.electroNP.treated.conc_mass_comp[0, "X_S"]))

                c_cost.append(
                    value(
                        units.convert(
                            m.fs.costing.total_capital_cost, to_units=units.USD_2018
                        )
                    )
                )
                o_cost.append(
                    value(
                        units.convert(
                            m.fs.costing.total_operating_cost,
                            to_units=units.USD_2018 / units.year,
                        )
                    )
                )
                lcow.append(
                    value(
                        units.convert(
                            m.fs.costing.LCOW, to_units=units.USD_2018 / units.m**3
                        )
                    )
                )
                e_cost.append(
                    value(
                        units.convert(
                            m.fs.electroNP.energy_electric_flow_mass,
                            to_units=units.kWh / units.kg,
                        )
                    )
                )

        except ValueError:
            pass
            # Feed.append("NaN")
            # S_A.append("NaN")
            # S_F.append("NaN")
            # S_I.append("NaN")
            # S_N2.append("NaN")
            # S_NH4.append("NaN")
            # S_NO3.append("NaN")
            # S_O2.append("NaN")
            # S_PO4.append("NaN")
            # X_AUT.append("NaN")
            # X_H.append("NaN")
            # X_I.append("NaN")
            # X_PAO.append("NaN")
            # X_PHA.append("NaN")
            # X_PP.append("NaN")
            # X_S.append("NaN")

            # c_cost.append("NaN")
            # o_cost.append("NaN")
            # lcow.append("NaN")
            # e_cost.append("NaN")

    df = pd.DataFrame(
        {
            "Feed": Feed,
            "S_A": S_A,
            "S_F": S_F,
            "S_I": S_I,
            "S_N2": S_N2,
            "S_NH4": S_NH4,
            "S_NO3": S_NO3,
            "S_O2": S_O2,
            "S_PO4": S_PO4,
            "X_AUT": X_AUT,
            "X_H": X_H,
            "X_I": X_I,
            "X_PAO": X_PAO,
            "X_PHA": X_PHA,
            "X_PP": X_PP,
            "X_S": X_S,
            "c_cost": c_cost,
            "o_cost": o_cost,
            "lcow": lcow,
            "e_cost": e_cost,
        }
    )
    df.to_csv("results.csv")

    return m


x = model()
