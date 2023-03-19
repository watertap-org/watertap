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
###############################################################################

import pytest
import pyomo.environ as pyo
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc
from pyomo.opt import SolverFactory

import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    large_residuals_set,
    variables_near_bounds_generator,
    total_constraints_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_diagnostics import DegeneracyHunter
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
    get_scaling_factor,
    set_scaling_factor,
)
from idaes.core import UnitModelCostingBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    DiffusivityCalculation,
)
from idaes.models.unit_models import Feed, Product
from watertap.unit_models.gac_rssct import (
    GAC_RSSCT,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

import numpy as np
import pandas as pd
import pyomo.contrib.parmest.parmest as parmest
import matplotlib.pyplot as plt
from os.path import join, abspath, dirname

__author__ = "Hunter Barber"

# obtain experimental RSSCT data from csv
source_name_list = [
    "OCWD",
    "Serrano Water District",
    "Anaheim",
    "Fullerton F-3A/1",
    "Fullerton F-5/1",
    "Santa Ana",
    "Tustin",
    "Orange",
    "Garden Grove",
    "IRWD",
]
# source_name = "Santa Ana"  # "OCWD" "Orange" "Serrano Water District"
data = pd.read_csv("F400_PFOA.csv")

# initial guess for regressed variables
guess_freund_k = 10
guess_freund_ninv = 0.8
guess_ds = 1e-15

# set up solver
solver = get_solver()

log = idaeslog.getSolveLogger("solver.demo")
log.setLevel(idaeslog.DEBUG)

pyo.SolverFactory.register("ipopt")(pyo.SolverFactory.get_class("ipopt-watertap"))


def main():

    global source_name

    data_regression = {}
    data_filtered = {}
    param_results = {}

    for source_name in source_name_list:

        print(f"--------------------------\t{source_name}\t--------------------------")
        # data set
        data_set = data[["data_iter", f"{source_name}_X", f"{source_name}_Y"]]
        data_filtered[f"{source_name}"] = data_filter(data_set)
        data_filtered_case = data_filtered[f"{source_name}"]

        # vars to estimate
        theta_names = [
            "fs.gac.freund_k",
            "fs.gac.freund_ninv",
            "fs.gac.ds",
        ]
        # vars initial values
        theta_values = pd.DataFrame(
            data=[[guess_freund_k, guess_freund_ninv, guess_ds]],
            columns=theta_names,
        )
    
        # sum of squared error function as objective
        expr_sf = 1e-10
    
        def SSE(model, data_filtered_case):
            expr = (float(data_filtered_case[f"{source_name}_X"]) - model.fs.gac.bed_volumes_treated) ** 2
            return expr * expr_sf
    
        print("--------------------------\tmodel regression\t--------------------------")
        # Create an instance of the parmest estimator
        pest = parmest.Estimator(
            parmest_regression,
            data_filtered_case,
            theta_names,
            SSE,
            tee=False,
        )
    
        pest.objective_at_theta(
            theta_values=theta_values,
            initialize_parmest_model=True,
        )
    
        # Parameter estimation
        obj, theta = pest.theta_est()
    
        print("The SSE at the optimal solution is %0.6f" % (obj / expr_sf))
        print("The values for the parameters are as follows:")
        for k, v in theta.items():
            print(k, "=", v)

        param_results[f"{source_name}"] = theta

        print("--------------------------\tmodel rebuild\t--------------------------")
        # rebuild model across CP to view regression results
        # theta = [guess_freund_k, guess_freund_ninv, guess_ds]
        df_iter = solve_regression(theta)

        data_regression[f"{source_name}"] = df_iter

    print("--------------------------\tshow plot\t--------------------------")
    plot_regression(data_regression, data_filtered)
    print(param_results)



def model_build():

    # build models
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # PFOA (C8HF15O2) molar volume by LeBas method
    mv_pfas = ((8 * 14.8) + (1 * 3.7) + (15 * 8.7) + (7.4 * 1 + 12 * 1)) * 1e-6
    m.fs.properties = MCASParameterBlock(
        solute_list=["PFOA"],
        mw_data={"H2O": 0.018, "PFOA": 0.41407},
        diffus_calculation=DiffusivityCalculation.HaydukLaudie,
        molar_volume_data={("Liq", "PFOA"): mv_pfas},
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 1000
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.gac = GAC_RSSCT(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="calculated",
        surface_diffusion_coefficient_type="fixed",
    )

    # touch properties and scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e15, index=("Liq", "PFOA")
    )
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 8.5e-8,  # 5.1 mL/min (OCWD, 2021)
            (
                "conc_mass_phase_comp",
                ("Liq", "PFOA"),
            ): 6.59e-9,  # 6.59 ng/L (OCWD, 2021)
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # streams
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.gac.b4 = 1

    # adsorption isotherm
    m.fs.gac.freund_k.fix(guess_freund_k)
    m.fs.gac.freund_ninv.fix(guess_freund_ninv)
    m.fs.gac.ds.fix(guess_ds)
    # gac particle specifications
    m.fs.gac.particle_dens_app.fix(962.89)
    m.fs.gac.particle_dia.fix(0.000127)
    # adsorber bed specifications
    m.fs.gac.ebct.fix(13.68)
    m.fs.gac.bed_voidage.fix(0.44)
    m.fs.gac.velocity_sup.fix(0.002209)
    # design spec
    m.fs.gac.conc_ratio_replace.fix(0.20)
    # parameters
    m.fs.gac.a0.fix(0.8)
    m.fs.gac.a1.fix(0)
    m.fs.gac.b0.fix(0.023)
    m.fs.gac.b1.fix(0.793673)
    m.fs.gac.b2.fix(0.039324)
    m.fs.gac.b3.fix(0.009326)
    m.fs.gac.b4.fix(0.08275)
    m.fs.gac.shape_correction_factor.fix()
    m.fs.gac.conc_ratio_start_breakthrough = 0.001

    return m


def parmest_regression(data):

    print(f'running regression case {int(data["data_iter"])}')
    m = model_build()

    # scaling
    model_scale(m)

    # initialization
    model_init(m, outlvl=idaeslog.ERROR)

    # re-fix to conc_ratio spec from data
    m.fs.gac.conc_ratio_replace.fix(float(data[f"{source_name}_Y"]))

    return m


def solve_regression(theta):

    # build model
    m = model_build()
    # scaling
    model_scale(m)
    # initialization
    model_init(m, outlvl=idaeslog.ERROR)

    # refix differing parameters
    m.fs.gac.freund_k.fix(theta[0])
    m.fs.gac.freund_ninv.fix(theta[1])
    m.fs.gac.ds.fix(theta[2])

    # initialization
    model_init(m, outlvl=idaeslog.ERROR)

    # check profile against pilot data
    conc_ratio_input = np.linspace(0.01, 0.95, 40)
    conc_ratio_list = []
    bed_volumes_treated_list = []

    for conc_ratio in conc_ratio_input:

        print("running conc_ratio", conc_ratio)

        # fix to conc_ratio_case
        m.fs.gac.conc_ratio_replace.fix(conc_ratio)

        # resolve initialized model
        model_solve(m, solver_log=False)

        # record variables of interest
        conc_ratio_list.append(m.fs.gac.conc_ratio_replace.value)
        bed_volumes_treated_list.append(m.fs.gac.bed_volumes_treated.value)

    df_iter = {
        "conc_ratio": conc_ratio_list,
        "bed_volumes_treated": bed_volumes_treated_list,
    }

    return df_iter


def plot_regression(data_regression, data_filtered):

    color_code = [
        '#ffc000',
        '#ed7d31',
        '#c00000',
        '#00b050',
        '#7030a0',
        '#0070c0',
        '#bdd7ee',
        '#000000',
        '#a5a5a5',
        '#f8cbad',
    ]

    # fig1
    fig1, axs = plt.subplots(
        nrows=5,
        ncols=2,
        sharex=True,
        sharey=True,
    )

    i = 0
    for ax in axs.flat:

        ax.plot(
            data[f"{source_name_list[i]}_X"],
            data[f"{source_name_list[i]}_Y"],
            "o",
            mec=color_code[i],
            mfc='None',
            label="raw data",
        )
        ax.plot(
            data_filtered[source_name_list[i]][f"{source_name_list[i]}_X"],
            data_filtered[source_name_list[i]][f"{source_name_list[i]}_Y"],
            "o",
            mec=color_code[i],
            mfc=color_code[i],
            label=f"filtered data",
        )
        ax.plot(
            data_regression[source_name_list[i]]["bed_volumes_treated"],
            data_regression[source_name_list[i]]["conc_ratio"],
            color_code[i],
            label="regression results"
        )
        ax.legend(
            title=f"{source_name_list[i]}",
            fontsize='x-small',
            loc='lower right',
        )

        if i == 8:
            ax.set_xlabel("bed volumes treated")
            ax.set_ylabel("concentration ratio")

        i = i+1

    fig1.subplots_adjust(
        left=0.1, bottom=0.1, right=0.9, top=0.9,
        wspace=0, hspace=0
    )

    plt.show()

    # fig2
    plt.figure(2)
    for j in range(10):

        plt.plot(
            data_regression[source_name_list[j]]["bed_volumes_treated"],
            data_regression[source_name_list[j]]["conc_ratio"],
            color_code[j],
            label=f"{source_name_list[j]}"
        )
        plt.plot(
            data[f"{source_name_list[j]}_X"],
            data[f"{source_name_list[j]}_Y"],
            "o",
            mec=color_code[j],
            mfc='None',
        )
        plt.plot(
            data_filtered[source_name_list[j]][f"{source_name_list[j]}_X"],
            data_filtered[source_name_list[j]][f"{source_name_list[j]}_Y"],
            "o",
            mec=color_code[j],
            mfc=color_code[j],
        )

    plt.xlabel("bed volumes treated")
    plt.ylabel("concentration ratio")
    plt.legend(
        fontsize='x-small',
        loc='lower right',
    )

    plt.show()


def model_scale(model):

    set_scaling_factor(model.fs.gac.ds, 1e17)
    set_scaling_factor(model.fs.gac.dg, 1e-6)
    set_scaling_factor(model.fs.gac.bed_diameter, 1e3)
    set_scaling_factor(model.fs.gac.bed_mass_gac, 1e3)
    set_scaling_factor(model.fs.gac.mass_adsorbed, 1e9)
    set_scaling_factor(model.fs.gac.gac_usage_rate, 1e11)

    iscale.calculate_scaling_factors(model)


def model_init(model, outlvl):

    optarg = solver.options
    model.fs.feed.initialize(optarg=optarg, outlvl=outlvl)
    propagate_state(model.fs.s01)
    model.fs.gac.initialize(optarg=optarg, outlvl=outlvl)


def model_solve(model, solver_log=True):

    # check model
    assert_units_consistent(model)  # check that units are consistent
    assert degrees_of_freedom(model) == 0

    # solve simulation
    if solver_log:
        with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
            term_cond = results.solver.termination_condition
            print("termination condition:", term_cond)
    else:
        results = solver.solve(model, tee=False)
        term_cond = results.solver.termination_condition
        print("termination condition:", term_cond)

    if not term_cond == "optimal":
        badly_scaled_var_list = badly_scaled_var_generator(model)
        print("----------------\tbadly_scaled_var_list\t----------------")
        for x in badly_scaled_var_list:
            print(f"{x[0].name}\t{x[0].value}\tsf: {get_scaling_factor(x[0])}")
        print("----------------\tvariables_near_bounds_list\t----------------")
        variables_near_bounds_list = variables_near_bounds_generator(model)
        for x in variables_near_bounds_list:
            print(f"{x.name}\t{x.value}")
        print("----------------\ttotal_constraints_set_list\t----------------")
        total_constraints_set_list = total_constraints_set(model)
        for x in total_constraints_set_list:
            residual = abs(pyo.value(x.body) - pyo.value(x.lb))
            if residual > 1e-8:
                print(f"{x}\t{residual}")


def data_filter(data):

    for i, row in data.iterrows():
        # skip last 2 iterations
        if i >= 5:
            continue

        slope1 = (data[f"{source_name}_Y"][i + 1] - data[f"{source_name}_Y"][i]) / (
            data[f"{source_name}_X"][i + 1] - data[f"{source_name}_X"][i]
        )
        slope2 = (data[f"{source_name}_Y"][i + 2] - data[f"{source_name}_Y"][i]) / (
            data[f"{source_name}_X"][i + 2] - data[f"{source_name}_X"][i]
        )

        if slope1 < -1e-6 or slope2 < 1e-6:
            data = data.drop(i)

    return data


if __name__ == "__main__":
    main()
