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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyomo.environ as pyo
import pyomo.contrib.parmest.parmest as parmest
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
import idaes.core.util.model_statistics as istat
import os as os

from math import floor
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.models.unit_models import Feed
from watertap.unit_models.gac import GAC
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    DiffusivityCalculation,
)

__author__ = "Hunter Barber"

solver = get_solver()
log = idaeslog.getSolveLogger("solver.demo")
log.setLevel(idaeslog.DEBUG)
pyo.SolverFactory.register("ipopt")(pyo.SolverFactory.get_class("ipopt-watertap"))


def main():

    global data, source_name, source_name_list, species_name, guess_freund_k, guess_freund_ninv, guess_ds, min_st_surrogate, throughput_surrogate

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
    species_name = "PFOA"
    data = pd.read_csv(
        "watertap/examples/flowsheets/pfas_treatment_analysis/gac_adsorption_data/F400_PFOA.csv"
    )

    # load surrogate objects
    min_st_surrogate = PysmoSurrogate.load_from_file(
        "watertap/examples/flowsheets/pfas_treatment_analysis/gac_surrogate/trained_surrogate_models/min_st_pysmo_surr_spline.json",
    )
    throughput_surrogate = PysmoSurrogate.load_from_file(
        "watertap/examples/flowsheets/pfas_treatment_analysis/gac_surrogate/trained_surrogate_models/throughput_pysmo_surr_linear.json",
    )

    # initial guess for regressed variables
    guess_freund_k = 10
    guess_freund_ninv = 0.8
    guess_ds = 1e-15

    terminal_len = os.get_terminal_size().columns
    dict_resolve_results, data_filtered, df_param_results = {}, {}, {}
    for source_name in source_name_list:

        middle_str = f"species: {species_name} / source: {source_name}"
        middle_len = len(middle_str) + 2
        end_len = floor((terminal_len - middle_len) / 2)
        end_str = "=" * end_len
        print(end_str, middle_str, end_str)
        # ---------------------------------------------------------------------
        # model regression

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
        expr_sf = 1e-6

        def SSE(model, data_filtered_case):
            expr = (
                float(data_filtered_case[f"{source_name}_X"])
                - model.fs.gac.bed_volumes_treated
            ) ** 2
            return expr * expr_sf

        middle_str = f"model regression"
        middle_len = len(middle_str) + 2
        end_len = floor((terminal_len - middle_len) / 2)
        end_str = "-" * end_len
        print(end_str, middle_str, end_str)
        # Create an instance of the parmest estimator
        pest = parmest.Estimator(
            parmest_regression,
            data_filtered_case,
            theta_names,
            SSE,
            tee=False,
            diagnostic_mode=False,
            solver_options={"bound_push": 1e-8},
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

        df_param_results[f"{source_name}"] = {
            "freund_k": theta[0],
            "freund_ninv": theta[1],
            "ds": theta[2],
        }

        middle_str = f"model resolve"
        middle_len = len(middle_str) + 2
        end_len = floor((terminal_len - middle_len) / 2)
        end_str = "-" * end_len
        print(end_str, middle_str, end_str)
        # rebuild model across CP to view regression results
        # theta = [guess_freund_k, guess_freund_ninv, guess_ds]
        dict_iter = solve_regression(theta)

        dict_resolve_results[f"{source_name}"] = dict_iter

    # save regression data
    df_regression_results = pd.DataFrame.from_dict(df_param_results)
    df_regression_results.to_csv(
        "watertap/examples/flowsheets/pfas_treatment_analysis/regression_results.csv"
    )

    middle_str = f"plotting regression to figure"
    middle_len = len(middle_str) + 2
    end_len = floor((terminal_len - middle_len) / 2)
    end_str = "=" * end_len
    print(end_str, middle_str, end_str)
    plot_regression(dict_resolve_results, data_filtered)
    print(df_param_results)


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
    m.fs.gac = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="calculated",
        surface_diffusion_coefficient_type="fixed",
        finite_elements_ss_approximation=5,
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
    deactivate_ss_calculations(m)
    # scaling
    model_scale(m)
    # initialization
    model_init(m, outlvl=idaeslog.ERROR)
    # switch to surrogate
    activate_surrogate(m)
    # re-fix to conc_ratio spec from data
    m.fs.gac.conc_ratio_replace.fix(float(data[f"{source_name}_Y"]))

    return m


def solve_regression(theta):

    # build model
    m = model_build()
    deactivate_ss_calculations(m)
    # scaling
    model_scale(m)
    # initialization
    model_init(m, outlvl=idaeslog.ERROR)
    # switch to surrogate input
    activate_surrogate(m)

    # refix differing parameters
    m.fs.gac.freund_k.fix(theta[0])
    m.fs.gac.freund_ninv.fix(theta[1])
    m.fs.gac.ds.fix(theta[2])

    # check profile against pilot data
    conc_ratio_input = np.linspace(0.005, 0.95, 50)
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

    dict_iter = {
        "conc_ratio": conc_ratio_list,
        "bed_volumes_treated": bed_volumes_treated_list,
    }

    return dict_iter


def plot_regression(dict_resolve_results, data_filtered):

    color_code = [
        "#ffc000",
        "#ed7d31",
        "#c00000",
        "#00b050",
        "#7030a0",
        "#0070c0",
        "#bdd7ee",
        "#000000",
        "#a5a5a5",
        "#f8cbad",
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
            mfc="None",
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
            dict_resolve_results[source_name_list[i]]["bed_volumes_treated"],
            dict_resolve_results[source_name_list[i]]["conc_ratio"],
            color_code[i],
            label="regression results",
        )
        ax.legend(
            title=f"{source_name_list[i]}",
            fontsize="x-small",
            loc="lower right",
        )

        if i == 8:
            ax.set_xlabel("bed volumes treated")
            ax.set_ylabel("concentration ratio")

        i = i + 1

    fig1.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0, hspace=0)

    plt.show()

    # fig2
    plt.figure(2)
    for j in range(10):

        plt.plot(
            dict_resolve_results[source_name_list[j]]["bed_volumes_treated"],
            dict_resolve_results[source_name_list[j]]["conc_ratio"],
            color_code[j],
            label=f"{source_name_list[j]}",
        )
        plt.plot(
            data[f"{source_name_list[j]}_X"],
            data[f"{source_name_list[j]}_Y"],
            "o",
            mec=color_code[j],
            mfc="None",
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
        fontsize="x-small",
        loc="lower right",
    )

    plt.show()


def model_scale(model):

    iscale.set_scaling_factor(model.fs.gac.bed_length, 1e2)
    iscale.set_scaling_factor(model.fs.gac.bed_diameter, 1e3)
    iscale.set_scaling_factor(model.fs.gac.bed_area, 1e5)
    iscale.set_scaling_factor(model.fs.gac.bed_volume, 1e6)
    iscale.set_scaling_factor(model.fs.gac.bed_mass_gac, 1e3)
    iscale.set_scaling_factor(model.fs.gac.mass_adsorbed, 1e9)
    iscale.set_scaling_factor(model.fs.gac.gac_usage_rate, 1e10)

    iscale.calculate_scaling_factors(model)


def model_init(model, outlvl):

    model.fs.gac.bed_length = 1e-2
    model.fs.gac.bed_diameter = 1e-3
    model.fs.gac.bed_area = 1e-5
    model.fs.gac.bed_volume = 1e-6
    model.fs.gac.bed_mass_gac = 1e-3
    model.fs.gac.mass_adsorbed = 1e-9
    model.fs.gac.gac_usage_rate = 1e-10

    optarg = solver.options
    model.fs.feed.initialize(optarg=optarg, outlvl=outlvl)
    propagate_state(model.fs.s01)
    model.fs.gac.initialize(optarg=optarg, outlvl=outlvl)


def deactivate_ss_calculations(m):

    # deactivate steady state equations
    m.fs.gac.eq_ele_throughput[:].deactivate()
    m.fs.gac.eq_ele_min_operational_time[:].deactivate()
    m.fs.gac.eq_ele_conc_ratio_replace[:].deactivate()
    m.fs.gac.eq_ele_operational_time[:].deactivate()
    m.fs.gac.eq_ele_conc_ratio_term[:].deactivate()
    m.fs.gac.eq_conc_ratio_avg.deactivate()

    # fix variables used in state equations
    m.fs.gac.ele_throughput[:].fix()
    m.fs.gac.ele_min_operational_time[:].fix()
    m.fs.gac.ele_conc_ratio_replace[:].fix()
    m.fs.gac.ele_operational_time[:].fix()
    m.fs.gac.ele_conc_ratio_term[:].fix()

    # fix conc_ratio_avg for mass balance
    m.fs.gac.conc_ratio_avg.fix(0.1)


def activate_surrogate(m):

    # deactivate empirical equations equations
    m.fs.gac.eq_min_number_st_cps.deactivate()
    m.fs.gac.eq_throughput.deactivate()

    # establish surrogates
    m.fs.min_st_surrogate = SurrogateBlock(concrete=True)
    m.fs.min_st_surrogate.build_model(
        min_st_surrogate,
        input_vars=[m.fs.gac.freund_ninv, m.fs.gac.N_Bi],
        output_vars=[m.fs.gac.min_N_St],
    )
    m.fs.throughput_surrogate = SurrogateBlock(concrete=True)
    m.fs.throughput_surrogate.build_model(
        throughput_surrogate,
        input_vars=[m.fs.gac.freund_ninv, m.fs.gac.N_Bi, m.fs.gac.conc_ratio_replace],
        output_vars=[m.fs.gac.throughput],
    )


def model_solve(model, solver_log=True):

    # check model
    assert_units_consistent(model)  # check that units are consistent
    assert istat.degrees_of_freedom(model) == 0

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

    # log problems of non-optimal solve
    if not term_cond == "optimal":
        badly_scaled_var_list = iscale.badly_scaled_var_generator(model)
        print("------------------      badly_scaled_var_list       ------------------")
        for x in badly_scaled_var_list:
            print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")
        print("------------------    variables_near_bounds_list    ------------------")
        variables_near_bounds_list = istat.variables_near_bounds_generator(model)
        for x in variables_near_bounds_list:
            print(f"{x.name}\t{x.value}")
        print("------------------    total_constraints_set_list    ------------------")
        istat.activated_constraints_set_list = istat.activated_constraints_set(model)
        for x in istat.activated_constraints_set_list:
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
