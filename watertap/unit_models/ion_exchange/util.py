import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from pyomo.environ import (
    ConcreteModel,
    Expression,
    value,
    Var,
    Block,
    SolverFactory,
    units as pyunits,
)
import pyomo.contrib.parmest.parmest as parmest

from copy import deepcopy

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from .ion_exchange_clark import IonExchangeClark
from .ion_exchange_cphsdm import IonExchangeCPHSDM

theta_title_dict = {
    "thomas_constant": "k$_{Th}$",
    "resin_max_capacity": "q$_{max}$",
    "freundlich_k": "K$_{Fr}$",
    "freundlich_ninv": "1/n",
    "surf_diff_coeff": "D$_S$",
}


def plot_theta(blk):

    print("\n\nPLOT THETA\n\n")

    if not all(hasattr(blk, attr) for attr in ["test_cnorms", "test_bvs", "R_sq"]):
        blk.test_theta()

    stat_str = f"Theta (R$^2$: {blk.R_sq:.3f})"
    theta_str = [f"{theta_title_dict[k]}: {v:.2e}" for k, v in blk.theta_dict.items()]
    textstr = "\n".join([stat_str] + theta_str)

    # textstr = "\n".join(
    #     [
    #         f"Theta (R$^2$: {blk.R_sq:.3f})",
    #         f"  k$_T$: {blk.theta_dict['mass_transfer_coeff']:.2e}",
    #         f"  n: {blk.theta_dict['freundlich_n']:.2f}",
    #         f"  BV50: {round(blk.theta_dict['bv_50'])}",
    #     ]
    # )
    boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)

    fig1, ax1 = plt.subplots(figsize=blk.figsize)

    ax1.scatter(
        blk.bv_pred_theta_test,
        blk.cnorm_pred_theta_test,
        color="blue",
        marker=".",
        label="Predicted",
    )
    ax1.plot(
        blk.test_bvs,
        blk.test_cnorms,
        color="red",
        linestyle="-.",
        alpha=0.1,
        label="Test Data",
    )

    if len(blk.bv_fail_theta_test) != 0:
        ax1.scatter(
            blk.bv_fail_theta_test,
            blk.cnorm_fail_theta_test,
            marker="x",
            color="black",
            label="Infeasible",
        )

    fig2, ax2 = plt.subplots(figsize=blk.figsize)
    ax2.plot(
        blk.bv_pred_theta_test,
        blk.cnorm_pred_theta_test,
        color="blue",
        marker=".",
        alpha=0.1,
        # label="Predicted",
    )
    ax2.scatter(
        blk.bv_pred_theta,
        blk.cnorm_pred_theta,
        marker=".",
        color="green",
        label="Predicted",
    )
    ax2.plot(
        blk.keep_bvs,
        blk.keep_cnorms,
        marker=".",
        color="red",
        alpha=0.2,
        label="Parmest Data",
    )

    ax2.plot(
        blk.all_bvs,
        blk.all_cnorms,
        marker=".",
        color="purple",
        alpha=0.2,
        label="All Data",
    )

    if len(blk.bv_fail_theta) != 0:
        ax2.scatter(
            blk.bv_fail_theta,
            blk.cnorm_fail_theta,
            marker="x",
            color="black",
            label="Infeasible",
        )

    # ylabe = blk.component.swapcase() + f" [{blk.conc_units_str}]"
    ylabe = "C/C$_0$, " + blk.component.title()
    title = (
        f"PARMEST Results\nCurve {blk.curve_id}:\n"
        + blk.component
        + ", "
        + blk.ref.replace("_", " ").title()
        + ", "
        + blk.resin.replace("_", " ").title()
        + ", "
        + f"{blk.ebct_min} min EBCT"
    )
    for ax in [ax1, ax2]:
        ax.plot(
            [
                0,
                max(
                    blk.keep_bvs
                    + blk.bv_pred_theta_test
                    + blk.bv_pred_theta
                    + blk.all_bvs
                ),
            ],
            [
                # blk.c0
                1.0
                for _ in [
                    0,
                    max(
                        blk.keep_bvs
                        + blk.bv_pred_theta_test
                        + blk.bv_pred_theta
                        + blk.all_bvs
                    ),
                ]
            ],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
            color="black",
        )

        ax.set_xlabel("BV")
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.01, 1.02])
        ax.text(
            0.04,
            0.95,
            textstr,
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=boxprops,
        )
        ax.legend(loc="lower right")

        plt.tight_layout()
    blk.all_figs["theta_test"] = {"fig": fig1, "ax": ax1}
    blk.all_figs["theta_real"] = {"fig": fig2, "ax": ax2}

    fig3, ax3 = plt.subplots(figsize=blk.figsize)

    ax3.scatter(
        blk.keep_bv_theta,
        blk.bv_pred_theta,
        marker=".",
        color="green",
        label="Predicted",
    )
    ax3.plot(blk.keep_bvs, blk.keep_bvs, linestyle=":", color="red", alpha=0.25)
    ax3.set_xlabel("Actual BV")
    ax3.set_ylabel("Predicted BV")
    ax3.set_title(f"PARMEST Results\nCurve {blk.curve_id} Parity")
    plt.tight_layout()
    blk.all_figs["theta_parity"] = {"fig": fig3, "ax": ax3}

    fig4, ax4 = plt.subplots(figsize=blk.figsize)
    # blk.bv_error = [
    #     pred - actual for (pred, actual) in zip(blk.bv_pred_theta, blk.keep_bv_theta)
    # ]
    # blk.bv_rel_error = [
    #     (pred - actual) / actual
    #     for (pred, actual) in zip(blk.bv_pred_theta, blk.keep_bv_theta)
    # ]

    ax4.scatter(blk.keep_bv_theta, blk.bv_rel_error, color="black")
    ax4.set_xlabel("Actual BV")
    ax4.set_ylabel("Relative Error")
    ax4.set_title(f"PARMEST Results\nCurve {blk.curve_id} Relative Error")
    plt.tight_layout()

    blk.all_figs["theta_relerror"] = {"fig": fig4, "ax": ax4}
    fig4.show()


def build_results_dict(
    blk,
    ix_blk=None,
    split_str="fs.ix.",
    components=[Var, Expression],
    skips=[
        "_ref",
        "traps",
        "process_flow",
        "regeneration_stream",
        "c_norm_avg",
        "H2O",
    ],
    keeps=["flow_vol_phase", "conc_mass_phase_comp"],
):
    blk.results_dict = {}
    blk.v_name_map = {}
    blk.curve_deets = [
        "curve_id",
        "ref",
        # "component",
        "resin",
        # "resin_type",
        # "resin_func_group",
        "target_component",
        # "conc_units",
        "ebct_min",
        "ebct", 
        "flow_in",
        "charge",
        "c0",
        "expr_sf",
        "c0_min_thresh",
        "c0_max_thresh",
        "max_zero",
        "min_one",
        "c_next_thresh",
        "cb50_min_thresh",
        "autoscale_fixed",
        "tc",
        # "obj",
        "data_filter",
        # "flag"
    ]

    if ix_blk is None:
        ix_blk = blk._get_ix_blk()
    # elif isinstance(blk, ConcreteModel):
    #     ix_blk = blk._get_ix_blk(m=blk)

    # assert isinstance(blk, IonExchange0D)

    for deet in blk.curve_deets:
        if hasattr(blk, deet):
            blk.results_dict[deet] = []
        else:
            continue

    for c in components:
        for v in ix_blk.component_objects(c):
            # print(v.name)
            v_name = v.name.split(split_str)[1]
            if any(True if s in v.name else False for s in keeps):
                # print(f"KEEP {v.name}")
                if v.is_indexed():
                    # idx = [*v._index_set]
                    for i, b in v.items():
                        # for i in idx:
                        v_name_i = v_name + f"_{i}"
                        blk.results_dict[b.name] = []
                        blk.v_name_map[b.name] = v_name_i
                else:
                    blk.results_dict[v.name] = []
                    blk.v_name_map[v.name] = v_name
                continue

            if any(True if s in v.name else False for s in skips):
                # print(f"SKIP {v.name}")
                continue
            if v.is_indexed():
                idx = [*v._index_set]
                for i in idx:
                    # print(f'KEEP INDEXED {repr(c), i, v.name}')
                    v_name_i = v_name + f"_{i}"
                    blk.results_dict[v[i].name] = []
                    blk.v_name_map[v[i].name] = v_name_i
            else:
                blk.results_dict[v.name] = []
                blk.v_name_map[v.name] = v_name

    for k in blk.initial_guess_dict.keys():
        blk.results_dict[f"{k}_ig"] = []

    for theta in blk.thetas:
        blk.results_dict[f"{theta}_theta"] = []

    

    # for k in blk.theta_initial_guess.keys():
    #     blk.results_dict[f"{k}_ig_theta"] = []

    for k in blk.calc_from_constr_dict.keys():
        blk.results_dict[f"{k}_cfc"] = []
    
    if isinstance(blk.set_bounds_dict, dict): 
        for k in blk.set_bounds_dict.keys():
            blk.results_dict[f"{k}_lb"] = []
            blk.results_dict[f"{k}_ub"] = []
    else:
        blk.results_dict["set_bounds"] = list()

    blk.results_dict["flag"] = []
    blk.results_dict["obj"] = []
    blk.results_dict["data_filter"] = []


def results_dict_append(blk, ix_blk=None, components=[Var, Expression]):
    if ix_blk is None:
        ix_blk = blk._get_ix_blk()
    # elif isinstance(blk, ConcreteModel):
    #     blk = blk._get_ix_blk(m=blk)
    # else:
    #     pass

    # assert isinstance(blk, IonExchange0D)

    for deet in blk.curve_deets:
        # print(deet)
        if hasattr(blk, deet):
            blk.results_dict[deet].append(getattr(blk, deet))

    # for c in components:
        # if c is Var:
    for v in ix_blk.component_objects(components):
        if v.is_indexed():
            # idx = [*v._index_set]
            for _, b in v.items():
                # for i in idx:
                if b.name in blk.results_dict.keys():
                    blk.results_dict[b.name].append(value(b))
        else:
            if v.name in blk.results_dict.keys():
                # print(v.name, v())
                blk.results_dict[v.name].append(value(v))

    for k, v in blk.initial_guess_dict.items():
        blk.results_dict[f"{k}_ig"].append(v)

    # for k, v in blk.theta_initial_guess.items():
    #     blk.results_dict[f"{k}_ig_theta"].append(v)

    for k, v in blk.calc_from_constr_dict.items():
        blk.results_dict[f"{k}_cfc"].append(v)
    
    if isinstance(blk.set_bounds_dict, dict): 
        for k, v in blk.set_bounds_dict.items():
            blk.results_dict[f"{k}_lb"].append(v[0])
            blk.results_dict[f"{k}_ub"].append(v[1])
    else:
        blk.results_dict["set_bounds"].append("all")

    if hasattr(blk, "theta_dict"):
        for k, v in blk.theta_dict.items():
            blk.results_dict[f"{k}_theta"].append(v)
        blk.results_dict["obj"].append(blk.obj)

    # blk.results_dict["flag"].append(blk.flag)
    # blk.results_dict["data_filter"].append(blk.data_filter)

def make_results_df(blk): 

    blk.results_dict_save = deepcopy(blk.results_dict)

    for old_key, new_key in blk.v_name_map.items():
        blk.results_dict_save[new_key] = blk.results_dict_save.pop(old_key)

    try:
        blk.df_results = pd.DataFrame.from_dict(blk.results_dict_save)
    except ValueError:  # there are some empty lists in results_dict
        col_length = len(
            blk.results_dict_save["curve_id"]
        )  # curve_id will always be populated
        for k, v in blk.results_dict_save.items():
            if len(v) != col_length:
                blk.results_dict_save[k] = list(None for _ in range(col_length))
        blk.df_results = pd.DataFrame.from_dict(blk.results_dict_save)

def save_results(blk, overwrite=True, results_filename=None):

    if blk.save_directory is None:
        raise ValueError("Must provide save directory to save results.")
    if not hasattr(blk, "df_results"):
        blk.make_results_df()
        
    if results_filename is None:
        blk.results_filename = (
            f"{blk.save_directory}/results/curve{blk.curve_id}_results.csv"
        )
    else:
        blk.results_filename = results_filename
    if overwrite:
        blk.df_results.to_csv(blk.results_filename, index=False)
    elif not overwrite:
        file_base = blk.results_filename.split(".csv")[0]
        append = 2
        while os.path.exists(blk.results_filename):
            blk.results_filename = file_base + f"_{append}.csv"
            append += 1
        blk.df_results.to_csv(blk.results_filename, index=False)


def save_figs(blk, overwrite=True, extension=None):

    if blk.save_directory is None:
        raise ValueError("Must provide save directory to save figures.")

    for figname, d in blk.all_figs.items():
        fig = d["fig"]
        fig_file = (
            f"{blk.save_directory}/figs/{figname}/curve{blk.curve_id}_{figname}.png"
        )
        if overwrite:
            fig.savefig(fig_file, bbox_inches="tight")
        elif not overwrite:
            x = fig_file.split(".png")[0]
            append = 2
            while os.path.exists(fig_file):
                fig_file = x + f"_{append}.png"
                append += 1
            fig.savefig(fig_file, bbox_inches="tight")

def make_output_df(blk):

    if blk.save_directory is None:
        raise ValueError("Must provide save directory to save output.")

    datas = [
        "bv50_bvs",
        "bv50_cbs",
        "bv_pred_ig_orig",
        "bv_pred_ig",
        "cb_pred_ig_orig",
        "cb_pred_ig",
        "keep_bvs",
        "keep_cbs",
        "excl_bvs",
        "excl_cbs",
        "all_bvs",
        "all_cbs",
        "keep_cnorms",
        "excl_cnorms",
        "bv_fail_ig",
        "cb_fail_ig",
        "bv_pred_theta_test",
        "cb_pred_theta_test",
        "bv_fail_theta_test",
        "cb_fail_theta_test",
        "bv_pred_theta",
        "bv_fail_theta",
        "keep_bv_theta",
        "cb_pred_theta",
        "cb_fail_theta",
        "bv_error",
        "bv_rel_error",
        "bv_abs_error",
        "bv_abs_rel_error",
    ]

    scalar_datas = [
        "c0_max_thresh",
        "c0_min_thresh", 
        "mean_absolute_error", 
        "sum_squared_error",
        "expr_sf", 
        "obj"]


    if hasattr(blk, "theta_dict"):
        tmp_dict = dict()
        tmp_dict[blk.curve_id] = blk.theta_dict
        df_theta = pd.DataFrame.from_dict(tmp_dict).T
        # df_theta["loading_rate"] = blk.loading_rate
        df_theta["curve_id"] = blk.curve_id
        df_theta["c0"] = blk.c0
        df_theta["target_component"] = blk.target_component
        df_theta["ref"] = blk.ref
        df_theta["sum_squared_error"] = blk.sum_squared_error
        df_theta["mean_absolute_error"] = blk.mean_absolute_error
        df_theta["obj"] = blk.obj
        df_theta["expr_sf"] = blk.expr_sf
        blk.df_theta = df_theta

    tmps = []
    for data in sorted(datas):
        # print(data)
        if hasattr(blk, data):
            tmp = pd.DataFrame.from_dict({data: getattr(blk, data)})
            tmps.append(tmp)

    blk.df_output = pd.concat(tmps, ignore_index=False, axis=1)

    for k, v in blk.initial_guess_dict.items():
        blk.df_output[f"{k}_ig"] = v

    for data in sorted(scalar_datas):
        print(data)
        if hasattr(blk, data):
            blk.df_output[data] = getattr(blk, data)
    # blk.df_output["c0_max_thresh"] = blk.c0_max_thresh
    # blk.df_output["c0_min_thresh"] = blk.c0_min_thresh
    # blk.df_output["use_all_data"] = blk.use_all_data
    # blk.df_output["use_this_data"] = blk.use_this_data

    for deet in blk.curve_deets:
        if hasattr(blk, deet):
            blk.df_output[deet] = getattr(blk, deet)

    if hasattr(blk, "theta_dict"):
        for k, v in blk.theta_dict.items():
            blk.df_output[f"{k}_theta"] = v
        blk.df_output["obj"] = blk.obj

def save_output(blk, overwrite=True):

    if not hasattr(blk, "df_output"):
        blk.make_output_df()

    output_file_base = f"{blk.save_directory}/output/curve{blk.curve_id}_output.csv"
    theta_file_base = f"{blk.save_directory}/theta/curve{blk.curve_id}_theta.csv"



    blk.df_theta.to_csv(theta_file_base, index=False)

    if overwrite:
        blk.df_output.to_csv(output_file_base, index=False)
    elif not overwrite:
        x = output_file_base.split(".csv")[0]
        append = 2
        while os.path.exists(output_file_base):
            output_file_base = x + f"_{append}.csv"
            append += 1
        blk.df_output.to_csv(output_file_base, index=False)

        

def plot_initial_guess(blk):

    
    # textstr = "\n".join([stat_str] + theta_str)

    fig, ax = plt.subplots(figsize=blk.figsize)
    ax.scatter(
        blk.bv_pred_ig,
        blk.cb_pred_ig,
        marker="o",
        color="green",
        label="Predicted",
    )
    ax.scatter(
        blk.bv_pred_ig_fake,
        blk.cb_pred_ig_fake,
        marker=".",
        color="blue",
        alpha=0.1,
        # label="Predicted",
    )
    if len(blk.bv_fail_ig) != 0:
        ax.scatter(
            blk.bv_fail_ig,
            blk.cb_fail_ig,
            marker="x",
            color="black",
            label="Infeasible",
        )

    ax.plot(
        blk.keep_bvs,
        blk.keep_cnorms,
        marker=".",
        color="red",
        alpha=0.25,
        label="Parmest Data",
    )
    ax.plot(
        blk.all_bvs,
        blk.all_cnorms,
        marker=".",
        color="purple",
        alpha=0.25,
        label="All Data",
    )
    ax.plot(
        [0, max(blk.bv_pred_ig + blk.bv_pred_ig_fake + blk.keep_bvs)],
        [1, 1],
        linestyle="-.",
        alpha=0.25,
        # label="Influent",
    )
    if blk.ix_model is IonExchangeClark:
        ax.plot(
            [0, max(blk.bv_pred_ig + blk.bv_pred_ig_fake + blk.keep_bvs)],
            [0.5, 0.5],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
        )

    title = (
        f"Initial Guess - Curve {blk.curve_id}:\n"
        + blk.component
        + " "
        + blk.ref.replace("_", " ").title()
        + " "
        + blk.resin.replace("_", " ").title()
        + f" EBCT = {blk.ebct_min} min"
    )
    ax.set_xlabel("BV")
    # ylabe = blk.component.swapcase() + f" [{blk.conc_units_str}]"
    ylabe = "C/C$_0$" + blk.component
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    ax.set_ylim([-0.01, 1.02])

    if blk.ix_model is IonExchangeClark:

        ax.scatter(
            [blk.initial_guess_dict["bv_50"]],
            [0.5],
            marker="*",
            color="springgreen",
            # label="BV50 Guess",
        )
        textstr = "\n".join(
            [
                f"Guess",
                f"  k$_T$: {blk.initial_guess_dict['mass_transfer_coeff']:.2e}",
                f"  $n$: {blk.initial_guess_dict['freundlich_n']:.2f}",
                f"  BV50: {round(blk.initial_guess_dict['bv_50'])}",
            ]
        )
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)
        ax.text(
            0.04,
            0.95,
            textstr,
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=boxprops,
        )
    if blk.ix_model is IonExchangeCPHSDM:
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)
        theta_str = [f"{theta_title_dict[k]}: {v:.2e}" for k, v in blk.initial_guess_dict.items()]
        theta_str = "\n".join(["Initial Guess:"] + theta_str)
        ax.text(
            0.04,
            0.95,
            theta_str,
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=boxprops,
        )
    ax.legend()
    plt.tight_layout()

    blk.ig_fig = fig

    blk.all_figs["initial_guess"] = {"fig": fig, "ax": ax}


def plot_estimate_bv50(blk):

    if blk.bv50_est < 10:
        blk.plot_initial_guess()
        return
    max_bv = max([blk._linear_inv(1)] + blk.all_bvs)
    cb_50 = 0.5

    title = (
        f"BV50 Estimate Results Compare - Curve {blk.curve_id}:\n"
        + blk.component.title()
        + " "
        + blk.ref.replace("_", " ").title()
        + " "
        + blk.resin.replace("_", " ").title()
        + f", EBCT = {blk.ebct_min} min"
    )
    try:
        textstr = "\n".join(
            [
                f"Guess",
                f"  k$_T$: {blk.initial_guess_dict['mass_transfer_coeff']:.2e}",
                f"  n: {blk.initial_guess_dict['freundlich_n']:.2f}",
                f"  BV50 (new): {round(blk.initial_guess_dict['bv_50'])}",
                f"  BV50 (orig): {round(blk.bv50_orig)}",
            ]
        )
    except:
        textstr = "\n".join(
            [
                f"Guess",
                f"  k$_T$: {blk.m0.fs.ix.mass_transfer_coeff():.2e}",
                f"  n: {blk.initial_guess_dict['freundlich_n']:.2f}",
                f"  BV50 (new): {round(blk.initial_guess_dict['bv_50'])}",
                f"  BV50 (orig): {round(blk.bv50_orig)}",
            ]
        )
    boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)
    fig, ax = plt.subplots(figsize=blk.figsize)
    ax.plot(
        blk.keep_bvs,
        blk.keep_cnorms,
        marker=".",
        color="red",
        alpha=0.2,
        label="Parmest Data",
    )
    ax.plot(
        blk.all_bvs,
        blk.all_cnorms,
        marker=".",
        color="purple",
        alpha=0.2,
        label="All Data",
    )
    ax.scatter(
        blk.bv_pred_ig_orig,
        blk.cb_pred_ig_orig,
        marker="o",
        color="green",
        label="Predicted (Orig)",
    )
    ax.scatter(
        blk.bv_pred_ig,
        blk.cb_pred_ig,
        marker="o",
        color="blue",
        label="Predicted (New)",
    )
    ax.scatter(
        blk.bv_pred_ig_fake,
        blk.cb_pred_ig_fake,
        marker=".",
        color="blue",
        alpha=0.1,
        # label="Predicted",
    )
    if len(blk.bv_fail_ig) != 0:
        ax.scatter(
            blk.bv_fail_ig,
            blk.cb_fail_ig,
            marker="x",
            color="black",
            label="Infeasible",
        )

    ax.plot(
        [0, max_bv],
        [1, 1],
        linestyle="-.",
        color="grey",
        alpha=0.25,
        # label="Influent",
    )
    ax.plot(
        [0, max_bv],
        [cb_50 for _ in [0, max_bv]],
        linestyle="-.",
        color="grey",
        alpha=0.25,
        # label="50% Influent",
    )

    xs = np.linspace(min(blk.bv50_bvs), max_bv, len(blk.bv50_bvs))
    ys = [blk._linear_fit(x) for x in xs]

    ax.plot(
        xs,
        ys,
        marker=".",
        color="purple",
        linestyle=":",
        alpha=0.5,
        label="Linear Fit",
    )
    ax.scatter(
        [blk.bv50_orig],
        [cb_50],
        marker="*",
        color="springgreen",
        label="BV50 Guess (Orig)",
    )
    ax.scatter(
        [blk.bv50_est],
        [cb_50],
        marker="*",
        color="fuchsia",
        label="BV50 Guess (New)",
    )

    ax.text(
        0.04,
        0.95,
        textstr,
        verticalalignment="top",
        transform=ax.transAxes,
        bbox=boxprops,
    )
    ax.set_xlabel("BV")
    ylabe = blk.component.title() + f" [{blk.conc_units_str}]"
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    ax.set_ylim([-0.01, 1.02])
    ax.set_xlim([0, max_bv])
    ax.legend(loc="lower right")
    plt.tight_layout()

    blk.bv50_fig = fig
    blk.all_figs["bv_50_est"] = {"fig": fig, "ax": ax}


def plot_curve(blk):

    modified_dropped = blk.dropped_points + blk.modified_points

    fig, ax = plt.subplots(figsize=blk.figsize)

    # ax.scatter(blk.)
    ax.plot(
        blk.keep_bvs,
        blk.keep_cnorms,
        color="red",
        marker=".",
        alpha=0.25,
        label="Parmest Data",
    )
    ax.plot(
        blk.all_bvs,
        blk.all_cnorms,
        color="purple",
        marker=".",
        alpha=0.25,
        label="Filtered Data",
    )
    ax.plot(
        [0, max(blk.all_bvs)],
        [1, 1],
        linestyle="-.",
        alpha=0.25,
        color="navy",
        # label="Influent",
    )
    ax.scatter(
        blk.excl_bvs, blk.excl_cnorms, marker="x", color="k", label="Excluded Data"
    )

    if blk.just_plot_curve:
        point_labels = blk.filtered_data.reset_index().index.to_list()
        for i, bv, cnorm in zip(point_labels, blk.filtered_data.bv.to_list(), blk.filtered_data.c_norm.to_list()):
            # ax.annotate(f"{cnorm:.3f}\n{i}")
            ax.text(bv, cnorm, i)

    # if len(modified_dropped) > 0:
    #     md = blk.input_data.loc[modified_dropped]
    #     ax.plot(md.bv, md.c_norm, marker="^", color="c", alpha=0.25, label="Modified/Dropped Data")
    #     for p in modified_dropped:
    #         if p in blk.dropped_points:
    #             print(f"Input data with point id {p} was dropped.")
    #         if p in blk.modified_points:
    #             print(f"Input data with point id {p} was modified:")
    #             s = blk.input_data.loc[p]
    #             f = blk.filtered_data.loc[p]
    #             e = s.eq(f)
    #             for i, b in zip(s.index, e):
    #                 if not b:
    #                     print(f"\t{i}: input {s[i]} --> filtered {f[i]}")

    # for c in i.

    title = (
        f"Data Used for Curve {blk.curve_id}\n"
        + blk.component.title()
        + ", "
        + blk.ref.replace("_", " ").title()
        + ", "
        + blk.resin.replace("_", " ").title()
        + f", EBCT = {blk.ebct_min} min"
    )
    ax.set_xlabel("BV")
    # ylabe = blk.component.swapcase()
    ylabe = "C/C$_0$, " + blk.component.title()
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    # ax.set_ylim([-0.05, 1.02])
    # ax.set_xlim([-5, max(blk.all_bvs) * 1.05])
    ax.legend()
    plt.tight_layout()

    blk.curve_fig = fig

    blk.all_figs["raw_curve"] = {"fig": fig, "ax": ax}


# def _just_plot_curve(blk):
#     tmp = blk.df_curve.reset_index()
#     _, ax = plt.subplots(figsize=blk.figsize)

#     ax.plot(tmp.bv, tmp.c_norm, marker=".")
#     ax.plot([0, tmp.bv.max()], [1, 1], linestyle=":")
#     ax.set_title(
#         f"Curve {blk.curve_id}: {blk.ref.replace('_', ' ').title()}, {blk.component.upper()}"
#     )

#     # ax.set_ylim([-0.01, 1.02])
#     ax.set_xlabel("BV")
#     ax.set_ylabel(f"C/C$_0$ [{blk.component.upper()}]")
#     blk.idx_labels = zip(
#         tmp.bv.to_list(), tmp.c_norm.to_list(), tmp.index.to_list()
#     )
#     for x, y, i in blk.idx_labels:
#         # ax.annotate(f"{i}\n{y}", (x, y), fontsize="small", rotation=-45)
#         ax.annotate(i, (x, y), fontsize="small")
#     plt.tight_layout()
#     print(f"INDEXES: {tmp.index.to_list()}")
