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

theta_title_dict = {
    "thomas_constant": "k$_{Th}$",
    "resin_max_capacity": "q$_{max}$"
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

    # ylabe = blk.compound.swapcase() + f" [{blk.conc_units_str}]"
    ylabe = "C/C$_0$, " + blk.compound.title()
    title = (
        f"Curve {blk.curve_id}:\n"
        + blk.compound.title()
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
    ax3.set_title(f"Curve {blk.curve_id} Parity")
    plt.tight_layout()
    blk.all_figs["theta_parity"] = {"fig": fig3, "ax": ax3}

    fig4, ax4 = plt.subplots(figsize=blk.figsize)
    blk.bv_error = [
        pred - actual for (pred, actual) in zip(blk.bv_pred_theta, blk.keep_bv_theta)
    ]
    blk.bv_rel_error = [
        (pred - actual) / actual
        for (pred, actual) in zip(blk.bv_pred_theta, blk.keep_bv_theta)
    ]

    ax4.scatter(blk.keep_bv_theta, blk.bv_rel_error, color="black")
    ax4.set_xlabel("Actual BV")
    ax4.set_ylabel("Relative Error")
    ax4.set_title(f"Curve {blk.curve_id} Relative Error")
    plt.tight_layout()

    blk.all_figs["theta_relerror"] = {"fig": fig4, "ax": ax4}
    fig4.show()


def build_results_dict(
    blk,
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
    check_keys=False,
):
    blk.results_dict = {}
    blk.v_name_map = {}
    blk.curve_deets = [
        "curve_id",
        "ref",
        "compound",
        "resin",
        "resin_type",
        "resin_func_group",
        "target_component",
        "conc_units",
        "ebct_min",
        "flow_in",
        "charge",
        "c0",
        # "expr_sf",
        "c0_min_thresh",
        "c0_max_thresh",
        "c_next_thresh",
        "cb50_min_thresh",
        "autoscale_fixed",
        "tc",
    ]

    if blk is None:
        blk = blk._get_ix_blk()
    elif isinstance(blk, ConcreteModel):
        blk = blk._get_ix_blk(m=blk)

    assert isinstance(blk, blk.ix_model)

    for deet in blk.curve_deets:
        blk.results_dict[deet] = []

    for c in components:
        for v in blk.component_objects(c):
            # print(v.name)
            v_name = v.name.split(split_str)[1]
            if any(True if s in v.name else False for s in keeps):
                # print(f"KEEP {v.name}")
                if v.is_indexed():
                    idx = [*v._index_set]
                    for i in idx:
                        v_name_i = v_name + f"_{i}"
                        blk.results_dict[v[i].name] = []
                        blk.v_name_map[v[i].name] = v_name_i
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

    for k in blk.initial_guess_theta_dict.keys():
        blk.results_dict[f"{k}_ig_theta"] = []

    for k in blk.calc_from_constr_dict.keys():
        blk.results_dict[f"{k}_cfc"] = []

    blk.results_dict["flag"] = []

    if check_keys:
        for k in blk.results_dict.keys():
            print(k)


def results_dict_append(
    blk=None,
    components=[Var, Expression],
    tmp_results_dict=None,
):

    if blk is None:
        blk = blk.ix
    elif isinstance(blk, ConcreteModel):
        blk = blk._get_ix_blk(m=blk)
    else:
        pass

    assert isinstance(blk, blk.ix_model)

    for deet in blk.curve_deets:
        # print(deet)
        blk.results_dict[deet].append(getattr(blk, deet))

    if tmp_results_dict is None:
        for c in components:
            # if c is Var:
            for v in blk.component_objects(c):
                if v.is_indexed():
                    idx = [*v._index_set]
                    for i in idx:
                        if v[i].name in blk.results_dict.keys():
                            blk.results_dict[v[i].name].append(v[i]())
                else:
                    if v.name in blk.results_dict.keys():
                        # print(v.name, v())
                        blk.results_dict[v.name].append(v())

        for k, v in blk.initial_guess_dict.items():
            blk.results_dict[f"{k}_ig"].append(v)

        for k, v in blk.initial_guess_theta_dict.items():
            blk.results_dict[f"{k}_ig_theta"].append(v)

        for k, v in blk.calc_from_constr_dict.items():
            blk.results_dict[f"{k}_cfc"].append(v)

        if hasattr(blk, "theta_dict"):
            for k, v in blk.theta_dict.items():
                blk.results_dict[f"{k}_theta"].append(v)
            blk.results_dict["obj"].append(blk.obj)

        blk.results_dict["flag"].append(blk.flag)

    else:
        for k, v in tmp_results_dict.items():
            blk.results_dict[k].append(v)


def save_results(blk, overwrite=False, results_filename=None):

    blk.results_dict_save = deepcopy(blk.results_dict)

    for old_key, new_key in blk.v_name_map.items():
        blk.results_dict_save[new_key] = blk.results_dict_save.pop(old_key)

    if results_filename is None:
        blk.results_filename = f"results/curve_id{blk.curve_id}_results.csv"
    else:
        blk.results_filename = results_filename

    blk.df_results = pd.DataFrame.from_dict(blk.results_dict_save)

    if overwrite:
        blk.df_results.to_csv(blk.results_filename, index=False)
    elif not overwrite:
        x = blk.results_filename.split(".csv")[0]
        append = 2
        while os.path.exists(blk.results_filename):
            blk.results_filename = x + f"_{append}.csv"
            append += 1
        blk.df_results.to_csv(blk.results_filename, index=False)


def save_figs(blk, overwrite=False, extension=None):

    # fig_file_base = f"figs/curve_id{blk.curve_id}_"

    for figname, d in blk.all_figs.items():
        # fig_file_base = f"figs/{figname}/curve_id{blk.curve_id}_"
        fig = d["fig"]
        fig_file = f"figs/{figname}/curve_id{blk.curve_id}_{figname}.png"
        if overwrite:
            fig.savefig(fig_file, bbox_inches="tight")
        elif not overwrite:
            x = fig_file.split(".png")[0]
            append = 2
            while os.path.exists(fig_file):
                fig_file = x + f"_{append}.png"
                append += 1
            fig.savefig(fig_file, bbox_inches="tight")


def save_output(blk, overwrite=False):

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
    ]

    data_file_base = f"output/curve_id{blk.curve_id}_output.csv"
    theta_file_base = f"theta/curve_id{blk.curve_id}_theta.csv"

    if hasattr(blk, "theta_dict"):
        tmp_dict = dict()
        tmp_dict[blk.curve_id] = blk.theta_dict
        df_theta = pd.DataFrame.from_dict(tmp_dict).T
        df_theta["loading_rate"] = blk.loading_rate
        df_theta["curve_id"] = blk.curve_id
        df_theta["c0"] = pyunits.convert(
            blk.c0 * blk.conc_units, to_units=pyunits.kg / pyunits.m**3
        )()
        df_theta.to_csv(theta_file_base, index=False)

    tmps = []
    for data in sorted(datas):
        if hasattr(blk, data):
            tmp = pd.DataFrame.from_dict({data: getattr(blk, data)})
            tmps.append(tmp)

    blk.df_data = pd.concat(tmps, ignore_index=False, axis=1)

    for k, v in blk.initial_guess_dict.items():
        blk.df_data[f"{k}_ig"] = v

    blk.df_data["c0_max_thresh"] = blk.c0_max_thresh
    blk.df_data["c0_min_thresh"] = blk.c0_min_thresh

    for deet in blk.curve_deets:
        if hasattr(blk, deet):
            blk.df_data[deet] = getattr(blk, deet)

        for k, v in blk.theta_dict.items():
            blk.df_data[f"{k}_theta"] = v
        blk.df_data["obj"] = blk.obj

    if overwrite:
        blk.df_data.to_csv(data_file_base, index=False)
    elif not overwrite:
        x = data_file_base.split(".csv")[0]
        append = 2
        while os.path.exists(data_file_base):
            data_file_base = x + f"_{append}.csv"
            append += 1
        blk.df_data.to_csv(data_file_base, index=False)


def plot_initial_guess(blk):

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
    ax.plot(
        [0, max(blk.bv_pred_ig + blk.bv_pred_ig_fake + blk.keep_bvs)],
        [0.5, 0.5],
        linestyle="-.",
        alpha=0.25,
        # label="Influent",
    )



    title = (
        f"Initial Guess - Curve {blk.curve_id}:\n"
        + blk.compound.title()
        + " "
        + blk.ref.replace("_", " ").title()
        + " "
        + blk.resin.replace("_", " ").title()
        + f" EBCT = {blk.ebct_min} min"
    )
    ax.set_xlabel("BV")
    # ylabe = blk.compound.swapcase() + f" [{blk.conc_units_str}]"
    ylabe = "C/C$_0$" + blk.compound.title()
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
        + blk.compound.title()
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
    ylabe = blk.compound.title() + f" [{blk.conc_units_str}]"
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    ax.set_ylim([-0.01, 1.02])
    ax.set_xlim([0, max_bv])
    ax.legend(loc="lower right")
    plt.tight_layout()

    blk.bv50_fig = fig
    blk.all_figs["bv_50_est"] = {"fig": fig, "ax": ax}


def plot_curve(blk):

    fig, ax = plt.subplots(figsize=blk.figsize)
    ax.scatter(blk.excl_bvs, blk.excl_cnorms, marker="x", color="k", label="Excluded data")
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
        label="All Data",
    )
    ax.plot(
        [0, max(blk.all_bvs)],
        [1, 1],
        linestyle="-.",
        alpha=0.25,
        label="Influent",
    )

    title = (
        f"Data Used for Curve {blk.curve_id}\n"
        + blk.compound.title()
        + ", "
        + blk.ref.replace("_", " ").title()
        + ", "
        + blk.resin.replace("_", " ").title()
        + f", EBCT = {blk.ebct_min} min"
    )
    ax.set_xlabel("BV")
    # ylabe = blk.compound.swapcase()
    ylabe = "C/C$_0$, " + blk.compound.title()
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    # ax.set_ylim([-0.05, 1.02])
    ax.set_xlim([-5, max(blk.all_bvs) * 1.05])
    ax.legend()
    plt.tight_layout()

    blk.curve_fig = fig

    blk.all_figs["raw_curve"] = {"fig": fig, "ax": ax}


# def _just_plot_curve(self):
#     tmp = self.df_curve.reset_index()
#     _, ax = plt.subplots(figsize=self.figsize)

#     ax.plot(tmp.bv, tmp.c_norm, marker=".")
#     ax.plot([0, tmp.bv.max()], [1, 1], linestyle=":")
#     ax.set_title(
#         f"Curve {self.curve_id}: {self.ref.replace('_', ' ').title()}, {self.compound.upper()}"
#     )

#     # ax.set_ylim([-0.01, 1.02])
#     ax.set_xlabel("BV")
#     ax.set_ylabel(f"C/C$_0$ [{self.compound.upper()}]")
#     self.idx_labels = zip(
#         tmp.bv.to_list(), tmp.c_norm.to_list(), tmp.index.to_list()
#     )
#     for x, y, i in self.idx_labels:
#         # ax.annotate(f"{i}\n{y}", (x, y), fontsize="small", rotation=-45)
#         ax.annotate(i, (x, y), fontsize="small")
#     plt.tight_layout()
#     print(f"INDEXES: {tmp.index.to_list()}")
