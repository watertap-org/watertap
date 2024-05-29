import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
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
# from watertap.unit_models import IonExchangeClark
from .ion_exchange_clark import IonExchangeClark

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
    PropertyPackageError,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.ion_exchange_0D import IonExchange0D
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from analysisWaterTAP.analysis_scripts.PFAS_ix_analysis.pfas_ix_data import (
    get_pfas_ix_resin_data,
)
from copy import deepcopy

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# curve_file = __location__ + "/ix_pfas_curves_MASTER.csv"
# curve_data = pd.read_csv(curve_file)
# (
#     mw_dict,
#     resin_fg_dict,
#     resin_type_dict,
#     resin_dens_dict,
#     resin_diam_dict,
#     ground_resins
# ) = get_pfas_ix_resin_data()

class IXParmest:
    def __init__(
        self,
        curve_id,
        ix_model=IonExchangeClark,
        data_file=None,
        resin_data=dict(),
        compound_data=dict(),
        df_curve=None,
        bv=None,
        cb=None,
        c_norm=None,
        isotherm="freundlich",
        dont_calc_effluent=True,
        initial_guess_dict=dict(),
        initial_guess_theta_dict=dict(),
        set_bounds_dict=dict(),
        calc_from_constr_dict=dict(),
        autoscale_fixed=True,
        scale_from_value=None,
        fix_vars=dict(),
        c0_min_thresh=0.01,  # for determining keep_bvs + keep_cbs
        c0_max_thresh=0.9999,  # for determining keep_bvs + keep_cbs
        c_next_thresh=1.0,  # for determining keep_bvs + keep_cbs
        cb50_min_thresh=None,  # all cb > cb50_min_thresh are used to make linear regression estimate bv_50
        figsize=(7, 5),
        regenerant="single_use",
        just_plot_curve=False,
        max_zero=1e-3,
        min_one=0.9999,
        use_all_data=False,
        use_this_data=None,
        **kwargs,
    ):
        self.curve_id = curve_id
        if df_curve is None:
            curve_data = pd.read_csv(data_file)
            self.df_curve = curve_data[curve_data.curve_id == self.curve_id].copy()
        else:
            # option to provide separate DataFrame with breakthrough curve_id info
            self.df_curve = df_curve
        self.df_curve.sort_values("bv", inplace=True)
        self.figsize = figsize
        self.bv = bv  # need value for initial model build
        self.cb = cb  # need value for initial model build
        self.c_norm = c_norm
        self.all_figs = dict()  # dict for storing all figs
        self.isotherm = isotherm  # CONFIG option for IonExchange0D
        self.dont_calc_effluent = dont_calc_effluent  # boolean to indicate if ss effluent concentrations are deactivated
        self.regenerant = regenerant
        self.just_plot_curve = just_plot_curve
        self.use_all_data = use_all_data
        self.use_this_data = use_this_data
        self.resin_data = resin_data
        self.compound_data = compound_data
        self.ix_model = ix_model

        if self.use_all_data and self.use_this_data:
            # can't use all data and some data
            # default to using "this" data if both are True
            self.use_all_data = False

        if initial_guess_dict == dict():
            self.initial_guess_dict = {  # default initial guesses for regressed parameters
                "mass_transfer_coeff": 0.1,  # by default, bv_50 has no initial guess. Default initial guess for bv_50 is avg of BV values from curve_id.
                "freundlich_n": 1.1,
            }
        else:
            self.initial_guess_dict = deepcopy(initial_guess_dict)

        self.set_bounds_dict = deepcopy(set_bounds_dict)
        self.calc_from_constr_dict = deepcopy(calc_from_constr_dict)

        if initial_guess_theta_dict == dict():
            # if wanted to provide a different initial guess for theta
            self.initial_guess_theta_dict = deepcopy(self.initial_guess_dict)

        if len(self.initial_guess_dict) == 0 and len(self.initial_guess_theta_dict) == 0:
            raise ValueError("need to provide initial guess")
        # if autoscale_fixed is True, will automatically scale fixed variables with 1 / value(var)
        self.autoscale_fixed = autoscale_fixed
        # min/max thresholds for filtering data to be used in parmest
        self.c0_min_thresh = c0_min_thresh
        self.c0_max_thresh = c0_max_thresh
        self.c_next_thresh = c_next_thresh

        if self.c0_min_thresh <= 0:
            self.max_zero = max_zero
            max_zero_i = self.df_curve.query("c_norm <= 0").index.max()
            self.df_curve.at[max_zero_i, "c_norm"] = max_zero
            self.df_curve.at[max_zero_i, "cb"] = (
                max_zero * self.df_curve.at[max_zero_i, "c0"]
            )
            self.c0_min_thresh = 0.9 * max_zero

        if self.c0_max_thresh >= 1:
            self.min_one = min_one
            min_one_i = self.df_curve.query("c_norm >= 1").index.to_list()
            for i in min_one_i:
                self.df_curve.at[i, "c_norm"] = min_one
                self.df_curve.at[i, "cb"] = min_one * self.df_curve.at[i, "c0"]
            self.c0_max_thresh = min_one

        if cb50_min_thresh is None:
            # default for the min threshold for estimating BV50 is the same
            self.cb50_min_thresh = self.c0_min_thresh
        else:
            self.cb50_min_thresh = cb50_min_thresh

        if scale_from_value is None:
            # scale these variables with 1 / value(var)
            # default is to use the same variables provided in calc_from_constr_dict dict
            # since these vars would have an initial value that is (in theory) very close to the solved value
            self.scale_from_value = [k for k in calc_from_constr_dict.keys()]
            self.scale_from_value += ["bv"]
        else:
            self.scale_from_value = scale_from_value

        self.get_curve_conditions()
        self.get_model_config()

        # if fix_vars == dict():
        #     self.fix_vars = None
        # else:
        #     self.fix_vars = fix_vars

        # if self.just_plot_curve:
        #     self._just_plot_curve()
        #     return
        # self.build_it()
        self.rebuild()  # initial build of model
        # self.m0 = self.m.clone()  # intact initial build of model
        # self.R_sq = np.nan
        # self.build_results_dict()  # build dict for storing results

        # self.print_curve_conditions()

    def run_all_things(self, plot_things=False, save_things=False, overwrite=True):
        """Run all methods to complete parmest run
        1. plot_curve(): Plot extracted data
        2. estimate_bv50(): Estimate BV50 for parmest initial guess
        3. plot_estimate_bv50(): Plot the BV50 estimate
        4. test_initial_guess(): Run the model over the BV/Cb values using the initial guess values
        5. plot_initial_guess(): Plot initial guess
        6. run_parmest(): Run parmest
        7. test_theta(): Test resulting estimated parameters
        8. plot_theta(): Plot parmest results
        9. Optionally save figures, output, and results
        """
        self.estimate_bv50()
        # clear_output(wait=True)
        # self.test_initial_guess()
        # clear_output(wait=True)
        self.run_parmest()
        # clear_output(wait=True)
        self.test_theta()
        # clear_output(wait=True)
        if plot_things:
            self.plot_curve()
            self.plot_estimate_bv50()
            # self.plot_initial_guess()
            self.plot_theta()
        if save_things:
            if plot_things:
                self.save_figs(overwrite=overwrite)
            self.save_output(overwrite=overwrite)
            self.save_results(overwrite=overwrite)

    def rebuild(self):
        self.m = self.build_it()
        self.ix = self._get_ix_blk()
        self._set_bounds()
        self._fix_initial_guess()
        self._calc_from_constr()
        self.scale_it()
        assert degrees_of_freedom(self.m) == 0
        try:
            self.ix.initialize()
        except InitializationError:
            print("Initial build of model failed.")
            print_infeasible_constraints(self.m)

    def estimate_bv50(self):
        def linear_fit(bv, slope, b):
            return slope * bv + b

        def linear_inv(cb, slope, b):
            return (cb - b) / slope

        if not hasattr(self, "bv_pred_ig"):
            self.test_initial_guess()

        self.bv50_bvs = [
            bv
            for bv, cnorm in zip(self.keep_bvs, self.keep_cnorms)
            if cnorm >= self.cb50_min_thresh
        ]

        self.bv50_cbs = [
            cnorm for cnorm in self.keep_cnorms if cnorm >= self.cb50_min_thresh
        ]

        cb_50 = 0.5

        self.bv50_fit_params, self.bv50_fit_cov = curve_fit(
            linear_fit, self.bv50_bvs, self.bv50_cbs
        )

        self.bv50_slope, self.bv50_int = [*self.bv50_fit_params]

        self.bv50_est = linear_inv(cb_50, self.bv50_slope, self.bv50_int)

        self.bv50_orig = self.initial_guess_dict["bv_50"]
        self.bv_pred_ig_orig = deepcopy(self.bv_pred_ig)
        self.cb_pred_ig_orig = deepcopy(self.cb_pred_ig)

        self.bv50_corr_matrix = np.corrcoef(
            self.bv50_bvs, [self._linear_fit(bv) for bv in self.bv50_bvs]
        )
        corr = self.bv50_corr_matrix[0, 1]
        self.bv50_R_sq = corr**2
        if self.bv50_est > 100:
            self.initial_guess_dict["bv_50"] = self.bv50_est
        else:
            self.initial_guess_dict["bv_50"] = self.bv50_orig
        self.initial_guess_theta_dict = deepcopy(self.initial_guess_dict)
        # rebuild model with new initial guess

        self.rebuild()
        self.test_initial_guess()
        # self.plot_estimate_bv50()

    def plot_estimate_bv50(self):

        if self.bv50_est < 10:
            self.plot_initial_guess()
            return
        max_bv = max([self._linear_inv(1)] + self.all_bvs)
        cb_50 = 0.5

        title = (
            f"BV50 Estimate Results Compare - Curve {self.curve_id}:\n"
            + self.compound.swapcase()
            + " "
            + self.ref.replace("_", " ").title()
            + " "
            + self.resin.replace("_", " ").title()
            + f", EBCT = {self.ebct_min} min"
        )
        try:
            textstr = "\n".join(
                [
                    f"Guess",
                    f"  k$_T$: {self.initial_guess_dict['mass_transfer_coeff']:.2e}",
                    f"  n: {self.initial_guess_dict['freundlich_n']:.2f}",
                    f"  BV50 (new): {round(self.initial_guess_dict['bv_50'])}",
                    f"  BV50 (orig): {round(self.bv50_orig)}",
                ]
            )
        except:
            textstr = "\n".join(
                [
                    f"Guess",
                    f"  k$_T$: {self.m0.fs.ix.mass_transfer_coeff():.2e}",
                    f"  n: {self.initial_guess_dict['freundlich_n']:.2f}",
                    f"  BV50 (new): {round(self.initial_guess_dict['bv_50'])}",
                    f"  BV50 (orig): {round(self.bv50_orig)}",
                ]
            )
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.plot(
            self.keep_bvs,
            self.keep_cnorms,
            marker=".",
            color="red",
            alpha=0.2,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cnorms,
            marker=".",
            color="purple",
            alpha=0.2,
            label="All Data",
        )
        ax.scatter(
            self.bv_pred_ig_orig,
            self.cb_pred_ig_orig,
            marker="o",
            color="green",
            label="Predicted (Orig)",
        )
        ax.scatter(
            self.bv_pred_ig,
            self.cb_pred_ig,
            marker="o",
            color="blue",
            label="Predicted (New)",
        )
        ax.scatter(
            self.bv_pred_ig_fake,
            self.cb_pred_ig_fake,
            marker=".",
            color="blue",
            alpha=0.1,
            # label="Predicted",
        )
        if len(self.bv_fail_ig) != 0:
            ax.scatter(
                self.bv_fail_ig,
                self.cb_fail_ig,
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

        xs = np.linspace(min(self.bv50_bvs), max_bv, len(self.bv50_bvs))
        ys = [self._linear_fit(x) for x in xs]

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
            [self.bv50_orig],
            [cb_50],
            marker="*",
            color="springgreen",
            label="BV50 Guess (Orig)",
        )
        ax.scatter(
            [self.bv50_est],
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
        ylabe = self.compound.swapcase() + f" [{self.conc_units_str}]"
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.01, 1.02])
        ax.set_xlim([0, max_bv])
        ax.legend(loc="lower right")
        plt.tight_layout()

        self.bv50_fig = fig
        self.all_figs["bv_50_est"] = {"fig": fig, "ax": ax}

    def filter_data(self):

        """
        Method to filter extracted breakthrough curve
        Loops through provided BVs and Cbs and keeps (bv, cb) pairs that meet all conditions:
        1. Current cb is greater than the last kept cb (last_cb)
        2. Current cb is greater than the minimum threshold
        3. Current cb is less than the maximum threshold
        4. Associated BV is not zero
        5. (if the current cb is not the last) The next cb (i + 1) is some threshold greater than the current cb
        """

        df = self.df_curve.copy()
        last_cb = -0.01
        self.keep_cbs = []
        self.keep_bvs = []
        self.excl_cbs = []
        self.excl_bvs = []
        self.all_bvs = df.bv.to_list()
        self.all_cbs = df.cb.to_list()
        self.all_cnorms = df.c_norm.to_list()

        if self.use_all_data:
            print("use all data")
            tmp = df[df.c_norm != 0].copy()
            self.keep_bvs = tmp.bv.to_list()
            self.keep_cbs = tmp.cb.to_list()
            assert 0 not in self.keep_cbs
        elif self.use_this_data is not None:
            print("use this data")
            tmp = df.reset_index(drop=True)
            tmp = tmp.loc[self.use_this_data].copy()
            self.keep_bvs = tmp.bv.to_list()
            self.keep_cbs = tmp.cb.to_list()
            assert 0 not in self.keep_cbs

        else:
            print("default filtering")

            for i, cb in enumerate(self.all_cbs):
                if i != len(self.all_cbs) - 1:
                    if (
                        # cb > last_cb
                        # and cb < self.c0
                        cb < self.c0
                        and cb >= self.c0_min_thresh * self.c0
                        and cb <= self.c0_max_thresh * self.c0
                        and self.all_bvs[i] > 0
                        and self.all_cbs[i + 1] > cb * self.c_next_thresh
                    ):
                        last_cb = cb
                        self.keep_cbs.append(cb)
                        self.keep_bvs.append(self.all_bvs[i])
                    else:
                        self.excl_cbs.append(cb)
                        self.excl_bvs.append(self.all_bvs[i])
                else:
                    if (
                        # cb > last_cb
                        # and cb < self.c0
                        cb < self.c0
                        and cb >= self.c0_min_thresh * self.c0
                        and cb <= self.c0_max_thresh * self.c0
                        and self.all_bvs[i] > 0
                    ):
                        last_cb = cb
                        self.keep_cbs.append(cb)
                        self.keep_bvs.append(self.all_bvs[i])
                    else:
                        self.excl_cbs.append(cb)
                        self.excl_bvs.append(self.all_bvs[i])

        self.keep_cnorms = [cb / self.c0 for cb in self.keep_cbs]
        self.excl_cnorms = [cb / self.c0 for cb in self.excl_cbs]

    def plot_curve(self):

        fig, ax = plt.subplots(figsize=self.figsize)
        ax.scatter(
            self.excl_bvs, self.excl_cbs, marker="x", color="k", label="Excluded data"
        )
        # ax.scatter(keep_bv, keep_cb, color="r", marker=".", label="Fitted data")
        ax.plot(
            self.keep_bvs,
            self.keep_cnorms,
            color="red",
            marker=".",
            alpha=0.25,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cnorms,
            color="purple",
            marker=".",
            alpha=0.25,
            label="All Data",
        )
        ax.plot(
            [0, max(self.all_bvs)],
            [1, 1],
            linestyle="-.",
            alpha=0.25,
            label="Influent",
        )

        title = (
            f"Data Used for Curve {self.curve_id}\n"
            + self.compound.swapcase()
            + ", "
            + self.ref.replace("_", " ").title()
            + ", "
            + self.resin.replace("_", " ").title()
            + f", EBCT = {self.ebct_min} min"
        )
        ax.set_xlabel("BV")
        # ylabe = self.compound.swapcase()
        ylabe = "C/C$_0$, " + self.compound.swapcase()
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.05, 1.02])
        ax.set_xlim([-5, max(self.all_bvs) * 1.05])
        ax.legend()
        plt.tight_layout()

        self.curve_fig = fig

        self.all_figs["raw_curve"] = {"fig": fig, "ax": ax}

    def test_initial_guess(self, min_cnorm_ig=0.05, max_cnorm_ig=0.95, num_pts=20):
        """
        Loops through the kept (bv, cb) pairs and runs the model
        with those values using initial_guess_dict
        """

        print("TEST INITIAL GUESS\n")
        ix = self.ix
        self._fix_initial_guess()
        self._calc_from_constr()

        self.test_cnorms0 = np.linspace(max_cnorm_ig, min_cnorm_ig, num_pts)
        self.bv_pred_ig_fake = []  # bv_pred_ig = Predicted BVs using initial guess
        # self.bv_fail_ig_fake = []
        self.cb_pred_ig_fake = []  # cb_pred_ig = Predicted Cbs using initial guess
        self.cb_fail_ig_fake = []

        self.bv_pred_ig = []  # bv_pred_ig = Predicted BVs using initial guess
        self.bv_fail_ig = []
        self.cb_pred_ig = []  # cb_pred_ig = Predicted Cbs using initial guess
        self.cb_fail_ig = []
        self.flag = "test_initial_guess_fake"

        m = self.m.clone()

        ix = m.fs.ix

        for cnorm in self.test_cnorms0:
            ix.c_norm.fix(cnorm)
            try:
                ix.initialize()
            except:
                # print_infeasible_constraints(ix)
                pass
            self.solve_it(m=m)
            if self.tc == "optimal":
                self.bv_pred_ig_fake.append(ix.bv())
                self.cb_pred_ig_fake.append(cnorm)
            elif self.tc != "optimal":
                # self.bv_fail_ig_fake.append(bv)
                self.cb_fail_ig_fake.append(cnorm)
            self.results_dict_append(blk=ix)

        self.flag = "test_initial_guess"
        m = self.m.clone()
        ix = m.fs.ix
        for (bv, cnorm) in zip(self.keep_bvs[::-1], self.keep_cnorms[::-1]):
            ix.bv.set_value(bv)
            ix.c_norm.fix(cnorm)

            try:
                ix.initialize()
            except:
                pass
            self.solve_it(m=m)
            if self.tc == "optimal":
                self.bv_pred_ig.append(ix.bv())
                self.cb_pred_ig.append(cnorm)
            elif self.tc != "optimal":
                self.bv_fail_ig.append(bv)
                self.cb_fail_ig.append(cnorm)
            self.results_dict_append(blk=ix)

    def plot_initial_guess(self):

        fig, ax = plt.subplots(figsize=self.figsize)
        ax.scatter(
            self.bv_pred_ig,
            self.cb_pred_ig,
            marker="o",
            color="green",
            label="Predicted",
        )
        ax.scatter(
            self.bv_pred_ig_fake,
            self.cb_pred_ig_fake,
            marker=".",
            color="blue",
            alpha=0.1,
            # label="Predicted",
        )
        if len(self.bv_fail_ig) != 0:
            ax.scatter(
                self.bv_fail_ig,
                self.cb_fail_ig,
                marker="x",
                color="black",
                label="Infeasible",
            )

        ax.plot(
            self.keep_bvs,
            self.keep_cnorms,
            marker=".",
            color="red",
            alpha=0.25,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cnorms,
            marker=".",
            color="purple",
            alpha=0.25,
            label="All Data",
        )
        ax.plot(
            [0, max(self.bv_pred_ig + self.bv_pred_ig_fake + self.keep_bvs)],
            [1, 1],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
        )
        ax.plot(
            [0, max(self.bv_pred_ig + self.bv_pred_ig_fake + self.keep_bvs)],
            [0.5, 0.5],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
        )

        ax.scatter(
            [self.initial_guess_dict["bv_50"]],
            [0.5],
            marker="*",
            color="springgreen",
            # label="BV50 Guess",
        )

        title = (
            f"Initial Guess - Curve {self.curve_id}:\n"
            + self.compound.swapcase()
            + " "
            + self.ref.replace("_", " ").title()
            + " "
            + self.resin.replace("_", " ").title()
            + f" EBCT = {self.ebct_min} min"
        )
        ax.set_xlabel("BV")
        # ylabe = self.compound.swapcase() + f" [{self.conc_units_str}]"
        ylabe = "C/C$_0$" + self.compound.swapcase()
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.01, 1.02])

        textstr = "\n".join(
            [
                f"Guess",
                f"  k$_T$: {self.initial_guess_dict['mass_transfer_coeff']:.2e}",
                f"  $n$: {self.initial_guess_dict['freundlich_n']:.2f}",
                f"  BV50: {round(self.initial_guess_dict['bv_50'])}",
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

        self.ig_fig = fig

        self.all_figs["initial_guess"] = {"fig": fig, "ax": ax}

    def get_curve_conditions(self):
        """
        Gets conditions for the breakthrough curve_id to use in model
        and other metadata
        """


        df = self.df_curve.copy()

        # TODO: insert check for columns that should all have same value

        self.ref = df.ref.iloc[0]
        self.c0 = df.c0.iloc[0]

        self.compound = df.compound.iloc[0]
        self.flow_in = df.flow_in.iloc[0] * pyunits.m**3 / pyunits.s # m3/s
        self.loading_rate = df.loading_rate.iloc[0] * pyunits.m / pyunits.s # m/s
        self.bed_depth = df.bed_depth.iloc[0] * pyunits.m
        self.bed_diam = df.bed_diam.iloc[0] * pyunits.m
        self.bed_vol_tot = df.bed_vol.iloc[0] * pyunits.m**3
        self.resin = df.resin.iloc[0]
        # self.resin_type = resin_type_dict[self.resin]
        # self.resin_func_group = resin_fg_dict[self.resin]
        self.ebct_min = df.ebct.iloc[0] # minutes
        self.ebct = pyunits.convert(
            self.ebct_min * pyunits.min, to_units=pyunits.second
        )
        self.conc_units_str = df.conc_units.iloc[0] # conc_mass_phase_comp units
        numer = getattr(pyunits, self.conc_units_str.split("/")[0])
        denom = getattr(pyunits, self.conc_units_str.split("/")[1])
        self.conc_units = numer / denom
        # self.charge = -1
        self.target_component = self.compound_data["name"]
        self.charge = self.compound_data["charge"]
        self.mw = self.compound_data["mw_comp"] * pyunits.kg / pyunits.mol
        self.diffusivity = self.compound_data["diffusivity"]
        # try:
        #     self.mw = mw_dict[self.target_component] * pyunits.kg / pyunits.mol
        # except:
        #     self.mw = mw_dict[self.compound] * pyunits.kg / pyunits.mol
        self.filter_data()

        if self.just_plot_curve:
            return

        try:
            if ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2) > 0:
                self.expr_sf = 1 / ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2)
        except:
            self.expr_sf = 1e-9

        # self.resin_density, self.resin_diam = self._get_resin_deets(self.resin)
        self.resin_density = self.resin_data["density"]
        self.resin_diam = self.resin_data["diameter"]
        self.bed_porosity = self.resin_data["porosity"]

    def get_model_config(self):

        mw_water = 0.018 * pyunits.kg / pyunits.mol
        rho = 1000 * pyunits.kg / pyunits.m**3

        self.ion_props = {
            "solute_list": [self.target_component],
            "diffusivity_data": {("Liq", self.target_component): self.diffusivity},
            "mw_data": {"H2O": mw_water, self.target_component: self.mw},
            "charge": {self.target_component: self.charge},
        }

        self.ix_config = {
            "target_component": self.target_component,
            "regenerant": self.regenerant,
        }

        if np.isnan(self.loading_rate):
            self.loading_rate = self.bed_depth / value(self.ebct)

        self.c0_mol_flow = pyunits.convert(
            (self.c0 * self.conc_units)
            / self.mw
            * (self.flow_in * pyunits.m**3 / pyunits.s),
            to_units=pyunits.mol / pyunits.s,
        )
        self.c0_mol_flow_sf = 1 / value(self.c0_mol_flow)
        self.water_mol_flow = pyunits.convert(
            ((self.flow_in * pyunits.m**3 / pyunits.s) * rho) / mw_water,
            to_units=pyunits.mol / pyunits.s,
        )
        self.water_mol_flow_sf = 1 / self.water_mol_flow()

        if self.cb is None:
            self.cb = self.keep_cbs[-1] * self.conc_units
        if self.c_norm is None:
            self.c_norm = self.keep_cnorms[-1]
        if self.bv is None:
            self.bv = self.keep_bvs[-1]

        if "bv_50" not in self.initial_guess_dict.keys():
            self.initial_guess_dict["bv_50"] = np.mean(self.keep_bvs)
        if "bv_50" not in self.initial_guess_theta_dict.keys():
            self.initial_guess_theta_dict["bv_50"] = np.mean(self.keep_bvs)

    def build_it(self):
        """
        Method used to build IX model
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**self.ion_props)

        self.ix_config["property_package"] = m.fs.properties

        m.fs.ix = ix = self.ix_model(**self.ix_config)

        # remove lower bounds for data fitting
        # small-scale experimental systems will likely fall outside these bounds
        ix.bed_depth.setlb(0)
        ix.ebct.setlb(0)
        ix.bv.setlb(0)
        ix.bv.setub(None)

        self._set_bounds(m=m)

        pf = ix.process_flow

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", self.c0_mol_flow_sf, index=("Liq", self.target_component)
        )

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", self.water_mol_flow_sf, index=("Liq", "H2O")
        )

        c_in = pyunits.convert(
            self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3
        )
        c_b = pyunits.convert(self.cb, to_units=pyunits.kg / pyunits.m**3)

        if self.c_norm is None:
            self.c_norm = value(c_b / c_in)

        pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(self.water_mol_flow)
        pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].fix(
            self.c0_mol_flow
        )
        pf.properties_in[0].pressure.fix(101325)
        pf.properties_in[0].temperature.fix(298)
        pf.properties_in[0].flow_vol_phase[...]
        pf.properties_out[0].flow_vol_phase[...]

        # Parameters set to force bed_depth == col_height
        ix.underdrain_h.set_value(0)
        ix.distributor_h.set_value(0)
        if self.regenerant != "single_use":
            ix.bed_expansion_frac_A.set_value(0)
            ix.bw_rate.set_value(0)
        ix.number_columns.fix(1)

        # ix.c_breakthru[self.target_component].fix(value(c_b))
        ix.c_norm[self.target_component].fix(self.c_norm)
        ix.resin_density.fix(self.resin_density)
        ix.resin_diam.fix(self.resin_diam)
        ix.bed_porosity.fix(self.bed_porosity)
        ix.bed_depth.fix(self.bed_depth)
        ix.loading_rate.fix(self.loading_rate)


        ix.bv.set_value(
            self.bv
        )  # since we want to predict BV, we use .set_value() rather than .fix()

        # if self.dont_calc_effluent:
            # Since we are fitting data to C/C0, we don't care about the steady-state effluent concentration, so deactivate them
        ix.eq_tb_traps.deactivate()
        ix.eq_c_traps.deactivate()
        ix.eq_traps.deactivate()
        ix.eq_c_norm_avg.deactivate()
        for _, v in ix.tb_traps.items():
            v.fix()
        for _, v in ix.c_traps.items():
            v.fix()
        for _, v in ix.traps.items():
            v.fix()
        ix.c_norm_avg.fix()
        ix.eq_mass_transfer_target_fr.deactivate()
        pf.mass_transfer_term[0, "Liq", self.target_component].fix()
        ix.eq_service_flow_rate.deactivate()

        # if isinstance(self.fix_vars, dict):
            # If fix_vars is provided
            # for k, v in self.fix_vars.items():
            #     getattr(m.fs.ix, k).fix(v)

        # ix.eq_c_breakthru.deactivate()

        return m

    def scale_it(self, m=None):
        """
        Scale the model
        """
        if m is None:
            m = self.m
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)

        target_component = ix.config.target_component
        prop_out = ix.process_flow.properties_out[0]
        set_scaling_factor(
            prop_out.flow_mol_phase_comp["Liq", target_component], self.c0_mol_flow_sf
        )
        set_scaling_factor(
            prop_out.flow_mol_phase_comp["Liq", "H2O"], self.water_mol_flow_sf
        )

        calculate_scaling_factors(m)

        if self.autoscale_fixed:
            for v in self._get_comp_list(ix, skip_list=["traps"]):
                # print(v)
                ixv = getattr(ix, v)
                if ixv.is_indexed():
                    idx = [*ixv.index_set()]
                    for i in idx:
                        if ixv[i].is_fixed():
                            if ixv[i]() == 0:
                                continue
                            # print(i, ixv.name, ixv[i](), 1 / ixv[i]())
                            set_scaling_factor(ixv[i], 1 / ixv[i]())
                else:
                    if ixv.is_fixed():
                        if ixv() == 0:
                            continue
                        set_scaling_factor(ixv, 1 / ixv())
                        # print(v, "fixed")

        if isinstance(self.scale_from_value, list):
            for v in self.scale_from_value:
                # print(v)
                ixv = getattr(ix, v)
                if ixv.is_indexed():
                    idx = [*ixv.index_set()]
                    for i in idx:
                        if ixv[i]() == 0:
                            continue
                        set_scaling_factor(ixv[i], 1 / ixv[i]())
                else:
                    if ixv() == 0:
                        continue
                    set_scaling_factor(ixv, 1 / ixv())

        constraint_scaling_transform(ix.eq_clark[target_component], 1e-2)

    def solve_it(self, m=None, solver="watertap", optarg=dict(), tee=False):

        """
        Solve the model
        Defaults to using watertap solver
        """

        if m is None:
            m = self.m

        if solver != "watertap":
            self.solver = SolverFactory("ipopt")
        else:
            self.solver = get_solver()
            for k, v in optarg.items():
                self.solver.options[k] = v
        try:
            # Should never break a model run if an error is thrown
            self.results = self.solver.solve(m, tee=tee)
            self.tc = self.results.solver.termination_condition
        except:
            self.tc = "fail"

    def run_parmest(self):
        """
        Run parmest
        """

        self.theta_names = [
            "fs.ix.mass_transfer_coeff",
            "fs.ix.freundlich_n",
            "fs.ix.bv_50",
        ]
        self.theta_input = dict(
            zip(self.theta_names, self.initial_guess_theta_dict.values())
        )
        self.theta_values = pd.DataFrame(
            data=[[v for v in self.initial_guess_theta_dict.values()]],
            columns=self.theta_names,
        )
        self.df_parmest = self.df_curve[self.df_curve.cb.isin(self.keep_cbs)]

        def parmest_regression(data):

            """
            Build function that is passed to parmest
            """

            cb = pyunits.convert(
                data.cb.to_list()[0] * self.conc_units,
                to_units=pyunits.kg / pyunits.m**3,
            )()
            cb_p = data.cb.to_list()[0] * self.conc_units
            c0_p = self.c0 * self.conc_units
            cx = value(cb_p / c0_p)
            bv = data.bv.to_list()[0]

            print(
                f"\nPARMEST FOR:\n\t"
                f"CURVE = {self.curve_id}\n\t"
                f"REF = {self.ref}\n\t"
                f"COMPOUND = {self.compound}\n\t"
                f"BV = {bv}\n\t"
                f"C0 = {c0_p}\n\t"
                f"CB = {cb_p}\n\t"
                f"CNORM = {cx}\n\t"
            )

            m_parmest = self.m.clone()
            m_parmest.fs.ix.c_norm.fix(cx)
            m_parmest.fs.ix.bv.set_value(bv)
            self._calc_from_constr(m=m_parmest)
            self.scale_it(m=m_parmest)
            try:
                m_parmest.fs.ix.initialize()
            except:
                pass
            # Stores the instantiated model that is passed to parmest for debugging outside of the class
            self.m_parmest = m_parmest.clone()
            return m_parmest

        def SSE(m, data):
            """
            Objective function for parmest to minimize
            We want the SSE of BV observed/predicted to be small
            """
            ix = self._get_ix_blk(m=m)
            expr = (float(data.bv) - ix.bv) ** 2
            return expr * self.expr_sf  # expr_sf is the scaling factor for the SSE

        self.pest = parmest.Estimator(
            parmest_regression,
            self.df_parmest,
            self.theta_names,
            SSE,
            tee=False,
            diagnostic_mode=False,
            solver_options={"max_iter": 10000},
        )

        self.pest.objective_at_theta(
            theta_values=self.theta_values, initialize_parmest_model=True
        )

        self.obj, self.theta = self.pest.theta_est()

        self.theta_dict = dict(zip(self.initial_guess_theta_dict.keys(), [*self.theta]))

        if hasattr(self, "results_dict"):
            for k, v in self.theta_dict.items():
                self.results_dict[f"{k}_theta"] = [
                    v for _ in range(len(self.results_dict["curve_id"]))
                ]
            self.results_dict["obj"] = [
                self.obj for _ in range(len(self.results_dict["curve_id"]))
            ]

    def test_theta(
        self,
        min_test_cbs=None,
        max_test_cbs=0.99,
        min_test_cnorms=0.05,
        max_test_cnorms=0.95,
        min_test_bvs=None,
        max_test_bvs=None,
        num_pts=20,
    ):
        """
        Test parmest param estimation
        First tested against regularly spaced BV and Cb values
        Then tested against extracted data
        """
        print("\n\nTEST THETA\n\n")

        if min_test_cbs is None:
            self.min_test_cbs = min(self.keep_cnorms)
        else:
            self.min_test_cbs = min_test_cbs
        if min_test_cnorms is None:
            self.min_test_cnorms = min(self.keep_cnorms)
        else:
            self.min_test_cnorms = min_test_cnorms
        if min_test_bvs is None:
            self.min_test_bvs = min(self.keep_bvs)
        else:
            self.min_test_bvs = min_test_bvs
        if max_test_bvs is None:
            self.max_test_bvs = max(self.all_bvs)
        else:
            self.max_test_bvs = max_test_bvs

        self.max_test_cbs = max_test_cbs
        self.max_test_cnorms = max_test_cnorms
        self.num_pts = num_pts

        self.test_cbs = [
            self.c0 * x
            for x in np.linspace(self.max_test_cbs, self.min_test_cbs, self.num_pts)
        ]
        self.test_cnorms = np.linspace(
            self.max_test_cnorms, self.min_test_cnorms, self.num_pts
        )
        self.test_bvs = np.linspace(self.max_test_bvs, self.min_test_bvs, self.num_pts)

        # Comparing against fake data
        self.bv_pred_theta_test = (
            []
        )  # bv_pred_theta_test = Predicted BV using theta and regularly spaced test data
        self.bv_fail_theta_test = []

        self.cnorm_pred_theta_test = (
            []
        )  # cb_pred_theta_test = Predicted Cb using theta and regularly spaced test data
        self.cnorm_fail_theta_test = []

        self.rebuild()
        self.m_theta = self.m.clone()
        self.ix = ix = self._get_ix_blk(m=self.m_theta)
        for k, v in self.theta_dict.items():
            ixv = getattr(ix, k)
            ixv.fix(v)

        self._calc_from_constr(m=self.m_theta)
        self.scale_it(m=self.m_theta)

        try:
            ix.initialize()
        except:
            # print_infeasible_constraints(ix)
            pass

        assert degrees_of_freedom(self.m_theta) == 0

        self.flag = "test_theta_fake"
        m = self.m_theta.clone()
        ix = m.fs.ix
        for bv, cnorm in zip(self.test_bvs, self.test_cnorms):

            ix.c_norm[self.target_component].fix(cnorm)
            ix.bv.set_value(bv)

            try:
                self._calc_from_constr(m=m)
                ix.initialize()
                ix.initialize()
            except:
                # print_infeasible_constraints(m)
                pass

            assert degrees_of_freedom(m) == 0
            self.solve_it(m=m)
            if self.tc != "optimal":
                self.bv_fail_theta_test.append(bv)
                self.cnorm_fail_theta_test.append(cnorm)
            elif self.tc == "optimal":
                self.bv_pred_theta_test.append(ix.bv())
                self.cnorm_pred_theta_test.append(cnorm)
            self.results_dict_append(blk=ix)

        self.bv_pred_theta = []
        self.bv_fail_theta = []
        self.keep_bv_theta = []
        self.cnorm_pred_theta = []
        self.cnorm_fail_theta = []

        self.flag = "test_theta_real"
        m = self.m_theta.clone()
        ix = m.fs.ix
        for bv, cnorm in zip(self.keep_bvs[::-1], self.keep_cnorms[::-1]):
            if cnorm == 0:
                cnorm = 1e-3

            ix.c_norm[self.target_component].fix(cnorm)
            ix.bv.set_value(bv)
            self._calc_from_constr(m=m)
            # self.scale_it(m=self.m_theta)
            try:
                ix.initialize()
                ix.initialize()
            except:
                # print_infeasible_constraints(m)
                pass

            assert degrees_of_freedom(m) == 0
            self.solve_it(m=m)
            if self.tc != "optimal":
                self.bv_fail_theta.append(bv)
                self.cnorm_fail_theta.append(cnorm)
            elif self.tc == "optimal":
                self.bv_pred_theta.append(ix.bv())
                self.keep_bv_theta.append(bv)
                self.cnorm_pred_theta.append(cnorm)
            self.results_dict_append(blk=ix)

        self.corr_matrix = np.corrcoef(self.keep_bv_theta, self.bv_pred_theta)
        corr = self.corr_matrix[0, 1]
        self.R_sq = corr**2
        # clear_output(wait=True)

    def plot_theta(self):

        print("\n\nPLOT THETA\n\n")

        if not all(hasattr(self, attr) for attr in ["test_cnorms", "test_bvs", "R_sq"]):
            self.test_theta()

        textstr = "\n".join(
            [
                f"Theta (R$^2$: {self.R_sq:.3f})",
                f"  k$_T$: {self.theta_dict['mass_transfer_coeff']:.2e}",
                f"  n: {self.theta_dict['freundlich_n']:.2f}",
                f"  BV50: {round(self.theta_dict['bv_50'])}",
            ]
        )
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)

        fig1, ax1 = plt.subplots(figsize=self.figsize)

        ax1.scatter(
            self.bv_pred_theta_test,
            self.cnorm_pred_theta_test,
            color="blue",
            marker=".",
            label="Predicted",
        )
        ax1.plot(
            self.test_bvs,
            self.test_cnorms,
            color="red",
            linestyle="-.",
            alpha=0.1,
            label="Test Data",
        )

        if len(self.bv_fail_theta_test) != 0:
            ax1.scatter(
                self.bv_fail_theta_test,
                self.cnorm_fail_theta_test,
                marker="x",
                color="black",
                label="Infeasible",
            )

        fig2, ax2 = plt.subplots(figsize=self.figsize)
        ax2.plot(
            self.bv_pred_theta_test,
            self.cnorm_pred_theta_test,
            color="blue",
            marker=".",
            alpha=0.1,
            # label="Predicted",
        )
        ax2.scatter(
            self.bv_pred_theta,
            self.cnorm_pred_theta,
            marker=".",
            color="green",
            label="Predicted",
        )
        ax2.plot(
            self.keep_bvs,
            self.keep_cnorms,
            marker=".",
            color="red",
            alpha=0.2,
            label="Parmest Data",
        )

        ax2.plot(
            self.all_bvs,
            self.all_cnorms,
            marker=".",
            color="purple",
            alpha=0.2,
            label="All Data",
        )

        if len(self.bv_fail_theta) != 0:
            ax2.scatter(
                self.bv_fail_theta,
                self.cnorm_fail_theta,
                marker="x",
                color="black",
                label="Infeasible",
            )

        # ylabe = self.compound.swapcase() + f" [{self.conc_units_str}]"
        ylabe = "C/C$_0$, " + self.compound.swapcase()
        title = (
            f"Curve {self.curve_id}:\n"
            + self.compound.swapcase()
            + ", "
            + self.ref.replace("_", " ").title()
            + ", "
            + self.resin.replace("_", " ").title()
            + ", "
            + f"{self.ebct_min} min EBCT"
        )
        for ax in [ax1, ax2]:
            ax.plot(
                [
                    0,
                    max(
                        self.keep_bvs
                        + self.bv_pred_theta_test
                        + self.bv_pred_theta
                        + self.all_bvs
                    ),
                ],
                [
                    # self.c0
                    1.0
                    for _ in [
                        0,
                        max(
                            self.keep_bvs
                            + self.bv_pred_theta_test
                            + self.bv_pred_theta
                            + self.all_bvs
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
        self.all_figs["theta_test"] = {"fig": fig1, "ax": ax1}
        self.all_figs["theta_real"] = {"fig": fig2, "ax": ax2}

        fig3, ax3 = plt.subplots(figsize=self.figsize)

        ax3.scatter(
            self.keep_bv_theta,
            self.bv_pred_theta,
            marker=".",
            color="green",
            label="Predicted",
        )
        ax3.plot(self.keep_bvs, self.keep_bvs, linestyle=":", color="red", alpha=0.25)
        ax3.set_xlabel("Actual BV")
        ax3.set_ylabel("Predicted BV")
        ax3.set_title(f"Curve {self.curve_id} Parity")
        plt.tight_layout()
        self.all_figs["theta_parity"] = {"fig": fig3, "ax": ax3}

        fig4, ax4 = plt.subplots(figsize=self.figsize)
        self.bv_error = [
            pred - actual
            for (pred, actual) in zip(self.bv_pred_theta, self.keep_bv_theta)
        ]
        self.bv_rel_error = [
            (pred - actual) / actual
            for (pred, actual) in zip(self.bv_pred_theta, self.keep_bv_theta)
        ]

        ax4.scatter(self.keep_bv_theta, self.bv_rel_error, color="black")
        ax4.set_xlabel("Actual BV")
        ax4.set_ylabel("Relative Error")
        ax4.set_title(f"Curve {self.curve_id} Relative Error")
        plt.tight_layout()

        self.all_figs["theta_relerror"] = {"fig": fig4, "ax": ax4}
        fig4.show()

    def print_curve_conditions(self):
        # clear_output(wait=False)
        print(f"\nCURVE: {self.curve_id}")
        print(f"  Reference = {self.ref.replace('_', ' ').title()}")
        print(f"  Compound = {self.compound.upper()}")
        print(f"  Resin = {self.resin.upper()}")
        print(f"  EBCT = {self.ebct_min} min \n")

    def build_results_dict(
        self,
        split_str="fs.ix.",
        components=[Var, Expression],
        blk=None,
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
        self.results_dict = {}
        self.v_name_map = {}
        self.curve_deets = [
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
            blk = self._get_ix_blk()
        elif isinstance(blk, ConcreteModel):
            blk = self._get_ix_blk(m=blk)

        assert isinstance(blk, IonExchange0D)

        for deet in self.curve_deets:
            self.results_dict[deet] = []

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
                            self.results_dict[v[i].name] = []
                            self.v_name_map[v[i].name] = v_name_i
                    else:
                        self.results_dict[v.name] = []
                        self.v_name_map[v.name] = v_name
                    continue

                if any(True if s in v.name else False for s in skips):
                    # print(f"SKIP {v.name}")
                    continue
                if v.is_indexed():
                    idx = [*v._index_set]
                    for i in idx:
                        # print(f'KEEP INDEXED {repr(c), i, v.name}')
                        v_name_i = v_name + f"_{i}"
                        self.results_dict[v[i].name] = []
                        self.v_name_map[v[i].name] = v_name_i
                else:
                    self.results_dict[v.name] = []
                    self.v_name_map[v.name] = v_name

        for k in self.initial_guess_dict.keys():
            self.results_dict[f"{k}_ig"] = []

        for k in self.initial_guess_theta_dict.keys():
            self.results_dict[f"{k}_ig_theta"] = []

        for k in self.calc_from_constr_dict.keys():
            self.results_dict[f"{k}_cfc"] = []

        self.results_dict["flag"] = []

        if check_keys:
            for k in self.results_dict.keys():
                print(k)

    def results_dict_append(
        self,
        blk=None,
        components=[Var, Expression],
        tmp_results_dict=None,
    ):

        if blk is None:
            blk = self.ix
        elif isinstance(blk, ConcreteModel):
            blk = self._get_ix_blk(m=blk)
        else:
            pass

        assert isinstance(blk, IonExchange0D)

        for deet in self.curve_deets:
            # print(deet)
            self.results_dict[deet].append(getattr(self, deet))

        if tmp_results_dict is None:
            for c in components:
                # if c is Var:
                for v in blk.component_objects(c):
                    if v.is_indexed():
                        idx = [*v._index_set]
                        for i in idx:
                            if v[i].name in self.results_dict.keys():
                                self.results_dict[v[i].name].append(v[i]())
                    else:
                        if v.name in self.results_dict.keys():
                            # print(v.name, v())
                            self.results_dict[v.name].append(v())

            for k, v in self.initial_guess_dict.items():
                self.results_dict[f"{k}_ig"].append(v)

            for k, v in self.initial_guess_theta_dict.items():
                self.results_dict[f"{k}_ig_theta"].append(v)

            for k, v in self.calc_from_constr_dict.items():
                self.results_dict[f"{k}_cfc"].append(v)

            if hasattr(self, "theta_dict"):
                for k, v in self.theta_dict.items():
                    self.results_dict[f"{k}_theta"].append(v)
                self.results_dict["obj"].append(self.obj)

            self.results_dict["flag"].append(self.flag)

        else:
            for k, v in tmp_results_dict.items():
                self.results_dict[k].append(v)

    def save_results(self, overwrite=False, results_filename=None):

        self.results_dict_save = deepcopy(self.results_dict)

        for old_key, new_key in self.v_name_map.items():
            self.results_dict_save[new_key] = self.results_dict_save.pop(old_key)

        if results_filename is None:
            self.results_filename = f"results/curve_id{self.curve_id}_results.csv"
        else:
            self.results_filename = results_filename

        self.df_results = pd.DataFrame.from_dict(self.results_dict_save)

        if overwrite:
            self.df_results.to_csv(self.results_filename, index=False)
        elif not overwrite:
            x = self.results_filename.split(".csv")[0]
            append = 2
            while os.path.exists(self.results_filename):
                self.results_filename = x + f"_{append}.csv"
                append += 1
            self.df_results.to_csv(self.results_filename, index=False)

    def save_figs(self, overwrite=False, extension=None):

        # fig_file_base = f"figs/curve_id{self.curve_id}_"

        for figname, d in self.all_figs.items():
            # fig_file_base = f"figs/{figname}/curve_id{self.curve_id}_"
            fig = d["fig"]
            fig_file = f"figs/{figname}/curve_id{self.curve_id}_{figname}.png"
            if overwrite:
                fig.savefig(fig_file, bbox_inches="tight")
            elif not overwrite:
                x = fig_file.split(".png")[0]
                append = 2
                while os.path.exists(fig_file):
                    fig_file = x + f"_{append}.png"
                    append += 1
                fig.savefig(fig_file, bbox_inches="tight")

    def save_output(self, overwrite=False):

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

        data_file_base = f"output/curve_id{self.curve_id}_output.csv"
        theta_file_base = f"theta/curve_id{self.curve_id}_theta.csv"

        if hasattr(self, "theta_dict"):
            tmp_dict = dict()
            tmp_dict[self.curve_id] = self.theta_dict
            df_theta = pd.DataFrame.from_dict(tmp_dict).T
            df_theta["loading_rate"] = self.loading_rate
            df_theta["curve_id"] = self.curve_id
            df_theta["c0"] = pyunits.convert(
                self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3
            )()
            df_theta.to_csv(theta_file_base, index=False)

        tmps = []
        for data in sorted(datas):
            if hasattr(self, data):
                tmp = pd.DataFrame.from_dict({data: getattr(self, data)})
                tmps.append(tmp)

        self.df_data = pd.concat(tmps, ignore_index=False, axis=1)

        for k, v in self.initial_guess_dict.items():
            self.df_data[f"{k}_ig"] = v

        self.df_data["c0_max_thresh"] = self.c0_max_thresh
        self.df_data["c0_min_thresh"] = self.c0_min_thresh

        for deet in self.curve_deets:
            if hasattr(self, deet):
                self.df_data[deet] = getattr(self, deet)

            for k, v in self.theta_dict.items():
                self.df_data[f"{k}_theta"] = v
            self.df_data["obj"] = self.obj

        if overwrite:
            self.df_data.to_csv(data_file_base, index=False)
        elif not overwrite:
            x = data_file_base.split(".csv")[0]
            append = 2
            while os.path.exists(data_file_base):
                data_file_base = x + f"_{append}.csv"
                append += 1
            self.df_data.to_csv(data_file_base, index=False)

    def _just_plot_curve(self):
        tmp = self.df_curve.reset_index()
        _, ax = plt.subplots(figsize=self.figsize)

        ax.plot(tmp.bv, tmp.c_norm, marker=".")
        ax.plot([0, tmp.bv.max()], [1, 1], linestyle=":")
        ax.set_title(
            f"Curve {self.curve_id}: {self.ref.replace('_', ' ').title()}, {self.compound.upper()}"
        )

        # ax.set_ylim([-0.01, 1.02])
        ax.set_xlabel("BV")
        ax.set_ylabel(f"C/C$_0$ [{self.compound.upper()}]")
        self.idx_labels = zip(
            tmp.bv.to_list(), tmp.c_norm.to_list(), tmp.index.to_list()
        )
        for x, y, i in self.idx_labels:
            # ax.annotate(f"{i}\n{y}", (x, y), fontsize="small", rotation=-45)
            ax.annotate(i, (x, y), fontsize="small")
        plt.tight_layout()
        print(f"INDEXES: {tmp.index.to_list()}")

    def _get_comp_list(self, blk, comp=Var, skip_list=[]):
        cs = []
        split_name = blk.name + "."
        skip_list += ["ref", "process_flow", "regeneration"]
        for c in blk.component_objects(comp):
            if any(s in c.name for s in skip_list):
                continue
            cs.append(c.name.split(split_name)[1])
        return cs

    # def _get_resin_deets(self, resin):
    #     return resin_dens_dict[resin], resin_diam_dict[resin]

    def _set_bounds(self, m=None):

        if m is None:
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)

        for k, v in self.set_bounds_dict.items():
            # print("set_bounds", k, v)
            ixv = getattr(ix, k)
            ixv.setlb(v[0])
            ixv.setub(v[1])

    def _fix_initial_guess(self, m=None):

        if m is None:
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)

        if isinstance(self.initial_guess_dict, dict):
            for k, v in self.initial_guess_dict.items():
                # print("initial_guess_dict", k, v)
                ixv = getattr(ix, k)
                ixv.fix(v)

    def _calc_from_constr(self, m=None):

        if m is None:
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)
        if not hasattr(self, "calc_from_constr_fail"):
            self.calc_from_constr_fail = dict()

        if isinstance(self.calc_from_constr_dict, dict):
            for k, v in self.calc_from_constr_dict.items():
                # print(k, "calc_from_constr_dict", v)
                ixv = getattr(ix, k)
                ixc = getattr(ix, v)
                if all(c.is_indexed() for c in [ixv, ixc]):
                    idx = [*ixv.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv[i]())
                        try:
                            calculate_variable_from_constraint(ixv[i], ixc[i])
                        except:
                            print(f"calc_from_constr_dict FAIL: {k, v, i, ixv[i]()}")
                            self.calc_from_constr_fail[k] = (i, v)
                elif ixv.is_indexed():
                    idx = [*ixv.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv[i]())
                        try:
                            calculate_variable_from_constraint(ixv[i], ixc)
                        except:
                            print(f"calc_from_constr_dict FAIL: {k, v, i, ixv[i]()}")
                            self.calc_from_constr_fail[k] = (i, v, ixv[i]())
                        # print(f'{k}[{i}] = {ixv[i]()}\n')

                elif ixc.is_indexed():
                    idx = [*ixc.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv())
                        try:
                            calculate_variable_from_constraint(ixv, ixc[i])
                        except:
                            print(f"calc_from_constr_dict FAIL: {k, v, i, ixv()}")
                            self.calc_from_constr_fail[k] = (i, v, ixv())
                        # print(f'{k}[{i}] = {ixv()}\n')

                else:
                    # print(k, v, i, ixv())
                    try:
                        calculate_variable_from_constraint(ixv, ixc)
                    except:
                        print(f"calc_from_constr_dict FAIL: {k, v, i, ixv()}")
                        self.calc_from_constr_fail[k] = (i, v, ixv())
                    # print(f'{k} = {ixv()}\n')

    def _get_ix_blk(self, m=None):
        if m is None:
            m = self.m
        for blk in m.fs.component_objects(Block):
            if isinstance(blk, self.ix_model):
                return blk

    def _linear_fit(self, bv):
        return self.bv50_slope * bv + self.bv50_int

    def _linear_inv(self, cb):
        return (cb - self.bv50_int) / self.bv50_slope


def main():
    set_bounds = {
        "mass_transfer_coeff": [0, None],
        "freundlich_n": [
            1.05,
            50,
        ],  # recommend keeping at least the lb to help in stable solves
        "service_flow_rate": [0, None],
        "loading_rate": [0, None],
        "ebct": [0, None],
        "resin_diam": [0, None],
    }

    calc_from_constr_dict = {
        "bed_vol_tot": "eq_bed_flow",
        "col_diam": "eq_bed_design",
        "col_height": "eq_col_height",
        "resin_surf_per_vol": "eq_resin_surf_per_vol",
        "service_flow_rate": "eq_service_flow_rate",
        "ebct": "eq_ebct",
        "N_Re": "eq_Re",
        "N_Sc": "eq_Sc",
        "N_Sh": "eq_Sh",
        "N_Pe_bed": "eq_Pe_bed",
        "N_Pe_particle": "eq_Pe_p",
    }
    calc_from_constr_dict = dict()
    initial_guess_dict = {
        "mass_transfer_coeff": 0.1,
        "freundlich_n": 3,
    }
    scale_from_value = [
        "bed_vol_tot",
        "col_diam",
        "col_height",
        "resin_surf_per_vol",
        "service_flow_rate",
        "ebct",
        "N_Re",
        "N_Sc",
        "N_Sh",
        "N_Pe_bed",
        "N_Pe_particle",
        "bv_50",
        "col_height_to_diam_ratio",
        "bv",
    ]
    curve_id = 2
    c = IXParmest(
        curve_id,
        regenerant="single_use",
        c0_min_thresh=0,
        initial_guess_dict=initial_guess_dict,
        set_bounds=set_bounds,
        calc_from_constr_dict=calc_from_constr_dict,
        scale_from_value=scale_from_value,
    )
    # c.estimate_bv50()
    # c.plot_estimate_bv50()

    c.run_all_things(plot_things=True, save_things=False)


if __name__ == "__main__":
    main()
    plt.show()
