import os
import math
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


from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.ion_exchange_0D import IonExchange0D
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)

from copy import deepcopy

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from IPython.display import clear_output

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
f = os.path.join(__location__, "ix_pfas_data.csv")
df_pfas = pd.read_csv(f)


mw_dict = {
    "pfba_-": 214.04e-3,
    "pfhpa_-": 314.05e-3,
    "pfhxa_-": 314.05e-3,
    "pfoa_-": 414.07e-3,
    "pfpea_-": 264.05e-3,
    "pfpra_-": 164.031e-3,
    "pfos_-": 500.13e-3,
    "pfbs_-": 300.1e-3,
    "pfhxs_-": 400.12e-3,
    "pfna_-": 464.08e-3,
    "pfda_-": 514.08e-3,
}

resin_fg_dict = {  # functional group
    "a694": "Complex Amino",
    "a600": "Quaternary Ammonium",
    "a520e": "Quaternary Ammonium",
    "ect_sorbix_a3f": "Unknown",  # no data, assumed
    "ect2_sorbix_lc4": "Unknown",  # no data, assumed
    "calgon_calres_2301": "N-Tributyl Amine",
    "amberlite_psr2": "N-Tributyl Amine",
    "aer1_ground": "Quaternary Amine",
    "aer1": "Quaternary Amine",
    "aer2_ground": "Quaternary Ammonium",
    "aer2": "Quaternary Ammonium",
    "amberlite_ira67": "Tertiary Amine",
    "amberlite_ira410": "Dimethyl Ethanol Amine",
    "a860": "Quaternary Ammonium",
    "amberlite_ira900": "Trimethyl Ammonium",
    "a592e": "Complex Amino",
}

resin_type_dict = {
    "a694": "Strong Base Type I Gel Polystyrene-Divinylbenzene",
    "a600": "Strong Base Type I Gel Polystyrene-Divinylbenzene",
    "a520e": "Strong Base Type I Macroporous Polystyrene-Divinylbenzene",
    "ect_sorbix_a3f": "Unknown",  # no data, assumed
    "ect2_sorbix_lc4": "Unknown",  # no data, assumed
    "calgon_calres_2301": "Macroporous",
    "amberlite_psr2": "Strong Base Type I Gel Polystyrene-Divinylbenzene",
    "aer1_ground": "Gel Polystyrene-Divinylbenzene",
    "aer1": "Gel Polystyrene-Divinylbenzene",
    "aer2_ground": "Macroporous Polystyrene-Divinylbenzene",
    "aer2": "Macroporous Polystyrene-Divinylbenzene",
    "amberlite_ira67": "Weak Base Gel Polyacrylic",
    "amberlite_ira410": "Strong Base Type II Gel Polystyrene-Divinylbenzene",
    "a860": "Strong Base Type I Macroporous Polyacrylic-Divinylbenzene",
    "amberlite_ira900": "Strong Base Type I Macroporous Polystyrene-Divinylbenzene",
    "a592e": "Strong Base Type I Macroporous Polystyrene-Divinylbenzene",
}

resin_dens_dict = {  # density of resin
    "a694": 0.72,
    "a600": 0.72,
    "a520e": 0.69,
    "ect_sorbix_a3f": 0.7,  # no data, assumed
    "ect2_sorbix_lc4": 0.7,  # no data, assumed
    "calgon_calres_2301": 0.67,
    "amberlite_psr2": 0.69,
    "aer1_ground": 0.7,  # no data, assumed
    "aer1": 0.7,  # no data, assumed
    "aer2_ground": 0.7,  # no data, assumed
    "aer2": 0.7,  # no data, assumed
    "amberlite_ira67": 0.65,
    "amberlite_ira410": 0.68,
    "a860": 0.71,
    "amberlite_ira900": 0.7,
    "a592e": 0.68,
}

resin_diam_dict = {  # diameter of resin bead
    "a694": 675e-6,
    "a600": 570e-6,
    "a520e": 650e-6,  # 300-1200 um range, uniformity coeff = 1.7
    "ect_sorbix_a3f": 675e-6,  # no data, assumed
    "ect2_sorbix_lc4": 675e-6,  # no data, assumed
    "calgon_calres_2301": 530e-6,  # "effective size" + std
    "amberlite_psr2": 700e-6,
    "aer1_ground": 0.26e-3,
    "aer1": 0.81e-3,
    "aer2_ground": 0.2e-3,
    "aer2": 0.73e-3,
    "amberlite_ira67": 625e-6,
    "amberlite_ira410": 675e-6,
    "a860": 650e-6,  # 300-1200 um range, uniformity coeff = 1.7
    "amberlite_ira900": 720e-6,
    "a592e": 675e-6,
}


class IXParmest:
    def __init__(
        self,
        curve,
        df_curve=None,
        bv=None,
        cb=None,
        isotherm="freundlich",
        dont_calc_effluent=True,
        fix_initial_guess=dict(),
        theta_initial_guess=dict(),
        set_bounds=dict(),
        calc_from_constr=dict(),
        autoscale_fixed=True,
        scale_from_value=None,
        fix_vars=dict(),
        c0_min_thresh=0.01,  # for determining keep_bvs + keep_cbs
        c0_max_thresh=0.99,  # for determining keep_bvs + keep_cbs
        c_next_thresh=1.05,  # for determining keep_bvs + keep_cbs
        cb50_min_thresh=None,  # all cb > cb50_min_thresh are used to make linear regression estimate bv_50
        figsize=(7, 5),
        **kwargs,
    ):
        self.curve = curve
        if df_curve is None:
            self.df_curve = df_pfas[df_pfas.curve == self.curve].copy()
        else:
            # option to provide separate DataFrame with breakthrough curve info
            self.df_curve = df_curve
        self.figsize = figsize
        self.bv = bv  # need value for initial model build
        self.cb = cb  # need value for initial model build
        self.all_figs = dict()  # dict for storing all figs
        self.isotherm = isotherm  # CONFIG option for IonExchange0D
        self.dont_calc_effluent = dont_calc_effluent  # boolean to indicate if ss effluent concentrations are deactivated
        if fix_initial_guess == dict():
            self.fix_initial_guess = {  # default initial guesses for regressed parameters
                "kinetic_param": 1e-6,  # by default, bv_50 has no initial guess. Default initial guess for bv_50 is avg of BV values from curve.
                "freundlich_n": 1.1,
            }
        else:
            self.fix_initial_guess = deepcopy(fix_initial_guess)
        if set_bounds == dict():
            # dict to set_bounds of IonExchange0D variables
            self.set_bounds = {
                "kinetic_param": [1e-10, None],
                "freundlich_n": [
                    1.05,
                    50,
                ],  # recommend keeping at least the lb to help in stable solves
                "service_flow_rate": [0, 500],
                "vel_bed": [0, None],
                "ebct": [0, None],
                "t_contact": [0, None],
            }
        else:
            self.set_bounds = deepcopy(set_bounds)
        if calc_from_constr == dict():
            # dict to calculate variables from constraint prior to initializaton
            self.calc_from_constr = {
                "c_norm": "eq_c_breakthru",
                "mass_transfer_coeff": "eq_mass_transfer_coeff",
                "t_breakthru": "eq_bv",
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
                # "bed_capacity_param": "eq_clark_2",
            }
        else:
            self.calc_from_constr = deepcopy(calc_from_constr)
        if theta_initial_guess == dict():
            # if wanted to provide a different initial guess for theta
            self.theta_initial_guess = deepcopy(self.fix_initial_guess)

        if len(self.fix_initial_guess) == 0 and len(self.theta_initial_guess) == 0:
            raise ValueError("need to provide initial guess")
        # if autoscale_fixed is True, will automatically scale fixed variables with 1 / value(var)
        self.autoscale_fixed = autoscale_fixed
        # min/max thresholds for filtering data to be used in parmest
        self.c0_min_thresh = c0_min_thresh
        self.c0_max_thresh = c0_max_thresh
        self.c_next_thresh = c_next_thresh

        if cb50_min_thresh is None:
            # default for the min threshold for estimating BV50 is the same
            self.cb50_min_thresh = self.c0_min_thresh
        else:
            self.cb50_min_thresh = cb50_min_thresh

        if scale_from_value is None:
            # scale these variables with 1 / value(var)
            # default is to use the same variables provided in calc_from_constr dict
            # since these vars would have an initial value that is (in theory) very close to the solved value
            self.scale_from_value = [k for k in calc_from_constr.keys()]
            self.scale_from_value += ["bv"]
        else:
            self.scale_from_value = scale_from_value

        self.get_curve_conditions()

        if fix_vars == dict():
            self.fix_vars = None
        else:
            self.fix_vars = fix_vars

        self.rebuild()  # initial build of model
        self.m0 = self.m.create_instance()  # intact initial build of model
        self.R_sq = np.nan
        self.build_results_dict()  # build dict for storing results

        self.print_curve_conditions()

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
        clear_output(wait=True)
        self.test_initial_guess()
        clear_output(wait=True)
        self.run_parmest()
        clear_output(wait=True)
        self.test_theta()
        clear_output(wait=True)
        if plot_things:
            self.plot_curve()
            self.plot_estimate_bv50()
            self.plot_initial_guess()
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

    def estimate_bv50(self):
        def linear_fit(bv, slope, b):
            return slope * bv + b

        def linear_inv(cb, slope, b):
            return (cb - b) / slope

        if not hasattr(self, "bv_pred_ig"):
            self.test_initial_guess()

        self.bv50_bvs = [
            bv
            for bv, cb in zip(self.keep_bvs, self.keep_cbs)
            if cb >= self.cb50_min_thresh * self.c0
        ]

        self.bv50_cbs = [
            cb for cb in self.keep_cbs if cb >= self.cb50_min_thresh * self.c0
        ]

        cb_50 = 0.5 * self.c0

        self.bv50_fit_params, self.bv50_fit_cov = curve_fit(
            linear_fit, self.bv50_bvs, self.bv50_cbs
        )

        self.bv50_slope, self.bv50_int = [*self.bv50_fit_params]

        self.bv50_est = linear_inv(cb_50, self.bv50_slope, self.bv50_int)

        self.bv50_orig = self.fix_initial_guess["bv_50"]
        self.bv_pred_ig_orig = deepcopy(self.bv_pred_ig)
        self.cb_pred_ig_orig = deepcopy(self.cb_pred_ig)

        self.bv50_corr_matrix = np.corrcoef(
            self.bv50_bvs, [self._linear_fit(bv) for bv in self.bv50_bvs]
        )
        corr = self.bv50_corr_matrix[0, 1]
        self.bv50_R_sq = corr**2
        if self.bv50_est > 100:
            self.fix_initial_guess["bv_50"] = self.bv50_est
        else:
            self.fix_initial_guess["bv_50"] = self.bv50_orig
        self.theta_initial_guess = deepcopy(self.fix_initial_guess)
        # rebuild model with new initial guess

        self.rebuild()
        self.test_initial_guess()
        # self.plot_estimate_bv50()

    def plot_estimate_bv50(self):

        if self.bv50_est < 100:
            self.plot_initial_guess()
            return
        max_bv = max([self._linear_inv(self.c0)] + self.all_bvs)
        cb_50 = 0.5 * self.c0

        title = (
            f"BV50 Estimate Results Compare - Curve {self.curve}:\n"
            + self.compound.swapcase()
            + " "
            + self.ref.replace("_", " ").title()
            + " "
            + self.resin.replace("_", " ").title()
            + f", EBCT = {self.ebct_min} min"
        )

        textstr = "\n".join(
            [
                f"Guess",
                f"  r: {self.fix_initial_guess['kinetic_param']:.2e}",
                f"  n: {self.fix_initial_guess['freundlich_n']:.2f}",
                f"  BV50 (new): {round(self.fix_initial_guess['bv_50'])}",
                f"  BV50 (orig): {round(self.bv50_orig)}",
            ]
        )
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.plot(
            self.keep_bvs,
            self.keep_cbs,
            marker=".",
            color="red",
            alpha=0.2,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cbs,
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
            [self.c0 for _ in [0, max_bv]],
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
        ax.set_ylim([-0.01, self.c0 * 1.02])
        ax.set_xlim([0, max_bv])
        ax.legend(loc="lower right")
        plt.tight_layout()

        self.bv50_fig = fig
        self.all_figs["bv_50"] = {"fig": fig, "ax": ax}

    def clean_bvs_cbs(self):

        """
        Method to filter extracted breakthrough
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

        for i, cb in enumerate(self.all_cbs):
            if i != len(self.all_cbs) - 1:
                if (
                    cb > last_cb
                    and cb < self.c0
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
                    cb > last_cb
                    and cb < self.c0
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
            self.keep_cbs,
            color="red",
            marker=".",
            alpha=0.25,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cbs,
            color="purple",
            marker=".",
            alpha=0.25,
            label="All Data",
        )
        ax.plot(
            [0, max(self.all_bvs)],
            [self.c0, self.c0],
            linestyle="-.",
            alpha=0.25,
            label="Influent",
        )

        title = (
            f"Data Used for Curve {self.curve}\n"
            + self.compound.swapcase()
            + ", "
            + self.ref.replace("_", " ").title()
            + ", "
            + self.resin.replace("_", " ").title()
            + f", EBCT = {self.ebct_min} min"
        )
        ax.set_xlabel("BV")
        ylabe = self.compound.swapcase()
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.05, self.c0 * 1.02])
        ax.set_xlim([-5, max(self.all_bvs) * 1.05])
        ax.legend()
        plt.tight_layout()

        self.curve_fig = fig

        self.all_figs["raw_curve"] = {"fig": fig, "ax": ax}

        # return plot_curve(self.curve)

    def test_initial_guess(self):
        """
        Loops through the kept (bv, cb) pairs and runs the model
        with those values using fix_initial_guess
        """
        self.flag = "test_initial_guess"
        ix = self.ix
        self._fix_initial_guess()
        self._calc_from_constr()
        self.bv_pred_ig = []  # bv_pred_ig = Predicted BVs using initial guess
        self.bv_fail_ig = []
        self.cb_pred_ig = []  # cb_pred_ig = Predicted Cbs using initial guess
        self.cb_fail_ig = []
        for (bv, cb) in zip(self.keep_bvs, self.keep_cbs):
            ix.c_breakthru.fix(cb * self.conc_units)
            ix.bv.set_value(bv)
            try:
                ix.initialize()
            except:
                # self._calc_from_constr()
                pass
            self.solve_it()
            if self.tc == "optimal":
                self.bv_pred_ig.append(ix.bv())
                self.cb_pred_ig.append(cb)
            elif self.tc != "optimal":
                self.bv_fail_ig.append(bv)
                self.cb_fail_ig.append(cb)
            self.results_dict_append()
            clear_output(wait=True)

        # self.plot_initial_guess()

    def plot_initial_guess(self):

        fig, ax = plt.subplots(figsize=self.figsize)
        ax.scatter(
            self.bv_pred_ig,
            self.cb_pred_ig,
            marker="o",
            color="green",
            label="Predicted",
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
            self.keep_cbs,
            marker=".",
            color="red",
            alpha=0.25,
            label="Parmest Data",
        )
        ax.plot(
            self.all_bvs,
            self.all_cbs,
            marker=".",
            color="purple",
            alpha=0.25,
            label="All Data",
        )
        ax.plot(
            [0, max(self.bv_pred_ig + self.keep_bvs)],
            [self.c0, self.c0],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
        )
        ax.plot(
            [0, max(self.bv_pred_ig + self.keep_bvs)],
            [self.c0 * 0.5, self.c0 * 0.5],
            linestyle="-.",
            alpha=0.25,
            # label="Influent",
        )

        ax.scatter(
            [self.fix_initial_guess["bv_50"]],
            [self.c0 * 0.5],
            marker="*",
            color="springgreen",
            # label="BV50 Guess",
        )
        # ax.plot(test_bvs, tes t_cbs, linestyle="-.", alpha=0.25, label="Test Data")
        title = (
            f"Initial Guess - Curve {self.curve}:\n"
            + self.compound.swapcase()
            + " "
            + self.ref.replace("_", " ").title()
            + " "
            + self.resin.replace("_", " ").title()
            + f" EBCT = {self.ebct_min} min"
        )
        ax.set_xlabel("BV")
        ylabe = self.compound.swapcase() + f" [{self.conc_units_str}]"
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.01, self.c0 * 1.02])

        textstr = "\n".join(
            [
                f"Guess",
                f"  r: {self.fix_initial_guess['kinetic_param']:.2e}",
                f"  n: {self.fix_initial_guess['freundlich_n']:.2f}",
                f"  BV50: {round(self.fix_initial_guess['bv_50'])}",
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

        # if save_it:
        #     file = f"../figs/theta_initial_guess_curve{curve}.png"
        #     x = file.split(".png")[0]
        #     append = 2
        #     while os.path.exists(file):
        #         file = x + f"_{append}.png"
        #         append += 1
        #     fig.savefig(file)

    def get_curve_conditions(self):
        """
        Gets conditions for the breakthrough curve to use in model
        and other metadata
        """
        df = self.df_curve.copy()
        for col in df.columns:
            if col in [
                "bv",
                "cb",
                "c_norm",
            ]:
                continue
            # if not len(df[col].unique()) == 1:
            #     raise Exception(f"{col} for curve {df.curve.iloc[0]} is not uniform")

        self.ref = df.ref.unique()[0]
        self.c0 = df.c0.unique()[0]

        self.compound = df.compound.unique()[0]
        self.flow_in = df.flow_in.unique()[0]
        self.vel_bed = df.vel_bed.unique()[0]
        self.bed_depth = df.bed_depth.unique()[0] * pyunits.m
        self.bed_vol_tot = df.bed_vol.unique()[0] * pyunits.m**3
        self.resin = df.resin.unique()[0]
        self.resin_type = resin_type_dict[self.resin]
        self.resin_func_group = resin_fg_dict[self.resin]
        self.ebct_min = df.ebct.unique()[0]
        self.ebct = pyunits.convert(
            self.ebct_min * pyunits.min, to_units=pyunits.second
        )
        self.conc_units_str = df.conc_units.unique()[0]
        numer = getattr(pyunits, self.conc_units_str.split("/")[0])
        denom = getattr(pyunits, self.conc_units_str.split("/")[1])
        self.conc_units = numer / denom
        self.charge = -1
        self.target_ion = self.compound + "_-"
        self.mw = mw_dict[self.target_ion] * pyunits.kg / pyunits.mol
        self.clean_bvs_cbs()
        if ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2) > 0:
            self.expr_sf = 1 / ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2)

        self.resin_dens, self.resin_diam = self._get_resin_deets(self.resin)

        self.ion_props = {
            "solute_list": [self.target_ion],
            "diffusivity_data": {("Liq", self.target_ion): 0.49e-9},
            "mw_data": {"H2O": 0.018, self.target_ion: self.mw},
            "charge": {self.target_ion: self.charge},
        }

        self.ix_config = {
            "target_ion": self.target_ion,
            "isotherm": self.isotherm,
        }

        if np.isnan(self.vel_bed):
            self.vel_bed = self.bed_depth / value(self.ebct)

        self.c0_mol_flow = pyunits.convert(
            (self.c0 * self.conc_units)
            / self.mw
            * (self.flow_in * pyunits.m**3 / pyunits.s),
            to_units=pyunits.mol / pyunits.s,
        )
        self.c0_mol_flow_sf = 1 / value(self.c0_mol_flow)

        if self.bv is None and self.cb is None:
            self.cb = self.keep_cbs[1] * self.conc_units
            self.bv = self.keep_bvs[1]

        if "bv_50" not in self.fix_initial_guess.keys():
            self.fix_initial_guess["bv_50"] = np.mean(self.keep_bvs)
        if "bv_50" not in self.theta_initial_guess.keys():
            self.theta_initial_guess["bv_50"] = np.mean(self.keep_bvs)

    def build_it(self):
        """
        Method used to build IX model
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**self.ion_props)

        self.ix_config["property_package"] = m.fs.properties

        m.fs.ix = ix = IonExchange0D(**self.ix_config)

        pf = ix.process_flow
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", self.c0_mol_flow_sf, index=("Liq", self.target_ion)
        )

        c_in = pyunits.convert(
            self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3
        )
        c_b = pyunits.convert(self.cb, to_units=pyunits.kg / pyunits.m**3)

        pf.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): self.flow_in,
                ("conc_mass_phase_comp", ("Liq", self.target_ion)): c_in,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        # Parameters set to force bed_depth == col_height
        ix.underdrain_h.set_value(0)
        ix.distributor_h.set_value(0)
        ix.bed_expansion_frac_A.set_value(0)
        ix.bw_rate.set_value(0)
        ix.number_columns.fix(1)

        ix.c_breakthru[self.target_ion].fix(value(c_b))
        ix.resin_bulk_dens.fix(self.resin_dens)
        ix.resin_diam.fix(self.resin_diam)
        ix.bed_porosity.fix(0.4)
        ix.bed_depth.fix(self.bed_depth)
        ix.vel_bed.fix(self.vel_bed)

        ix.bv.setlb(0)
        ix.bv.setub(None)
        ix.bv.set_value(
            self.bv
        )  # since we want to predict BV, we use .set_value() rather than .fix()

        if self.dont_calc_effluent:
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
            ix.eq_mass_transfer_regen.deactivate()
            pf.mass_transfer_term[0, "Liq", self.target_ion].fix()

        if isinstance(self.fix_vars, dict):
            # If fix_vars is provided
            for k, v in self.fix_vars.items():
                getattr(m.fs.ix, k).fix(v)

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

        target_ion = ix.config.target_ion
        prop_out = ix.process_flow.properties_out[0]
        prop_regen = ix.regeneration_stream[0]
        set_scaling_factor(
            prop_out.flow_mol_phase_comp["Liq", target_ion], self.c0_mol_flow_sf
        )
        set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], 1)
        set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", target_ion], 1e15)

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

        set_scaling_factor(ix.c_traps, 1)
        # Still not clear if constraint_scaling_transform is helping or hurting here...
        constraint_scaling_transform(ix.eq_clark_2[target_ion], 1e6)
        constraint_scaling_transform(ix.eq_clark_1[target_ion], 1e6)
        constraint_scaling_transform(ix.eq_mass_transfer_regen[target_ion], 1e7)
        constraint_scaling_transform(ix.eq_mass_transfer_target_fr[target_ion], 1e7)

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

        self.theta_names = ["fs.ix.kinetic_param", "fs.ix.freundlich_n", "fs.ix.bv_50"]
        self.theta_input = dict(
            zip(self.theta_names, self.theta_initial_guess.values())
        )
        self.theta_values = pd.DataFrame(
            data=[[v for v in self.theta_initial_guess.values()]],
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
                f"CURVE = {self.curve}\n\t"
                f"REF = {self.ref}\n\t"
                f"COMPOUND = {self.compound}\n\t"
                f"BV = {bv}\n\t"
                f"C0 = {c0_p}\n\t"
                f"CB = {cb_p}\n\t"
                f"CNORM = {cx}\n\t"
            )

            m_parmest = self.m.create_instance()
            m_parmest.fs.ix.c_breakthru.fix(cb)
            m_parmest.fs.ix.bv.set_value(bv)
            self._calc_from_constr(m=m_parmest)
            self.scale_it(m=m_parmest)
            try:
                m_parmest.fs.ix.initialize()
            except:
                pass
            # Stores the instantiated model that is passed to parmest for debugging outside of the class
            self.m_parmest = m_parmest
            return m_parmest

        def SSE(m, data):
            """
            Objective function for parmest to minimize
            We want the SSE of BV observed/predicted to be small
            """
            expr = (float(data.bv) - m.fs.ix.bv) ** 2
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

        self.theta_dict = dict(zip(self.theta_initial_guess.keys(), [*self.theta]))

        if hasattr(self, "results_dict"):
            for k, v in self.theta_dict.items():
                self.results_dict[f"{k}_theta"] = [
                    v for _ in range(len(self.results_dict["curve"]))
                ]
            self.results_dict["obj"] = [
                self.obj for _ in range(len(self.results_dict["curve"]))
            ]

    def test_theta(
        self,
        min_test_cbs=None,
        max_test_cbs=0.99,
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
        if min_test_bvs is None:
            self.min_test_bvs = min(self.keep_bvs)
        else:
            self.min_test_bvs = min_test_bvs
        if max_test_bvs is None:
            self.max_test_bvs = max(self.all_bvs)
        else:
            self.max_test_bvs = max_test_bvs

        self.max_test_cbs = max_test_cbs
        self.num_pts = num_pts

        self.test_cbs = [
            self.c0 * x
            for x in np.linspace(self.min_test_cbs, self.max_test_cbs, self.num_pts)
        ]
        self.test_bvs = np.linspace(self.min_test_bvs, self.max_test_bvs, self.num_pts)

        # Comparing against fake data
        self.bv_pred_theta_test = (
            []
        )  # bv_pred_theta_test = Predicted BV using theta and regularly spaced test data
        self.bv_fail_theta_test = []
        self.cb_pred_theta_test = (
            []
        )  # cb_pred_theta_test = Predicted Cb using theta and regularly spaced test data
        self.cb_fail_theta_test = []

        self.m_theta = self.m.create_instance()
        self.ix = ix = self._get_ix_blk(m=self.m_theta)
        for k, v in self.theta_dict.items():
            ixv = getattr(ix, k)
            ixv.fix(v)

        self._calc_from_constr(m=self.m_theta)
        self.scale_it(m=self.m_theta)

        try:
            ix.initialize()
        except:
            pass

        assert degrees_of_freedom(self.m_theta) == 0

        self.flag = "test_theta_fake"
        for bv, cb in zip(self.test_bvs, self.test_cbs):
            ix.c_breakthru[self.target_ion].fix(
                pyunits.convert(
                    cb * self.conc_units,
                    to_units=pyunits.kg / pyunits.m**3,
                )()
            )
            calculate_variable_from_constraint(
                ix.c_norm[self.target_ion], ix.eq_c_breakthru[self.target_ion]
            )
            ix.bv.set_value(bv)
            # self._calc_from_constr(m=self.m_theta)
            # self.scale_it(m=self.m_theta)
            ix.initialize()
            # try:
            #     ix.initialize()
            # except:
            #     pass

            assert degrees_of_freedom(self.m_theta) == 0
            self.solve_it(m=self.m_theta)
            if self.tc != "optimal":
                self.bv_fail_theta_test.append(bv)
                self.cb_fail_theta_test.append(cb)
            elif self.tc == "optimal":
                self.bv_pred_theta_test.append(ix.bv())
                self.cb_pred_theta_test.append(cb)
            self.results_dict_append()

        self.bv_pred_theta = []
        self.bv_fail_theta = []
        self.keep_bv_theta = []
        self.cb_pred_theta = []
        self.cb_fail_theta = []

        self.flag = "test_theta_real"
        for bv, cb in zip(self.keep_bvs, self.keep_cbs):
            if cb == 0:
                cb = 1e-3 * self.c0
            ix.c_breakthru[self.target_ion].fix(
                pyunits.convert(
                    cb * self.conc_units,
                    to_units=pyunits.kg / pyunits.m**3,
                )()
            )
            calculate_variable_from_constraint(
                ix.c_norm[self.target_ion], ix.eq_c_breakthru[self.target_ion]
            )
            ix.bv.set_value(bv)
            # self._calc_from_constr(m=self.m_theta)
            # self.scale_it(m=self.m_theta)
            try:
                ix.initialize()
            except:
                pass

            assert degrees_of_freedom(self.m_theta) == 0
            self.solve_it(m=self.m_theta)
            if self.tc != "optimal":
                self.bv_fail_theta.append(bv)
                self.cb_fail_theta.append(cb)
            elif self.tc == "optimal":
                self.bv_pred_theta.append(ix.bv())
                self.keep_bv_theta.append(bv)
                self.cb_pred_theta.append(cb)
            self.results_dict_append()

        self.corr_matrix = np.corrcoef(self.keep_bv_theta, self.bv_pred_theta)
        corr = self.corr_matrix[0, 1]
        self.R_sq = corr**2
        clear_output(wait=True)

    def plot_theta(self):

        print("\n\nPLOT THETA\n\n")

        if not all(hasattr(self, attr) for attr in ["test_cbs", "test_bvs", "R_sq"]):
            self.test_theta()

        textstr = "\n".join(
            [
                f"Theta (R2: {self.R_sq:.3f})",
                f"  r: {self.theta_dict['kinetic_param']:.2e}",
                f"  n: {self.theta_dict['freundlich_n']:.2f}",
                f"  BV50: {round(self.theta_dict['bv_50'])}",
            ]
        )
        boxprops = dict(boxstyle="round", facecolor="wheat", alpha=0.25)

        fig1, ax1 = plt.subplots(figsize=self.figsize)
        # ax1, ax2 = [*axs]

        ax1.scatter(
            self.bv_pred_theta_test,
            self.cb_pred_theta_test,
            color="blue",
            marker=".",
            label="Predicted",
        )
        ax1.plot(
            self.test_bvs,
            self.test_cbs,
            color="red",
            linestyle="-.",
            alpha=0.1,
            label="Test Data",
        )

        if len(self.bv_fail_theta_test) != 0:
            ax1.scatter(
                self.bv_fail_theta_test,
                self.cb_fail_theta_test,
                marker="x",
                color="black",
                label="Infeasible",
            )

        fig2, ax2 = plt.subplots(figsize=self.figsize)
        ax2.plot(
            self.bv_pred_theta_test,
            self.cb_pred_theta_test,
            color="blue",
            marker=".",
            alpha=0.1,
            # label="Predicted",
        )
        ax2.scatter(
            self.bv_pred_theta,
            self.cb_pred_theta,
            marker=".",
            color="green",
            label="Predicted",
        )
        ax2.plot(
            self.keep_bvs,
            self.keep_cbs,
            marker=".",
            color="red",
            alpha=0.2,
            label="Parmest Data",
        )

        ax2.plot(
            self.all_bvs,
            self.all_cbs,
            marker=".",
            color="purple",
            alpha=0.2,
            label="All Data",
        )

        if len(self.bv_fail_theta) != 0:
            ax2.scatter(
                self.bv_fail_theta,
                self.cb_fail_theta,
                marker="x",
                color="black",
                label="Infeasible",
            )

        ylabe = self.compound.swapcase() + f" [{self.conc_units_str}]"
        title = (
            f"Curve {self.curve}:\n"
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
                    self.c0
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
                label="Influent",
                color="black",
            )

            ax.set_xlabel("BV")
            ax.set_ylabel(ylabe)
            ax.set_title(title)
            ax.set_ylim([-0.01, self.c0 * 1.02])
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

        # file = f"../figs/theta_test_curve{curve}.png"
        # fig1.savefig(file)

        # fig2, axs2 = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        # ax3, ax4 = [*axs2]

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
        ax3.set_title(f"Curve {self.curve} Parity")
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

        # ax4.scatter(keep_bvs, bv_pred2, marker=".", color="green", label="Predicted")
        ax4.scatter(self.keep_bv_theta, self.bv_rel_error, color="black")
        ax4.set_xlabel("Actual BV")
        ax4.set_ylabel("Relative Error")
        ax4.set_title(f"Curve {self.curve} Relative Error")
        plt.tight_layout()

        self.all_figs["theta_relerror"] = {"fig": fig4, "ax": ax4}
        fig4.show()
        # file = f"../figs/theta_parity_curve{curve}.png"
        # fig2.savefig(file)

    def print_curve_conditions(self):
        clear_output(wait=False)
        print(f"\nCURVE: {self.curve}")
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
            "curve",
            "ref",
            "compound",
            "resin",
            "resin_type",
            "resin_func_group",
            "target_ion",
            "conc_units",
            "ebct_min",
            "flow_in",
            "charge",
            "c0",
            "expr_sf",
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

        for k in self.fix_initial_guess.keys():
            self.results_dict[f"{k}_ig"] = []

        for k in self.theta_initial_guess.keys():
            self.results_dict[f"{k}_ig_theta"] = []

        for k in self.calc_from_constr.keys():
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

            for k, v in self.fix_initial_guess.items():
                self.results_dict[f"{k}_ig"].append(v)

            for k, v in self.theta_initial_guess.items():
                self.results_dict[f"{k}_ig_theta"].append(v)

            for k, v in self.calc_from_constr.items():
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
            self.results_filename = f"results/curve{self.curve}_results.csv"
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

    def save_figs(self, overwrite=False, figs_filename=None):

        fig_file_base = f"figs/curve{self.curve}_"

        for figname, d in self.all_figs.items():
            fig = d["fig"]
            fig_file = fig_file_base + f"{figname}.png"
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

        data_file_base = f"output/curve{self.curve}_output.csv"

        tmps = []
        for data in sorted(datas):
            if hasattr(self, data):
                tmp = pd.DataFrame.from_dict({data: getattr(self, data)})
                tmps.append(tmp)

        self.df_data = pd.concat(tmps, ignore_index=False, axis=1)

        for k, v in self.fix_initial_guess.items():
            self.df_data[f"{k}_ig"] = v

        self.df_data["c0_max_thresh"] = self.c0_max_thresh
        self.df_data["c0_min_thresh"] = self.c0_min_thresh

        for deet in self.curve_deets:
            if hasattr(self, deet):
                self.df_data[deet] = getattr(self, deet)

        if hasattr(self, "theta_dict"):
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

    def _get_comp_list(self, blk, comp=Var, skip_list=[]):
        cs = []
        split_name = blk.name + "."
        skip_list += ["ref", "process_flow", "regeneration"]
        for c in blk.component_objects(comp):
            if any(s in c.name for s in skip_list):
                continue
            cs.append(c.name.split(split_name)[1])
        return cs

    def _get_resin_deets(self, resin):
        return resin_dens_dict[resin], resin_diam_dict[resin]

    def _set_bounds(self):

        for k, v in self.set_bounds.items():
            # print("set_bounds", k, v)
            ixv = getattr(self.ix, k)
            ixv.setlb(v[0])
            ixv.setub(v[1])

    def _fix_initial_guess(self):

        if isinstance(self.fix_initial_guess, dict):
            for k, v in self.fix_initial_guess.items():
                # print("fix_initial_guess", k, v)
                ixv = getattr(self.ix, k)
                ixv.fix(v)

    def _calc_from_constr(self, m=None):
        def calc_bed_capacity_param():
            c0 = ix.process_flow.properties_in[0].conc_mass_phase_comp[
                "Liq", self.target_ion
            ]()
            cb = ix.c_breakthru[self.target_ion]()
            return value(
                ((c0 / cb) ** (ix.freundlich_n() - 1) - 1)
                / math.exp(
                    (-ix.kinetic_param() * ix.bed_depth() * ix.bv()) / ix.vel_bed()
                )
            )

        if m is None:
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)
        if not hasattr(self, "calc_from_constr_fail"):
            self.calc_from_constr_fail = dict()

        if isinstance(self.calc_from_constr, dict):
            for k, v in self.calc_from_constr.items():
                # print(k, "calc_from_constr", v)
                ixv = getattr(ix, k)
                ixc = getattr(ix, v)
                if all(c.is_indexed() for c in [ixv, ixc]):
                    idx = [*ixv.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv[i]())
                        try:
                            calculate_variable_from_constraint(ixv[i], ixc[i])
                        except:
                            print(f"calc_from_constr FAIL: {k, v, i, ixv[i]()}")
                            self.calc_from_constr_fail[k] = (i, v)
                elif ixv.is_indexed():
                    idx = [*ixv.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv[i]())
                        try:
                            calculate_variable_from_constraint(ixv[i], ixc)
                        except:
                            print(f"calc_from_constr FAIL: {k, v, i, ixv[i]()}")
                            self.calc_from_constr_fail[k] = (i, v, ixv[i]())
                        # print(f'{k}[{i}] = {ixv[i]()}\n')

                elif ixc.is_indexed():
                    idx = [*ixc.index_set()]
                    for i in idx:
                        # print(k, v, i, ixv())
                        try:
                            calculate_variable_from_constraint(ixv, ixc[i])
                        except:
                            print(f"calc_from_constr FAIL: {k, v, i, ixv()}")
                            self.calc_from_constr_fail[k] = (i, v, ixv())
                        # print(f'{k}[{i}] = {ixv()}\n')

                else:
                    # print(k, v, i, ixv())
                    try:
                        calculate_variable_from_constraint(ixv, ixc)
                    except:
                        print(f"calc_from_constr FAIL: {k, v, i, ixv()}")
                        self.calc_from_constr_fail[k] = (i, v, ixv())
                    # print(f'{k} = {ixv()}\n')
        if "bed_capacity_param" not in self.calc_from_constr.keys():
            try:
                self.bed_capacity_param_calc = calc_bed_capacity_param()
                ix.bed_capacity_param.set_value(calc_bed_capacity_param())
            except ZeroDivisionError:
                ix.bed_capacity_param.set_value(1e-8)

    def _get_ix_blk(self, m=None):
        if m is None:
            m = self.m
        for blk in m.fs.component_objects(Block):
            if isinstance(blk, IonExchange0D):
                return blk

    def _linear_fit(self, bv):
        return self.bv50_slope * bv + self.bv50_int

    def _linear_inv(self, cb):
        return (cb - self.bv50_int) / self.bv50_slope


if __name__ == "__main__":
    curve = 2
    fix_initial_guess = {
        "kinetic_param": 5e-6,
        "freundlich_n": 2,
    }
    c = IXParmest(curve, fix_intial_guess=fix_initial_guess)
    c.run_all_things()
    # plt.show()
