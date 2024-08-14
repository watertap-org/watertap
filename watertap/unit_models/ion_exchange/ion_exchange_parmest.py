import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from scipy.optimize import curve_fit

from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Block,
    Expression,
    Objective,
    SolverFactory,
    Suffix,
    ComponentUID,
    maximize,
    minimize,
    units as pyunits,
)
import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    get_scaling_factor,
    constraint_scaling_transform,
    set_variable_scaling_from_current_value,
    badly_scaled_var_generator,
    report_scaling_issues,
    unset_scaling_factor,
)
from idaes.core.util.exceptions import InitializationError
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate

from watertap.core.util.model_diagnostics.infeasible import (
    print_infeasible_constraints,
    print_variables_close_to_bounds,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.ion_exchange.util import (
    plot_initial_guess,
    plot_theta,
    plot_estimate_bv50,
    plot_curve,
    build_results_dict,
    results_dict_append,
    save_results,
    make_results_df,
    save_output,
    make_output_df,
    save_figs,
)
from watertap.core.solvers import get_solver
from .ion_exchange_clark import IonExchangeClark
from .ion_exchange_thomas import IonExchangeThomas
from .ion_exchange_cphsdm import IonExchangeCPHSDM

# thetas = ["freundlich_n", "mass_transfer_coeff", "bv_50"]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

min_st_surrogate = PysmoSurrogate.load_from_file(
    # surr_dir + "/min_st_pysmo_surr_linear.json",
    "/Users/ksitterl/Documents/Python/nawi-analysis/NAWI-analysis/analysis_waterTAP/analysisWaterTAP/analysis_scripts/pfas_ix_cphsdm/surrogates/min_st_pysmo_surr_linear.json"
    # "/Users/ksitterl/Documents/Python/nawi-analysis/NAWI-analysis/analysis_waterTAP/analysisWaterTAP/analysis_scripts/pfas_ix_2/trained_surrogate_models/min_st_pysmo_surr_spline.json"
)
throughput_surrogate = PysmoSurrogate.load_from_file(
    "/Users/ksitterl/Documents/Python/nawi-analysis/NAWI-analysis/analysis_waterTAP/analysisWaterTAP/analysis_scripts/pfas_ix_cphsdm/surrogates/throughput_pysmo_surr_linear.json"
)
mw_water = 0.018 * pyunits.kg / pyunits.mol


class IXParmest:
    def __init__(
        self,
        curve_id=None,
        input_data=None,
        ix_model=IonExchangeCPHSDM,
        regenerant="single_use",
        diff_calculation="HaydukLaudie",
        data_file=None,
        resin_data=dict(),
        component_data=dict(),
        bv=None,
        cb=None,
        c_norm=None,
        thetas=None,
        theta_names=None,
        fix_vars_dict=dict(),
        initial_guess_dict=dict(),
        initial_guess_theta_dict=dict(),
        set_bounds_dict=dict(),
        calc_from_constr_dict=dict(),
        autoscale_fixed=True,
        scale_from_value=None,
        parmest_kwargs=dict(),
        c0_min_thresh=0.01,  # for determining keep_bvs + keep_cbs
        c0_max_thresh=0.9999,  # for determining keep_bvs + keep_cbs
        c_next_thresh=1.0,  # for determining keep_bvs + keep_cbs
        cb50_min_thresh=None,  # all cb > cb50_min_thresh are used to make linear regression estimate bv_50
        max_zero=1e-3,  # replace all values of c_norm < c0_min_thresh with this value
        min_one=0.9999,  # replace all values of c_norm > c0_max_thresh with this value
        use_all_data=False,  # eliminate filtering data and use all points
        use_this_data=None,  # list of indexes in filtered_data corresponding to data points to use
        figsize=(7, 5),
        just_plot_curve=False,
        save_directory=None,
        solver=None,
        expr_sf=None,
        rho=1000,  # kg/m3
        c0_mol_flow_sf=None,
        water_mol_flow_sf=None,
        tee=False,
        scale_var_dict=dict(),
        scale_constr_dict=dict()
    ):

        if input_data is None and data_file is None:
            raise ValueError(
                "No input data provided."
                "Must provide either DataFrame of data via input_data keyword"
                " or path to input file via data_file keyword."
            )

        if input_data is None:
            curve_data = pd.read_csv(data_file)
            if curve_id is None:
                assert len(curve_data.curve_id.unique()) == 1
                self.curve_id = curve_data.curve_id.unique()[0]
            else:
                self.curve_id = curve_id
            self.input_data = curve_data[curve_data.curve_id == self.curve_id].copy()
        else:
            # option to provide separate DataFrame with breakthrough curve_id info
            self.input_data = input_data
            # assume the provided DataFrame only contains one curve
            assert len(self.input_data.curve_id.unique()) == 1
            # override curve_id if provided
            self.curve_id = self.input_data.curve_id.unique()[0]

        self.input_data["point_id"] = [x for x in range(len(self.input_data))]
        self.input_data.set_index("point_id", inplace=True)
        self.filtered_data = deepcopy(self.input_data)
        self.filtered_data.sort_values("bv", inplace=True)
        # self.filtered_data.dropna(inplace=True)
        self.ix_model = ix_model
        self.regenerant = regenerant
        self.fix_vars_dict = fix_vars_dict

        self.figsize = figsize
        self.bv = bv  # need value for initial model build
        self.cb = cb  # need value for initial model build
        self.c_norm = c_norm
        self.rho = rho * pyunits.kg / pyunits.m**3

        self.resin_data = resin_data
        self.component_data = component_data
        self.parmest_kwargs = parmest_kwargs
        self.just_plot_curve = just_plot_curve
        self.use_all_data = use_all_data
        self.use_this_data = use_this_data
        self.diff_calculation = diff_calculation
        self.c0_mol_flow_sf = c0_mol_flow_sf
        self.tee = tee
        self.water_mol_flow_sf = water_mol_flow_sf

        self.scale_var_dict = scale_var_dict
        self.scale_constr_dict = scale_constr_dict

        if solver is None:
            self.solver = get_solver()
        else:
            self.solver = solver

        self.all_figs = dict()  # dict for storing all figs

        if self.use_all_data and self.use_this_data:
            # can't use all data and some data
            # default to using "this" data if both are True
            self.use_all_data = False

        if initial_guess_dict == dict():
            if self.ix_model is IonExchangeClark:
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

        if (
            len(self.initial_guess_dict) == 0
            and len(self.initial_guess_theta_dict) == 0
        ):
            raise ValueError("need to provide initial guess")
        # if autoscale_fixed is True, will automatically scale fixed variables with 1 / value(var)
        self.autoscale_fixed = autoscale_fixed
        # min/max thresholds for filtering data to be used in parmest
        self.c0_min_thresh = c0_min_thresh
        self.c0_max_thresh = c0_max_thresh
        self.c_next_thresh = c_next_thresh
        self.expr_sf = expr_sf

        if thetas is None:
            if ix_model is IonExchangeClark:
                print("ix_model is  Clark")
                self.thetas = ["freundlich_n", "mass_transfer_coeff", "bv_50"]
            if ix_model is IonExchangeThomas:
                print("ix_model is Thomas")
                self.thetas = ["thomas_constant", "resin_max_capacity"]
            if ix_model is IonExchangeCPHSDM:
                print("ix_model is CPHSDM")
                # raise ValueError("")
                self.thetas = ["freundlich_k", "freundlich_ninv", "resin_max_capacity"]
        else:
            self.thetas = thetas

        if theta_names is None:
            print("theta_names is None")
            self.theta_names = self.thetas
        else:
            self.theta_names = theta_names

        if self.c0_min_thresh <= 0:
            self.max_zero = max_zero
            max_zero_i = self.filtered_data.query("c_norm <= 0").index.max()
            self.filtered_data.at[max_zero_i, "c_norm"] = max_zero
            self.filtered_data.at[max_zero_i, "cb"] = (
                max_zero * self.filtered_data.at[max_zero_i, "c0"]
            )
            self.c0_min_thresh = 0.9 * max_zero
        if 0 in self.filtered_data.c_norm.to_list():
            self.filtered_data = self.filtered_data.query("c_norm > 0")

        if self.c0_max_thresh >= 1:
            self.min_one = min_one
            min_one_i = self.filtered_data.query("c_norm >= 1").index.to_list()
            for i in min_one_i:
                self.filtered_data.at[i, "c_norm"] = min_one
                self.filtered_data.at[i, "cb"] = (
                    min_one * self.filtered_data.at[i, "c0"]
                )
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

        self.save_directory = save_directory

        self.dropped_points = [
            p for p in self.input_data.index if p not in self.filtered_data.index
        ]
        self.modified_points = [
            p
            for p in self.input_data.index
            if p not in self.dropped_points
            and not self.filtered_data.loc[p].equals(self.input_data.loc[p])
        ]

        self.get_curve_conditions()
        if self.just_plot_curve:
            self.plot_curve()
            tmp = self.input_data.reset_index()
            print(tmp.index.to_list())
            return

        self.get_model_config()
        # return

        self.rebuild()  # initial build of model
        # self.m0 = self.m.clone()  # intact initial build of model
        # self.R_sq = np.nan
        self.build_results_dict()  # build dict for storing results

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
        # self.test_initial_guess()
        self.run_parmest()
        self.test_theta()
        # if plot_things:
        #     self.plot_curve()
        #     self.plot_estimate_bv50()
        # self.plot_initial_guess()
        # self.plot_theta()
        # if save_things:
        #     if plot_things:
        #         self.save_figs(overwrite=overwrite)
        #     self.save_output(overwrite=overwrite)
        #     self.save_results(overwrite=overwrite)

    def activate_surrogate(self, m=None):

        if m is None:
            m = self.m
        ix = m.fs.ix
        ix.eq_min_number_st_cps.deactivate()
        ix.eq_throughput.deactivate()


        m.fs.min_st_surrogate = SurrogateBlock(concrete=True)
        m.fs.min_st_surrogate.build_model(
            min_st_surrogate,
            input_vars=[ix.freundlich_ninv, ix.N_Bi],
            # input_vars=[ix.freundlich_ninv, ix.Bi],
            output_vars=[ix.min_N_St],
        )
        m.fs.throughput_surrogate = SurrogateBlock(concrete=True)
        m.fs.throughput_surrogate.build_model(
            throughput_surrogate,
            input_vars=[ix.freundlich_ninv, ix.N_Bi, ix.c_norm],
            # input_vars=[ix.freundlich_ninv, ix.Bi, ix.c_norm],
            output_vars=[ix.throughput],
        )
        # ix.N_Bi.setlb(0)

    def rebuild(self):
        self.m0 = self.build_init()
        self.m = self.m0.clone()
        self._set_bounds()
        print("set_bounds")
        # self.m.fs.ix.N_Bi.display()
        self._calc_from_constr()
        self.scale_it()
        print("scale_it")
        # self.m.fs.ix.N_Bi.display()
        self.m.fs.ix.initialize()
        print("initialize")
        # self.m.fs.ix.N_Bi.display()
        self.activate_surrogate()
        print("activate_surrogate")
        # self.m.fs.ix.N_Bi.display()
        self.solve_it()
        print(f"initial standard build termination {self.tc}")
        self.build_design()
        print("build_design")
        # self.m.fs.ix.N_Bi.display()
        self.solve_it()
        print(f"initial design build termination {self.tc}")
        print(f"dof for design build = {degrees_of_freedom(self.m)}")

    def rebuild_test_theta(self):
        self.m0 = self.build_init()
        self.m = self.m0.clone()
        self._set_bounds()
        self._calc_from_constr()
        self.scale_it()
        self.m.fs.ix.initialize()
        self.activate_surrogate()
        # self.m.fs.ix.N_Bi.display()
        self.solve_it()
        print(f"initial standard build termination {self.tc}")
        self.build_design()
        for k, v in self.theta_dict.items():
            ixv = getattr(self.m.fs.ix, k)
            ixv.fix(v)
        self.solve_it()
        print(f"initial design build termination {self.tc}")
        print(f"dof for design build = {degrees_of_freedom(self.m)}")

    def test_initial_guess(
        self, m=None, min_cnorm_ig=0.05, max_cnorm_ig=0.95, num_pts=7
    ):
        """
        Loops through the kept (bv, cb) pairs and runs the model
        with those values using initial_guess_dict
        """
        print("\nTesting initial guess\n")

        # ix = self.ix
        # self._fix_initial_guess()
        # self._calc_from_constr()

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

        if m is None:
            m = self.m.clone()

        ix = m.fs.ix

        for cnorm in self.test_cnorms0:
            print(f"testing initial guess with test data:\n\tC/C0 = {cnorm}...")
            ix.c_norm.fix(cnorm)
            self.solve_it(m=m)
            if self.tc == "optimal":
                self.bv_pred_ig_fake.append(ix.bv())
                self.cb_pred_ig_fake.append(cnorm)
            elif self.tc != "optimal":
                # self.bv_fail_ig_fake.append(bv)
                self.cb_fail_ig_fake.append(cnorm)
            self.results_dict_append(ix_blk=ix)

        self.flag = "test_initial_guess"
        if m is None:
            m = self.m.clone()
        ix = m.fs.ix
        for bv, cnorm in zip(self.keep_bvs[::-1], self.keep_cnorms[::-1]):
            print(
                f"testing initial guess with actual data:\n\tBV = {bv}\n\tC/C0 = {cnorm}..."
            )
            ix.bv.set_value(bv)
            ix.c_norm.fix(cnorm)
            self.solve_it(m=m)
            print(f"termination {self.tc}")
            if self.tc == "optimal":
                self.bv_pred_ig.append(ix.bv())
                self.cb_pred_ig.append(cnorm)
            elif self.tc != "optimal":
                self.bv_fail_ig.append(bv)
                self.cb_fail_ig.append(cnorm)
            self.results_dict_append(ix_blk=ix)

    def get_curve_conditions(self):
        """
        Gets conditions for the breakthrough curve_id to use in model
        and other metadata
        """

        df = self.filtered_data.copy()

        # TODO: insert check for columns that should all have same value

        self.ref = df.ref.iloc[0]
        self.c0 = float(df.c0.iloc[0])
        # * pyunits.kg / pyunits.m**3 # kg/m3

        self.component = df.component.iloc[0]
        self.flow_in = float(df.flow_in.iloc[0])  # m3/s
        self.loading_rate = float(df.loading_rate.iloc[0])  # m/s
        self.bed_depth = float(df.bed_depth.iloc[0]) * pyunits.m
        self.bed_diameter = float(df.bed_diameter.iloc[0]) * pyunits.m
        self.bed_volume_total = float(df.bed_volume.iloc[0]) * pyunits.m**3
        self.bed_volume = float(df.bed_volume.iloc[0]) * pyunits.m**3
        self.resin = df.resin.iloc[0]
        self.ebct_min = float(df.ebct_min.iloc[0])  # minutes
        self.ebct = float(df.ebct.iloc[0])
        # self.conc_units_str = df.conc_units.iloc[0]  # conc_mass_phase_comp units
        # numer = getattr(pyunits, self.conc_units_str.split("/")[0])
        # denom = getattr(pyunits, self.conc_units_str.split("/")[1])
        # self.conc_units = numer / denom
        # self.c0_kg_m3 = pyunits.convert(self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3)

        self.target_component = self.component_data["name"]
        self.charge = self.component_data["charge"]
        self.mw = self.component_data["mw_comp"] * pyunits.kg / pyunits.mol
        if "molar_vol" in self.component_data.keys():
            self.molar_vol = self.component_data["molar_vol"]
        self.diffusivity = self.component_data["diffusivity"]
        self.filter_data()

        if self.expr_sf is None:
            try:
                if ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2) > 0:
                    self.expr_sf = 1 / ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2)
            except:
                self.expr_sf = 1e-9

        self.resin_density = self.resin_data["density"]
        self.resin_diam = self.resin_data["diameter"]
        self.resin_porosity = self.resin_data["porosity"]
        # self.bed_porosity = self.resin_data["porosity"]

    def get_model_config(self):

        if self.diff_calculation != "HaydukLaudie":
            self.ion_props = {
                "solute_list": [self.target_component],
                "diffusivity_data": {("Liq", self.target_component): self.diffusivity},
                "mw_data": {"H2O": mw_water, self.target_component: self.mw},
                "charge": {self.target_component: self.charge},
            }
        else:
            self.ion_props = {
                "solute_list": [self.target_component],
                "mw_data": {"H2O": mw_water, self.target_component: self.mw},
                "diffus_calculation": self.diff_calculation,
                "charge": {self.target_component: self.charge},
                "molar_volume_data": {("Liq", self.target_component): self.molar_vol},
            }

        self.ix_config = {
            "target_component": self.target_component,
            "regenerant": self.regenerant,
        }

        if np.isnan(self.loading_rate):
            self.loading_rate = self.bed_depth / value(self.ebct)

        self.c0_mol_flow = pyunits.convert(
            (self.c0 * pyunits.kg / pyunits.m**3)
            / self.mw
            * (self.flow_in * pyunits.m**3 / pyunits.s),
            to_units=pyunits.mol / pyunits.s,
        )
        if self.c0_mol_flow_sf is None:
            self.c0_mol_flow_sf = 1 / value(self.c0_mol_flow)
        self.water_mol_flow = pyunits.convert(
            ((self.flow_in * pyunits.m**3 / pyunits.s) * self.rho) / mw_water,
            to_units=pyunits.mol / pyunits.s,
        )
        if self.water_mol_flow_sf is None:
            self.water_mol_flow_sf = 1 / value(self.water_mol_flow)

        if self.cb is None:
            self.cb = self.keep_cbs[-1]  # kg/m3
        if self.c_norm is None:
            self.c_norm = self.keep_cnorms[-1]
        if self.bv is None:
            self.bv = self.keep_bvs[-1]

        if self.ix_model is IonExchangeClark:
            if "bv_50" not in self.initial_guess_dict.keys():
                self.initial_guess_dict["bv_50"] = np.mean(self.keep_bvs)
            if "bv_50" not in self.initial_guess_theta_dict.keys():
                self.initial_guess_theta_dict["bv_50"] = np.mean(self.keep_bvs)

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

        df = self.filtered_data.copy()
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
            tmp = tmp[tmp.c_norm < 1].copy()
            self.keep_bvs = [bv for bv in tmp.bv.to_list() if not np.isnan(bv)]
            self.keep_cbs = [cb for cb in tmp.cb.to_list() if not np.isnan(cb)]
            self.data_filter = "use_all_data"
            # assert 0 not in self.keep_cbs
        elif self.use_this_data is not None:
            print("use this data")
            tmp = df.reset_index(drop=True)
            tmp = tmp.loc[self.use_this_data].copy()
            self.keep_bvs = [bv for bv in tmp.bv.to_list() if not np.isnan(bv)]
            self.keep_cbs = [cb for cb in tmp.cb.to_list() if not np.isnan(cb)]
            self.data_filter = str(self.use_this_data)
            # assert 0 not in self.keep_cbs

        else:
            print("default filtering")
            self.data_filter = "default_filtering"

            for i, cb in enumerate(self.all_cbs):
                if np.isnan(cb):
                    continue
                # if cb > df.c0.iloc[0]:
                #     self.excl_cbs.append(cb)
                #     self.excl_bvs.append(self.all_bvs[i])
                #     continue

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
        assert 0 not in self.keep_cbs

        self.keep_cnorms = [cb / self.c0 for cb in self.keep_cbs if not np.isnan(cb)]
        self.excl_cnorms = [cb / self.c0 for cb in self.excl_cbs if not np.isnan(cb)]

    def build_init(self):
        """
        Method used to build IX model
        """

        self.flow_in_init = 2e-8
        self.c0_init = 4e-9
        self.mw_init = 0.5 * pyunits.kg / pyunits.mol

        if self.diff_calculation != "HaydukLaudie":
            self.ion_props_init = {
                "solute_list": [self.target_component],
                "diffusivity_data": {("Liq", self.target_component): self.diffusivity},
                "mw_data": {"H2O": mw_water, self.target_component: self.mw},
                "charge": {self.target_component: self.charge},
            }
        else:
            self.ion_props_init = {
                "solute_list": [self.target_component],
                "mw_data": {"H2O": mw_water, self.target_component: self.mw},
                "diffus_calculation": self.diff_calculation,
                "charge": {self.target_component: self.charge},
                "molar_volume_data": {("Liq", self.target_component): self.molar_vol},
            }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**self.ion_props)
        # m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
        # m.fs.properties.dens_mass_const = 1000

        self.ix_config["property_package"] = m.fs.properties

        self.c0_mol_flow_init = pyunits.convert(
            (self.c0_init * pyunits.kg / pyunits.m**3)
            / self.mw
            * (self.flow_in_init * pyunits.m**3 / pyunits.s),
            to_units=pyunits.mol / pyunits.s,
        )
        self.c0_mol_flow_sf_init = 1 / value(self.c0_mol_flow_init)
        self.water_mol_flow_init = pyunits.convert(
            ((self.flow_in_init * pyunits.m**3 / pyunits.s) * self.rho) / mw_water,
            to_units=pyunits.mol / pyunits.s,
        )
        self.water_mol_flow_sf_init = 1 / self.water_mol_flow_init()

        m.fs.ix = ix = self.ix_model(**self.ix_config)

        self._set_bounds(m=m)

        pf = ix.process_flow

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp",
            # self.c0_mol_flow_sf_init,
            1e15,
            index=("Liq", self.target_component),
        )

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp",
            1e5,
            # self.water_mol_flow_sf_init,
            index=("Liq", "H2O"),
        )
        calculate_scaling_factors(m)
        # pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)
        # ix.regeneration_stream[0].flow_mol_phase_comp["Liq", self.target_component].fix(0)
        pf.properties_in[0].pressure.fix(101325)
        pf.properties_in[0].temperature.fix(298)
        pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].set_value(
            self.water_mol_flow_init
        )
        pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].set_value(
            self.c0_mol_flow_init
        )

        pf.properties_in[0].conc_mass_phase_comp
        pf.properties_in[0].flow_vol_phase
        pf.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): self.flow_in_init,
                ("conc_mass_phase_comp", ("Liq", self.target_component)): self.c0_init,
            },
            hold_state=True,
        )
        # report_scaling_issues(m, prefix="\n")

        if self.ix_model is IonExchangeClark:
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
            ix.bed_diameter.fix(self.bed_diameter)

            # since we want to predict BV, we use .set_value() rather than .fix()
            ix.bv.set_value(self.bv)
            ix.loading_rate.set_value(self.loading_rate)
            ix.bed_volume_total.set_value(self.bed_volume_total)

            # Since we are fitting data to C/C0,
            # we don't care about the steady-state effluent concentration, so deactivate them
            # if self.ix_model is IonExchangeClark:
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
            ix.eq_mass_transfer_term.deactivate()

        if self.ix_model is IonExchangeCPHSDM:
            ix.eq_isothermal_regen_stream.deactivate()
            ix.eq_Pe_bed.deactivate()
            ix.eq_Pe_p.deactivate()
            ix.eq_service_flow_rate.deactivate()
            ix.eq_column_height.deactivate()
            ix.eq_bed_design.deactivate()

            ix.number_columns.fix(1)
            ix.freundlich_ninv.fix(0.9)
            ix.freundlich_k.fix(1)
            ix.surf_diff_coeff.fix(1e-13)
            ix.c_norm.fix(0.5)

            ix.shape_correction_factor.fix(1)
            ix.tortuosity.fix(1)
            ix.shape_correction_factor.fix(1)
            # ix.resin_density_app.fix(700)
            ix.resin_density.fix(700)
            ix.ebct.fix(140)
            ix.bed_porosity.fix(0.4)
            ix.loading_rate.fix(0.003)
            ix.resin_porosity.fix(0.4)
            ix.resin_diam.fix(500e-6)
# 				
            ix.b0.set_value(0.86989825)
            ix.b1.set_value(0.15652075)
            ix.b2.set_value(0.53451525)
            ix.b3.set_value(0.00182175)
            ix.b4.set_value(0.15304425)
            ix.a0.set_value(6.07)
            ix.a1.set_value(21.31)
            ix.add_ss_approximation()
            self.deactivate_ss_approximation(m=m)
# 6.07E+00	2.13E+01
            # ix.process_flow.mass_transfer_term[(0.0, "Liq", self.target_component)].fix(1e-15)
        for k, v in self.fix_vars_dict.items():
            print(f"\nfixing {k} to {v} for design...")
            ixv = getattr(ix, k)
            ixv.fix(v)
        return m

    def build_design(self, m=None):

        if m is None:
            m = self.m

        self.ix = ix = m.fs.ix
        pf = ix.process_flow
        # m.fs.properties.set_default_scaling(
        #     "flow_mol_phase_comp",
        #     self.c0_mol_flow_sf,
        #     # 1e15,
        #     index=("Liq", self.target_component),
        # )

        # m.fs.properties.set_default_scaling(
        #     "flow_mol_phase_comp",
        #     # 1e5,
        #     self.water_mol_flow_sf,
        #     index=("Liq", "H2O"),
        # )
        pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].unfix()
        pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].unfix()

        pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].set_value(
            self.water_mol_flow
        )
        pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].set_value(
            self.c0_mol_flow
        )
        set_scaling_factor(
            pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component],
            self.c0_mol_flow_sf,
        )
        set_scaling_factor(
            pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"],
            self.water_mol_flow_sf,
        )
        pf.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): self.flow_in,
                ("conc_mass_phase_comp", ("Liq", self.target_component)): self.c0,
            },
            hold_state=True,
        )

        ix.loading_rate.fix(self.loading_rate)
        ix.ebct.fix(self.ebct)
        ix.resin_porosity.fix(self.resin_porosity)
        ix.resin_diam.fix(self.resin_diam)
        # ix.resin_density_app.fix(self.resin_density)
        ix.resin_density.fix(self.resin_density)

        self._fix_initial_guess(m=m)
        self._calc_from_constr(m=m)

        print(f"dof for design build = {degrees_of_freedom(m)}")

        for k, v in self.fix_vars_dict.items():
            print(f"\nfixing {k} to {v} for design...")
            ixv = getattr(ix, k)
            ixv.fix(v)

    def build_parmest(self):

        self.m_pe = m = self.build_init()
        ix = m.fs.ix
        self._set_bounds(m=self.m_pe)
        self._calc_from_constr(m=self.m_pe)
        self.scale_it(m=self.m_pe)
        ix.initialize()
        self.activate_surrogate(m=self.m_pe)
        # self.m_pe.fs.ix.N_Bi.display()
        self.solve_it(m=self.m_pe)
        self.build_design(m=self.m_pe)
        self.solve_it(m=self.m_pe)
        set_scaling_factor(self.m_pe.fs.ix.bv, 1e-5)

        return m

    def deactivate_ss_approximation(self, m=None):
        if m is None:
            m = self.m
        ix = m.fs.ix
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
        ix.c_norm_avg.fix(0.1)
        ix.eq_mass_transfer_term.deactivate()
        ix.process_flow.mass_transfer_term[0, "Liq", self.target_component].fix(0)

    def scale_it(self, m=None):
        """
        Scale the model
        """
        if m is None:
            m = self.m

        ix = m.fs.ix
        for c in ix.component_objects(Var):
            unset_scaling_factor(c)

        set_scaling_factor(
            ix.process_flow.properties_in[0].flow_mol_phase_comp[
                "Liq", self.target_component
            ],
            self.c0_mol_flow_sf,
        )
        set_scaling_factor(
            ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"],
            self.water_mol_flow_sf,
        )
        set_scaling_factor(ix.bed_depth, 1e2)
        set_scaling_factor(ix.bed_diameter, 1e3)
        set_scaling_factor(ix.bed_area, 1e5)
        set_scaling_factor(ix.bed_volume, 1e6)
        set_scaling_factor(ix.c_eq, 1e5)

        calculate_scaling_factors(m)

        for v, sf in self.scale_var_dict.items():
            ixv = getattr(ix, v)
            if ixv.is_indexed():
                for _, ixvi in ixv.items():
                    set_scaling_factor(ixvi, sf)
            else:
                set_scaling_factor(ixv, sf)

        for c, sf in self.scale_constr_dict.items():
            ixc = getattr(ix, c)
            if ixc.is_indexed():
                for _, ixci in ixc.items():
                    constraint_scaling_transform(ixci, sf)
            else:
                constraint_scaling_transform(ixc, sf)

        if self.ix_model is IonExchangeClark:
            constraint_scaling_transform(ix.eq_clark[self.target_component], 1e-2)

        if self.ix_model is IonExchangeThomas:
            constraint_scaling_transform(ix.eq_thomas[self.target_component], 1e-2)

        # constraint_scaling_transform(ix.eq_gnielinski[self.target_component], 1e-4)
        # constraint_scaling_transform(
        #     ix.eq_solute_dist_param[self.target_component], 1e-2
        # )
        # constraint_scaling_transform(ix.eq_minimum_breakthrough_time, 1e-4)
        # constraint_scaling_transform(ix.eq_Bi, 1e7)
        # constraint_scaling_transform(ix.eq_Sc[self.target_component], 1e-4)

    def solve_it(self, m=None, optarg=dict()):
        """
        Solve the model
        """

        if m is None:
            m = self.m
        for k, v in optarg.items():
            self.solver.options[k] = v
        try:
            # Should never break a model run if an error is thrown
            self.results = self.solver.solve(m, tee=self.tee)
            self.tc = self.results.solver.termination_condition
        except:
            print_infeasible_constraints(m)
            print_variables_close_to_bounds(m)
            self.tc = "fail"

    def build_parmest_experiment_list(self):
        self.experiment_list = list()
        exp_dict = {
            "exp_num": [x for x in range(len(self.keep_bvs))],
            "cb": self.keep_cbs,
            "bv": self.keep_bvs,
            "c_norm": self.keep_cnorms,
        }
        self.df_exp = pd.DataFrame.from_dict(exp_dict).set_index("exp_num")
        self.df_exp["c0"] = self.c0
        # self.df_exp = self.filtered_data[self.filtered_data.bv.isin(self.keep_bvs)].copy()
        # self.df_exp.reset_index(inplace=True, drop=True)
        # for i in self.filtered_data.index:
        for i in self.df_exp.index:
            self.experiment_list.append(IXExperiment(self, self.df_exp, i))

    def run_parmest(self):
        """
        Run parmest
        """

        print(f"\nRunning parmest...")

        def SSE(model):
            """
            Sum of squared error between `experiment_output` model and data values
            """
            expr = sum(
                (y - y_hat) ** 2 for y, y_hat in model.experiment_outputs.items()
            )
            return expr * self.expr_sf

        self.build_parmest_experiment_list()
        self.pestimator = parmest.Estimator(
            self.experiment_list, obj_function=SSE, **self.parmest_kwargs
        )
        self.obj, self.parmest_theta = self.pestimator.theta_est()
        self.theta_dict = dict()
        for t in self.theta_names:
            self.theta_dict[t] = self.parmest_theta[f"fs.ix.{t}"]

    def test_theta(
        self,
        m=None,
        min_test_cbs=None,
        max_test_cbs=0.99,
        min_test_cnorms=0.05,
        max_test_cnorms=0.95,
        min_test_bvs=None,
        max_test_bvs=None,
        num_pts=9,
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
        # bv_pred_theta_test = Predicted BV using theta and regularly spaced test data
        self.bv_pred_theta_test = []
        self.bv_fail_theta_test = []

        # cb_pred_theta_test = Predicted Cb using theta and regularly spaced test data
        self.cnorm_pred_theta_test = []
        self.cnorm_fail_theta_test = []
        if m is None:
            self.rebuild_test_theta()
            m = self.m
        # self.ix = ix = self._get_ix_blk(m=self.m_theta)
        # ix = m.fs.ix
        # ix.initialize()
        # self.activate_surrogate(m=m)
        # self.solve_it(m=m)
        # self.build_design(m=m)
        # self.solve_it(m=m)
        # for k, v in self.theta_dict.items():
        #     ixv = getattr(m.fs.ix, k)
        #     ixv.fix(v)
        # m0 = m.clone()
        ix = m.fs.ix

        # assert degrees_of_freedom(self.m_theta) == 0

        self.flag = "test_theta_fake"

        for bv, cnorm in zip(self.test_bvs, self.test_cnorms):

            ix.c_norm[self.target_component].fix(cnorm)
            ix.bv.set_value(bv)
            assert degrees_of_freedom(m) == 0
            self.solve_it(m=m)
            if self.tc != "optimal":
                self.bv_fail_theta_test.append(bv)
                self.cnorm_fail_theta_test.append(cnorm)
            elif self.tc == "optimal":
                self.bv_pred_theta_test.append(ix.bv())
                self.cnorm_pred_theta_test.append(cnorm)
            self.results_dict_append(ix_blk=ix)

        self.bv_pred_theta = []
        self.bv_fail_theta = []
        self.keep_bv_theta = []
        self.cnorm_pred_theta = []
        self.cnorm_fail_theta = []

        self.flag = "test_theta_real"
        # m = self.m_theta.clone()
        # ix = m.fs.ix
        if m is None:
            self.rebuild_test_theta()
            m = self.m
        ix = m.fs.ix
        for bv, cnorm in zip(self.keep_bvs[::-1], self.keep_cnorms[::-1]):
            if cnorm == 0:
                cnorm = 1e-3

            ix.c_norm[self.target_component].fix(cnorm)
            ix.bv.set_value(bv)

            assert degrees_of_freedom(m) == 0
            self.solve_it(m=m)
            if self.tc != "optimal":
                self.bv_fail_theta.append(bv)
                self.cnorm_fail_theta.append(cnorm)
            elif self.tc == "optimal":
                self.bv_pred_theta.append(ix.bv())
                self.keep_bv_theta.append(bv)
                self.cnorm_pred_theta.append(cnorm)
            self.results_dict_append(ix_blk=ix)

        self.bv_error = [
            pred - actual
            for (pred, actual) in zip(self.bv_pred_theta, self.keep_bv_theta)
        ]
        self.bv_rel_error = [
            (pred - actual) / actual
            for (pred, actual) in zip(self.bv_pred_theta, self.keep_bv_theta)
        ]
        self.bv_abs_rel_error = [abs(e) for e in self.bv_rel_error]
        self.bv_abs_error = [abs(e) for e in self.bv_error]
        self.mean_absolute_error = np.mean(self.bv_abs_error)
        self.sum_squared_error = sum(
            [
                (pred - actual) ** 2
                for (pred, actual) in zip(self.bv_pred_theta, self.keep_bv_theta)
            ]
        )
        self.corr_matrix = np.corrcoef(self.keep_bv_theta, self.bv_pred_theta)
        corr = self.corr_matrix[0, 1]
        self.R_sq = corr**2

    def build_results_dict(self, **kwargs):
        build_results_dict(self, **kwargs)

    def results_dict_append(self, **kwargs):
        results_dict_append(self, **kwargs)

    def make_output_df(self):
        make_output_df(self) 

    def save_output(self):
        save_output(self)

    def make_results_df(self):
        make_results_df(self)

    def save_results(self):
        save_results(self)

    def save_figs(self):
        save_figs(self)

    def plot_initial_guess(self):
        plot_initial_guess(self)

    def plot_theta(self):
        plot_theta(self)

    def plot_estimate_bv50(self):
        plot_estimate_bv50(self)

    def plot_curve(self):
        plot_curve(self)

    def _get_comp_list(self, blk, comp=Var, skip_list=[]):
        cs = []
        split_name = blk.name + "."
        skip_list += ["ref", "process_flow", "regeneration"]
        for c in blk.component_objects(comp):
            if any(s in c.name for s in skip_list):
                continue
            # cs.append(c.name.split(split_name)[1])
            yield c.name.split(split_name)[1]
        # return cs

    def _set_bounds(self, m=None):

        if m is None:
            m = self.m
        # else:
        #     ix = self._get_ix_blk(m=m)
        ix = m.fs.ix
        if self.set_bounds_dict == "all":
            for v in ix.component_objects(Var):
                # print(v.name)
                v.setlb(None)
                v.setub(None)
                # v.display()
                # if "N_Bi" in v.name:
                #     v.display()
        else:
            for k, v in self.set_bounds_dict.items():
                ixv = getattr(ix, k)
                ixv.setlb(v[0])
                ixv.setub(v[1])

    def _fix_initial_guess(self, m=None):

        if m is None:
            ix = self.ix
        else:
            ix = self._get_ix_blk(m=m)

        for k, v in self.initial_guess_dict.items():
            ixv = getattr(ix, k)
            ixv.fix(v)

    def _calc_from_constr(self, m=None):

        if m is None:
            m = self.m
        # else:
        #     ix = self._get_ix_blk(m=m)
        ix = m.fs.ix
        # print("calc_from_constr")

        for v, c in self.calc_from_constr_dict.items():
            # print("\t", v, c)

            ixv = getattr(ix, v)
            ixc = getattr(ix, c)
            if ixc.is_indexed():
                for i, ic in ixc.items():
                    if ixv.is_indexed():
                        calculate_variable_from_constraint(ixv[i], ic)
                    else:
                        calculate_variable_from_constraint(ixv, ic)

            else:
                calculate_variable_from_constraint(ixv, ixc)

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


class IXExperiment(Experiment):
    def __init__(self, ix_parmest_obj, data, experiment_number):
        self.data = data
        self.experiment_number = experiment_number
        self.data_i = data.loc[experiment_number]
        self.ix_parmest_obj = ix_parmest_obj
        self.model = None
        self.thetas = self.ix_parmest_obj.thetas

    def create_model(self):
        self.model = self.ix_parmest_obj.build_parmest()

    def finalize_model(self):

        model = self.model
        cb = float(self.data_i.cb)
        c0 = float(self.data_i.c0)
        bv = float(self.data_i.bv)
        self.c_norm = c_norm = cb / c0
        model.fs.ix.c_norm.fix(c_norm)
        model.fs.ix.bv.set_value(bv)
        print(f"initializing parmest model for {c_norm}...\n")
        self.ix_parmest_obj.solve_it(m=self.model)
        print(f"parmest build termination {self.ix_parmest_obj.tc}")
        for k in self.ix_parmest_obj.initial_guess_dict.keys():
            ixv = getattr(self.model.fs.ix, k)
            ixv.unfix()
        print(f"\ndof = {degrees_of_freedom(self.model)}")
        # assert len(self.thetas) == degrees_of_freedom(self.model)
        # self.ix_parmest_obj.scale_it(m=self.model)

        # try:
        #     self.model.fs.ix.initialize()
        # except:
        #     pass

        return model

    def label_model(self):
        m = self.model

        m.experiment_outputs = Suffix(direction=Suffix.LOCAL)
        m.experiment_outputs.update(
            [
                (
                    m.fs.ix.c_norm[m.fs.ix.config.target_component],
                    self.data_i["c_norm"],
                ),
                (m.fs.ix.bv, self.data_i["bv"]),
            ]
        )
        m.unknown_parameters = Suffix(direction=Suffix.LOCAL)
        m.unknown_parameters.update(
            [
                (k, ComponentUID(k))
                for k in [getattr(m.fs.ix, theta) for theta in self.thetas]
            ]
        )
        return m

    def get_labeled_model(self):
        m = self.create_model()
        m = self.finalize_model()
        m = self.label_model()

        return m


##### OLD SCALING ROUTINE


# set_scaling_factor(
#     prop_out.flow_mol_phase_comp["Liq", target_component], self.c0_mol_flow_sf
# )
# set_scaling_factor(
#     prop_out.flow_mol_phase_comp["Liq", "H2O"], self.water_mol_flow_sf
# )

# set_variable_scaling_from_current_value(m, descend_into=True, overwrite=True)

# if self.autoscale_fixed:
#     for v in self._get_comp_list(ix, skip_list=["traps"]):
#         ixv = getattr(ix, v)
#         if ixv.is_indexed():
#             for i, vv in ixv.items():
#                 if vv.is_fixed():
#                     if value(vv) == 0:
#                         continue
#                     sf = 1 / value(vv)
#                     set_scaling_factor(vv, sf)
#         else:
#             if ixv.is_fixed():
#                 if value(ixv) == 0:
#                     continue
#                 sf = 1 / value(ixv)
#                 set_scaling_factor(ixv, sf)

# if isinstance(self.scale_from_value, list):
#     for v in self.scale_from_value:
#         # print(v)
#         ixv = getattr(ix, v)
#         if ixv.is_indexed():
#             for i, vv in ixv.items():
#                 if vv.is_fixed():
#                     if value(vv) == 0:
#                         continue
#                     sf = 1 / value(vv)
#                     if get_scaling_factor(vv) is None:
#                         set_scaling_factor(vv, sf)
#         else:
#             if value(ixv) == 0:
#                 continue
#             sf = 1 / value(ixv)
#             if get_scaling_factor(ixv) is None:
#                 set_scaling_factor(ixv, sf)
# if self.ix_model is IonExchangeCPHSDM:
#     bed_area = self.bed_volume() / self.bed_depth()
#     vel_inter = self.loading_rate / 0.4
#     t_contact = self.ebct * 0.4
#     calculate_scaling_factors(m)
# set_variable_scaling_from_current_value(ix.process_flow, descend_into=True, overwrite=True)
#     set_scaling_factor(ix.bed_volume, 1 / self.bed_volume())
#     # set_scaling_factor(ix.breakthrough_time, 1e-6)
#     # set_scaling_factor(ix.bv, 1e-5)
# set_scaling_factor(ix.c_eq, 1e5)
#     # set_scaling_factor(ix.min_breakthrough_time, 1e-6)
#     # set_scaling_factor(ix.spdfr, 1)
#     # set_scaling_factor(ix.bed_area, 1 / bed_area)
#     # set_scaling_factor(ix.bed_diameter, 1 / self.bed_diameter())
#     # set_scaling_factor(ix.solute_dist_param, 1e-4)
#     # set_scaling_factor(ix.throughput, 1e3)
#     # set_scaling_factor(ix.vel_inter, 1 / vel_inter)
#     # set_scaling_factor(ix.t_contact, 1 / t_contact)
#     # set_scaling_factor(ix.resin_density_app, 1 / self.resin_density)
#     # set_scaling_factor(ix.c_norm, 1 / self.c_norm)
#     # set_scaling_factor(ix.Bi, 0.1)
