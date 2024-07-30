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
    save_output,
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
    # surr_dir + "/throughput_pysmo_surr_linear.json",
    "/Users/ksitterl/Documents/Python/nawi-analysis/NAWI-analysis/analysis_waterTAP/analysisWaterTAP/analysis_scripts/pfas_ix_cphsdm/surrogates/throughput_pysmo_surr_linear.json"
)
mw_water = 0.018 * pyunits.kg / pyunits.mol
class IXParmest:
    def __init__(
        self,
        curve_id=None,
        input_data=None,
        ix_model=IonExchangeClark,
        build_func=None,
        build_func_kwargs=dict(),
        regenerant="single_use",
        diff_calculation="HaydukLaudie",
        data_file=None,
        resin_data=dict(),
        compound_data=dict(),
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
        rho=1000, #kg/m3
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
        self.ix_model = ix_model
        self.regenerant = regenerant
        self.fix_vars_dict = fix_vars_dict

        self.figsize = figsize
        self.bv = bv  # need value for initial model build
        self.cb = cb  # need value for initial model build
        self.c_norm = c_norm
        self.rho = rho * pyunits.kg / pyunits.m**3

        self.resin_data = resin_data
        self.compound_data = compound_data
        self.parmest_kwargs = parmest_kwargs
        self.just_plot_curve = just_plot_curve
        self.use_all_data = use_all_data
        self.use_this_data = use_this_data
        self.diff_calculation = diff_calculation
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
            return

        self.get_model_config()

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


    def ix_cphsdm_init(self, m=None):

        if m is None:
            m = self.m
        assert self.ix_model is IonExchangeCPHSDM
        self.activate_surrogate(m=m)

        # return m 
        # try:
        self.solve_it(m=m)
        print(f"termination with surrogates {self.tc}")
        # except:
        print_infeasible_constraints(m)



    def rebuild(self):
        self.m0 = self.build_init()

        self.m = self.m0.clone()
        # self.m.fs.ix.initialize()

        # self.m = self.build_func(**self.build_func_kwargs)
        # self.ix = self._get_ix_blk()
        self._set_bounds()
        # self._fix_initial_guess()
        self._calc_from_constr()
        self.scale_it()
        # print(f"dof = {degrees_of_freedom(self.m)}")
        # # assert degrees_of_freedom(self.m) == 0
        # try:
        #     self.ix.initialize()
        # # except InitializationError:
        # except:
        #     print("Initial build of model failed.\nTry a different initial guess.")
        #     print_infeasible_constraints(self.m)
        #     print_variables_close_to_bounds(self.m)
        
        # if self.ix_model is IonExchangeCPHSDM:
        #     self.ix_cphsdm_init()


    def test_initial_guess(self, min_cnorm_ig=0.05, max_cnorm_ig=0.95, num_pts=20):
        """
        Loops through the kept (bv, cb) pairs and runs the model
        with those values using initial_guess_dict
        """
        print("\nTesting initial guess\n")

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
        m = self.m.clone()
        ix = m.fs.ix
        for bv, cnorm in zip(self.keep_bvs[::-1], self.keep_cnorms[::-1]):
            print(f"testing initial guess with actual data:\n\tBV = {bv}\n\tC/C0 = {cnorm}...")
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
        self.c0 = df.c0.iloc[0]

        self.compound = df.compound.iloc[0]
        self.flow_in = float(df.flow_in.iloc[0])  # m3/s
        self.loading_rate = float(df.loading_rate.iloc[0])  # m/s
        self.bed_depth = float(df.bed_depth.iloc[0]) * pyunits.m
        self.bed_diameter = float(df.bed_diam.iloc[0]) * pyunits.m
        self.bed_volume_total = float(df.bed_vol.iloc[0]) * pyunits.m**3
        self.resin = df.resin.iloc[0]
        self.ebct_min = float(df.ebct.iloc[0])  # minutes
        self.ebct = pyunits.convert(
            self.ebct_min * pyunits.min, to_units=pyunits.second
        )
        self.conc_units_str = df.conc_units.iloc[0]  # conc_mass_phase_comp units
        numer = getattr(pyunits, self.conc_units_str.split("/")[0])
        denom = getattr(pyunits, self.conc_units_str.split("/")[1])
        self.conc_units = numer / denom
        self.c0_kg_m3 = pyunits.convert(self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3)

        self.target_component = self.compound_data["name"]
        self.charge = self.compound_data["charge"]
        self.mw = self.compound_data["mw_comp"] * pyunits.kg / pyunits.mol
        if "molar_vol" in self.compound_data.keys():
            self.molar_vol = self.compound_data["molar_vol"]
        self.diffusivity = self.compound_data["diffusivity"]
        self.filter_data()

        try:
            if ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2) > 0:
                self.expr_sf = 1 / ((max(self.keep_bvs) - min(self.keep_bvs)) ** 2)
        except:
            self.expr_sf = 1e-9

        self.resin_density = self.resin_data["density"]
        self.resin_diam = self.resin_data["diameter"]
        self.bed_porosity = self.resin_data["porosity"]

    def get_model_config(self):

        mw_water = 0.018 * pyunits.kg / pyunits.mol
        rho = 1000 * pyunits.kg / pyunits.m**3
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
            self.keep_bvs = tmp.bv.to_list()
            self.keep_cbs = tmp.cb.to_list()
            # assert 0 not in self.keep_cbs
        elif self.use_this_data is not None:
            print("use this data")
            tmp = df.reset_index(drop=True)
            tmp = tmp.loc[self.use_this_data].copy()
            self.keep_bvs = tmp.bv.to_list()
            self.keep_cbs = tmp.cb.to_list()
            # assert 0 not in self.keep_cbs

        else:
            print("default filtering")

            for i, cb in enumerate(self.all_cbs):
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

        self.keep_cnorms = [cb / self.c0 for cb in self.keep_cbs]
        self.excl_cnorms = [cb / self.c0 for cb in self.excl_cbs]
    
    
    def build_init(self):
        """
        Method used to build IX model
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**self.ion_props)
        # m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
        # m.fs.properties.dens_mass_const = 1000

        self.ix_config["property_package"] = m.fs.properties
        self.flow_in_init = 2e-8
        self.c0_init = 4e-9
        self.mw_init = 0.5 * pyunits.kg / pyunits.mol

        self.c0_mol_flow_init = pyunits.convert(
            (self.c0_init * pyunits.kg / pyunits.m**3)
            / self.mw_init
            * (self.flow_in_init * pyunits.m**3 / pyunits.s),
            to_units=pyunits.mol / pyunits.s,
        )
        self.c0_mol_flow_sf_init = 1 / value(self.c0_mol_flow)
        self.water_mol_flow_init = pyunits.convert(
            ((self.flow_in * pyunits.m**3 / pyunits.s) * self.rho) / mw_water,
            to_units=pyunits.mol / pyunits.s,
        )
        self.water_mol_flow_sf_init = 1 / self.water_mol_flow()

        m.fs.ix = ix = self.ix_model(**self.ix_config)

        self._set_bounds(m=m)

        pf = ix.process_flow

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp",
            self.c0_mol_flow_sf_init,
            index=("Liq", self.target_component),
        )

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", self.water_mol_flow_sf_init, index=("Liq", "H2O")
        )

        # c_in = pyunits.convert(
        #     self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3
        # )
        # c_b = pyunits.convert(self.cb, to_units=pyunits.kg / pyunits.m**3)

        # if self.c_norm is None:
        #     self.c_norm = value(c_b / c_in)

        pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(self.water_mol_flow_init)
        pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].fix(
            self.c0_mol_flow_init
        )
        pf.properties_in[0].pressure.fix(101325)
        pf.properties_in[0].temperature.fix(298)
        pf.properties_in[0].flow_vol_phase[...]
        pf.properties_out[0].flow_vol_phase[...]

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

            for k, v in self.fix_vars_dict.items():
                ixv = getattr(ix, k)
                ixv.fix(v)

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

            ix.eq_Pe_bed.deactivate()
            ix.eq_Pe_p.deactivate()
            ix.eq_service_flow_rate.deactivate()
            ix.eq_column_height.deactivate()
            ix.eq_bed_design.deactivate()

            ix.number_columns.fix(1)
            ix.freundlich_ninv.fix(0.9)
            ix.freundlich_k.fix(1)
            ix.surf_diff_coeff.fix(8e-15)
            ix.c_norm.fix(0.5)

            ix.shape_correction_factor.fix(1)
            ix.tortuosity.fix(1)
            ix.shape_correction_factor.fix(1)
            ix.resin_density_app.fix(700)
            ix.ebct.fix(140)
            # ix.ebct.fix(10)
            ix.bed_porosity.fix(0.4)
            ix.loading_rate.fix(0.003)
            ix.resin_porosity.fix(0.4)
            # ix.bed_volume = 1e-6
            ix.resin_diam.fix(500e-6)
            # ix.loading_rate.fix(0.001)
            # ix.bed_depth.fix(0.001)

            # ix.process_flow.mass_transfer_term[(0.0, "Liq", self.target_component)].fix(1e-15)
        pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)


        return m

    def build_parmest(self):
            pass
        

        # if self.ix_model is IonExchangeCPHSDM:

        #     ix.eq_Pe_bed.deactivate()
        #     ix.eq_Pe_p.deactivate()
        #     ix.eq_service_flow_rate.deactivate()
        #     ix.eq_column_height.deactivate()
        #     # ix.eq_Bi.deactivate()
        #     ix.eq_bed_design.deactivate()
        #     # ix.eq_surf_diff_coeff.deactivate()

            
# initial_guess_dict = dict(
#     freundlich_k=1,
#     freundlich_ninv=1,
#     surf_diff_coeff=8e-15
# )
# fix_vars_dict = dict(
#     bed_porosity=0.4,
#     ebct=10,
#     resin_density_app=700,

# )
            # ix.freundlich_ninv.fix(1)
            # ix.number_columns.fix(1)

            # ix.freundlich_ninv.fix(0.9)
            # ix.freundlich_k.fix(1)
            # ix.surf_diff_coeff.fix(8e-15)

            # ix.shape_correction_factor.fix(1)
            # ix.tortuosity.fix(1)
            # ix.shape_correction_factor.fix(1)
            # ix.resin_density_app.fix(700)
            # ix.ebct.fix(10)
            # ix.bed_porosity.fix(0.4)



            # ix.resin_density_app.fix(self.resin_density)
            # ix.resin_diam.fix(self.resin_diam)
            # ix.ebct.fix(self.ebct)
            # ix.bed_porosity.fix(0.4)
            # ix.loading_rate.fix(self.loading_rate)
            # ix.c_norm.fix(self.c_norm)
            # ix.resin_porosity.fix(self.resin_data["particle_porosity"])
            # pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)
            # ix.initialize()
            # activate_surrogate(m)

            # self.solve_it(m=m)

            # print(f"initial solve: {self.tc}")

            # activate_surrogate(m)
            # self.solve_it(m=m)
            # print(f"solve with surrogates: {self.tc}")
            # ix.a0.set_value(3.12)
            # ix.a1.set_value(12.8)
            # ix.b0.set_value(0.8789)
            # ix.b1.set_value(0.15652)
            # ix.b2.set_value(0.5345)
            # ix.b3.set_value(0.00182175)
            # ix.b4.set_value(0.153044)



            # ix.process_flow.mass_transfer_term[(0.0, "Liq", self.target_component)].fix(1e-15)
        # pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)

        # for var, val in self.fix_vars_dict.items():
        #     ixv = getattr(ix, var)
        #     ixv.fix(val)
        # return m
    


    def scale_it(self, m=None):
        """
        Scale the model
        """
        if m is None:
            m = self.m
            # ix = self.ix
        # else:
        #     ix = self._get_ix_blk(m=m)
        ix = m.fs.ix


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
                ixv = getattr(ix, v)
                if ixv.is_indexed():
                    for i, vv in ixv.items():
                        if vv.is_fixed():
                            if value(vv) == 0:
                                continue
                            sf = 1 / value(vv)
                            set_scaling_factor(vv, sf)
                else:
                    if ixv.is_fixed():
                        if value(ixv) == 0:
                            continue
                        sf = 1 / value(ixv)
                        set_scaling_factor(ixv, sf)

        if isinstance(self.scale_from_value, list):
            for v in self.scale_from_value:
                # print(v)
                ixv = getattr(ix, v)
                if ixv.is_indexed():
                    for i, vv in ixv.items():
                        if vv.is_fixed():
                            if value(vv) == 0:
                                continue
                            sf = 1 / value(vv)
                            if get_scaling_factor(vv) is None:
                                set_scaling_factor(vv, sf)
                else:
                    if value(ixv) == 0:
                        continue
                    sf = 1 / value(ixv)
                    if get_scaling_factor(ixv) is None:
                        set_scaling_factor(ixv, sf)

        if self.ix_model is IonExchangeClark:
            constraint_scaling_transform(ix.eq_clark[target_component], 1e-2)

        if self.ix_model is IonExchangeThomas:
            constraint_scaling_transform(ix.eq_thomas[target_component], 1e-2)
        
        calculate_scaling_factors(m)
        set_scaling_factor(ix.bed_volume, 1e4)
        set_scaling_factor(ix.breakthrough_time, 1e-4)
        set_scaling_factor(ix.bv, 1e-5)
        set_scaling_factor(ix.c_eq, 1e6)
        set_scaling_factor(ix.min_breakthrough_time, 1e-4)
        set_scaling_factor(ix.spdfr, 1e2)
        set_scaling_factor(ix.bed_area, 1e5)
        set_scaling_factor(ix.solute_dist_param, 1e-4)
        set_scaling_factor(ix.throughput, 1e-3)
        set_scaling_factor(ix.vel_inter, 1e2)
        set_scaling_factor(ix.t_contact, 1e2)
        set_scaling_factor(ix.resin_density, 1 / self.resin_density)
        set_scaling_factor(ix.c_norm, 10)
        set_scaling_factor(ix.Bi, 0.1)



    def solve_it(self, m=None, solver="watertap", optarg=dict(), tee=False):
        """
        Solve the model
        Defaults to using watertap solver
        """

        if m is None:
            m = self.m

        # if solver != "watertap":
        #     self.solver = SolverFactory("ipopt")
        # else:
        #     self.solver = get_solver()
        for k, v in optarg.items():
            self.solver.options[k] = v
        try:
            # Should never break a model run if an error is thrown
            self.results = self.solver.solve(m, tee=tee)
            self.tc = self.results.solver.termination_condition
        except:
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

        self.build_parmest_experiment_list()
        self.pestimator = parmest.Estimator(
            self.experiment_list, obj_function="SSE", **self.parmest_kwargs
        )
        self.obj, self.parmest_theta = self.pestimator.theta_est()
        self.theta_dict = dict()
        for t in self.theta_names:
            self.theta_dict[t] = self.parmest_theta[f"fs.ix.{t}"]

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
        # bv_pred_theta_test = Predicted BV using theta and regularly spaced test data
        self.bv_pred_theta_test = []
        self.bv_fail_theta_test = []

        # cb_pred_theta_test = Predicted Cb using theta and regularly spaced test data
        self.cnorm_pred_theta_test = []
        self.cnorm_fail_theta_test = []

        self.rebuild()
        self.m_theta = self.m.clone()
        self.ix = ix = self._get_ix_blk(m=self.m_theta)
        for k, v in self.theta_dict.items():
            ixv = getattr(ix, k)
            ixv.fix(v)

        self._calc_from_constr(m=self.m_theta)
        self.scale_it(m=self.m_theta)
        # activate_surrogate(self.m_theta)
        # self.ix_cphsdm_init(m=self.m_theta)

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
            self.results_dict_append(ix_blk=ix)

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
            self.results_dict_append(ix_blk=ix)

        self.corr_matrix = np.corrcoef(self.keep_bv_theta, self.bv_pred_theta)
        corr = self.corr_matrix[0, 1]
        self.R_sq = corr**2

    def build_results_dict(self, **kwargs):
        build_results_dict(self, **kwargs)

    def results_dict_append(self, **kwargs):
        results_dict_append(self, **kwargs)

    def save_output(self):
        save_output(self)

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
                v.setlb(None)
                v.setub(None)
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

        for v, c in self.calc_from_constr_dict.items():
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
        # if self.ix_parmest_obj.ix_model is IonExchangeClark:
        #     self.thetas = ["freundlich_n", "mass_transfer_coeff", "bv_50"]
        # if self.ix_parmest_obj.ix_model is IonExchangeThomas:
        #     self.thetas = ["freundlich_n", "mass_transfer_coeff", "bv_50"]

    def create_model(self):
        self.model = self.ix_parmest_obj.build_it()

    def finalize_model(self):

        model = self.model
        cb = value(
            pyunits.convert(
                float(self.data_i.cb) * self.ix_parmest_obj.conc_units,
                to_units=pyunits.kg / pyunits.m**3,
            )
        )

        c0 = value(
            pyunits.convert(
                float(self.data_i.c0) * self.ix_parmest_obj.conc_units,
                to_units=pyunits.kg / pyunits.m**3,
            )
        )
        bv = float(self.data_i.bv)
        self.c_norm = c_norm = cb / c0
        model.fs.ix.c_norm.fix(c_norm)
        model.fs.ix.bv.set_value(bv)
        self.ix_parmest_obj._set_bounds(m=self.model)
        self.ix_parmest_obj._fix_initial_guess(m=self.model)
        self.ix_parmest_obj._calc_from_constr(m=self.model)
        # self.ix_parmest_obj.scale_it(m=self.model)
        print(f"initializing parmest model for {c_norm}...\n")
        self.model.fs.ix.initialize()
        self.ix_parmest_obj.ix_cphsdm_init(m=self.model)
        for k in self.ix_parmest_obj.initial_guess_dict.keys():
            ixv = getattr(self.model.fs.ix, k)
            ixv.unfix()
        print(f"\ndof = {degrees_of_freedom(self.model)}")
        # assert len(self.thetas) == degrees_of_freedom(self.model)
        self.ix_parmest_obj.scale_it(m=self.model)

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

    # def estimate_bv50(self):
    #     def linear_fit(bv, slope, b):
    #         return slope * bv + b

    #     def linear_inv(cb, slope, b):
    #         return (cb - b) / slope

    #     if not self.ix_model is IonExchangeClark:
    #         raise ValueError("Can only use estimate_bv50 method if using Clark model.")

    #     if not hasattr(self, "bv_pred_ig"):
    #         self.test_initial_guess()

    #     self.bv50_bvs = [
    #         bv
    #         for bv, cnorm in zip(self.keep_bvs, self.keep_cnorms)
    #         if cnorm >= self.cb50_min_thresh
    #     ]

    #     self.bv50_cbs = [
    #         cnorm for cnorm in self.keep_cnorms if cnorm >= self.cb50_min_thresh
    #     ]

    #     cb_50 = 0.5

    #     self.bv50_fit_params, self.bv50_fit_cov = curve_fit(
    #         linear_fit, self.bv50_bvs, self.bv50_cbs
    #     )

    #     self.bv50_slope, self.bv50_int = [*self.bv50_fit_params]

    #     self.bv50_est = linear_inv(cb_50, self.bv50_slope, self.bv50_int)

    #     self.bv50_orig = self.initial_guess_dict["bv_50"]
    #     self.bv_pred_ig_orig = deepcopy(self.bv_pred_ig)
    #     self.cb_pred_ig_orig = deepcopy(self.cb_pred_ig)

    #     self.bv50_corr_matrix = np.corrcoef(
    #         self.bv50_bvs, [self._linear_fit(bv) for bv in self.bv50_bvs]
    #     )
    #     corr = self.bv50_corr_matrix[0, 1]
    #     self.bv50_R_sq = corr**2
    #     if self.bv50_est > 100:
    #         self.initial_guess_dict["bv_50"] = self.bv50_est
    #     else:
    #         self.initial_guess_dict["bv_50"] = self.bv50_orig
    #     self.initial_guess_theta_dict = deepcopy(self.initial_guess_dict)
    #     # rebuild model with new initial guess

    #     self.rebuild()
    #     self.test_initial_guess()
        # self.plot_estimate_bv50()
# ================ OLD BUILD FUNCTION

#     def build_it(self):
#         """
#         Method used to build IX model
#         """
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(**self.ion_props)

#         self.ix_config["property_package"] = m.fs.properties

#         m.fs.ix = ix = self.ix_model(**self.ix_config)

#         self._set_bounds(m=m)

#         pf = ix.process_flow

#         m.fs.properties.set_default_scaling(
#             "flow_mol_phase_comp", 1, index=("Liq", "H2O")
#         )
#         m.fs.properties.set_default_scaling(
#             "flow_mol_phase_comp",
#             self.c0_mol_flow_sf,
#             index=("Liq", self.target_component),
#         )

#         m.fs.properties.set_default_scaling(
#             "flow_mol_phase_comp", self.water_mol_flow_sf, index=("Liq", "H2O")
#         )

#         c_in = pyunits.convert(
#             self.c0 * self.conc_units, to_units=pyunits.kg / pyunits.m**3
#         )
#         c_b = pyunits.convert(self.cb, to_units=pyunits.kg / pyunits.m**3)

#         if self.c_norm is None:
#             self.c_norm = value(c_b / c_in)

#         pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(self.water_mol_flow)
#         pf.properties_in[0].flow_mol_phase_comp["Liq", self.target_component].fix(
#             self.c0_mol_flow
#         )
#         pf.properties_in[0].pressure.fix(101325)
#         pf.properties_in[0].temperature.fix(298)
#         pf.properties_in[0].flow_vol_phase[...]
#         pf.properties_out[0].flow_vol_phase[...]

#         if self.ix_model is IonExchangeClark:
#             # Parameters set to force bed_depth == col_height
#             ix.underdrain_h.set_value(0)
#             ix.distributor_h.set_value(0)
#             if self.regenerant != "single_use":
#                 ix.bed_expansion_frac_A.set_value(0)
#                 ix.bw_rate.set_value(0)
#             ix.number_columns.fix(1)

#             # ix.c_breakthru[self.target_component].fix(value(c_b))
#             ix.c_norm[self.target_component].fix(self.c_norm)


#             ix.resin_density.fix(self.resin_density)
#             ix.resin_diam.fix(self.resin_diam)
#             ix.bed_porosity.fix(self.bed_porosity)
#             ix.bed_depth.fix(self.bed_depth)
#             ix.bed_diameter.fix(self.bed_diameter)

#             # since we want to predict BV, we use .set_value() rather than .fix()
#             ix.bv.set_value(self.bv)
#             ix.loading_rate.set_value(self.loading_rate)
#             ix.bed_volume_total.set_value(self.bed_volume_total)

#             for k, v in self.fix_vars_dict.items():
#                 ixv = getattr(ix, k)
#                 ixv.fix(v)

#             # Since we are fitting data to C/C0,
#             # we don't care about the steady-state effluent concentration, so deactivate them
#             # if self.ix_model is IonExchangeClark:
#             ix.eq_tb_traps.deactivate()
#             ix.eq_c_traps.deactivate()
#             ix.eq_traps.deactivate()
#             ix.eq_c_norm_avg.deactivate()
#             for _, v in ix.tb_traps.items():
#                 v.fix()
#             for _, v in ix.c_traps.items():
#                 v.fix()
#             for _, v in ix.traps.items():
#                 v.fix()
#             ix.c_norm_avg.fix()
#             ix.eq_mass_transfer_term.deactivate()

#         if self.ix_model is IonExchangeCPHSDM:

#             ix.eq_Pe_bed.deactivate()
#             ix.eq_Pe_p.deactivate()
#             ix.eq_service_flow_rate.deactivate()
#             ix.eq_column_height.deactivate()
#             # ix.eq_Bi.deactivate()
#             ix.eq_bed_design.deactivate()
#             # ix.eq_surf_diff_coeff.deactivate()

            
# # initial_guess_dict = dict(
# #     freundlich_k=1,
# #     freundlich_ninv=1,
# #     surf_diff_coeff=8e-15
# # )
# # fix_vars_dict = dict(
# #     bed_porosity=0.4,
# #     ebct=10,
# #     resin_density_app=700,

# # )
#             # ix.freundlich_ninv.fix(1)
#             ix.number_columns.fix(1)

#             ix.freundlich_ninv.fix(0.9)
#             ix.freundlich_k.fix(1)
#             ix.surf_diff_coeff.fix(8e-15)

#             ix.shape_correction_factor.fix(1)
#             ix.tortuosity.fix(1)
#             ix.shape_correction_factor.fix(1)
#             ix.resin_density_app.fix(700)
#             ix.ebct.fix(10)
#             ix.bed_porosity.fix(0.4)



#             # ix.resin_density_app.fix(self.resin_density)
#             # ix.resin_diam.fix(self.resin_diam)
#             # ix.ebct.fix(self.ebct)
#             # ix.bed_porosity.fix(0.4)
#             # ix.loading_rate.fix(self.loading_rate)
#             # ix.c_norm.fix(self.c_norm)
#             # ix.resin_porosity.fix(self.resin_data["particle_porosity"])
#             # pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)
#             # ix.initialize()
#             # activate_surrogate(m)

#             # self.solve_it(m=m)

#             # print(f"initial solve: {self.tc}")

#             # activate_surrogate(m)
#             # self.solve_it(m=m)
#             # print(f"solve with surrogates: {self.tc}")
#             # ix.a0.set_value(3.12)
#             # ix.a1.set_value(12.8)
#             # ix.b0.set_value(0.8789)
#             # ix.b1.set_value(0.15652)
#             # ix.b2.set_value(0.5345)
#             # ix.b3.set_value(0.00182175)
#             # ix.b4.set_value(0.153044)



#             # ix.process_flow.mass_transfer_term[(0.0, "Liq", self.target_component)].fix(1e-15)
#         pf.mass_transfer_term[0, "Liq", self.target_component].fix(0)

#         for var, val in self.fix_vars_dict.items():
#             ixv = getattr(ix, var)
#             ixv.fix(val)
#         return m

### ============ OLD PARMEST ROUTINE
### Deprecated for Pyomo 6.7.2

# self.theta_names = [
#     "fs.ix.mass_transfer_coeff",
#     "fs.ix.freundlich_n",
#     "fs.ix.bv_50",
# ]
# self.theta_input = dict(
#     zip(self.theta_names, self.initial_guess_theta_dict.values())
# )
# self.theta_values = pd.DataFrame(
#     data=[[v for v in self.initial_guess_theta_dict.values()]],
#     columns=self.theta_names,
# )
# self.df_parmest = self.filtered_data[self.filtered_data.cb.isin(self.keep_cbs)]

# def parmest_regression(data):
#     """
#     Build function that is passed to parmest
#     """

#     cb = pyunits.convert(
#         data.cb.to_list()[0] * self.conc_units,
#         to_units=pyunits.kg / pyunits.m**3,
#     )()
#     cb_p = data.cb.to_list()[0] * self.conc_units
#     c0_p = self.c0 * self.conc_units
#     cx = value(cb_p / c0_p)
#     bv = data.bv.to_list()[0]

#     print(
#         f"\nPARMEST FOR:\n\t"
#         f"CURVE = {self.curve_id}\n\t"
#         f"REF = {self.ref}\n\t"
#         f"COMPOUND = {self.compound}\n\t"
#         f"BV = {bv}\n\t"
#         f"C0 = {c0_p}\n\t"
#         f"CB = {cb_p}\n\t"
#         f"CNORM = {cx}\n\t"
#     )

#     m_parmest = self.m.clone()
#     m_parmest.fs.ix.c_norm.fix(cx)
#     m_parmest.fs.ix.bv.set_value(bv)
#     self._calc_from_constr(m=m_parmest)
#     self.scale_it(m=m_parmest)
#     try:
#         m_parmest.fs.ix.initialize()
#     except:
#         pass
#     # Stores the instantiated model that is passed to parmest for debugging outside of the class
#     self.m_parmest = m_parmest.clone()
#     return m_parmest

# def SSE(m, data):
#     """
#     Objective function for parmest to minimize
#     We want the SSE of BV observed/predicted to be small
#     """
#     ix = self._get_ix_blk(m=m)
#     expr = (float(data.bv.iloc[0]) - ix.bv) ** 2
#     return expr * self.expr_sf  # expr_sf is the scaling factor for the SSE

# self.pest = parmest.Estimator(
#     parmest_regression,
#     self.df_parmest,
#     self.theta_names,
#     SSE,
#     tee=False,
#     diagnostic_mode=False,
#     solver_options={"max_iter": 10000},
#     **self.parmest_kwargs,
# )

# self.pest.objective_at_theta(
#     theta_values=self.theta_values, initialize_parmest_model=True
# )

# self.obj, self.theta = self.pest.theta_est()

# self.theta_dict = dict(zip(self.initial_guess_theta_dict.keys(), [*self.theta]))

# if hasattr(self, "results_dict"):
#     for k, v in self.theta_dict.items():
#         self.results_dict[f"{k}_theta"] = [
#             v for _ in range(len(self.results_dict["curve_id"]))
#         ]
#     self.results_dict["obj"] = [
#         self.obj for _ in range(len(self.results_dict["curve_id"]))
#     ]
