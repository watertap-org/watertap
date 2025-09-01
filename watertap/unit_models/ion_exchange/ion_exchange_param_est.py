import math
import pandas as pd
import numpy as np
from collections import defaultdict

import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    Var,
    Param,
    Constraint,
    ComponentUID,
    Suffix,
    assert_optimal_termination,
    check_optimal_termination,
    value,
    units as pyunits,
)

from idaes.core.util.model_diagnostics import *
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *

from watertap.core.solvers import get_solver


__author__ = "Kurban Sitterley"

solver = get_solver()


def update_rd(
    b,
    results_dict,
    components=[Var, Expression, Param],
):
    appended = list()

    for c in components:
        for v in b.component_objects(c):
            if v.is_indexed():
                for _, vi in v.items():
                    if vi.name in appended:
                        continue
                    results_dict[vi.name].append(value(vi))
                    appended.append(vi.name)
            else:
                if v.name in appended:
                    continue
                results_dict[v.name].append(value(v))
                appended.append(v.name)
    return results_dict


class IXParamEst:

    def __init__(
        self,
        data=None,
        build_function=None,
        obj_function=None,
        filter_data_function=None,
        initial_guess=dict(),
        thetas=list(),
        xlabel="c_norm",
        ylabel="bv",
        build_kwargs=dict(),
        filter_data_kwargs=dict(),
        parmest_kwargs=dict(),
    ):
        self.data = data
        self.build_function = build_function
        self.filter_data_function = filter_data_function
        self.obj_function = obj_function
        self.initial_guess = initial_guess

        self.thetas = thetas
        self.xlabel = xlabel
        self.ylabel = ylabel

        self.build_kwargs = build_kwargs
        self.filter_data_kwargs = filter_data_kwargs
        self.parmest_kwargs = parmest_kwargs

        self._vars_validated = False

    def build_model(self):
        self.model = self.build_function(initial_guess=self.initial_guess, **self.build_kwargs)
        if not self._vars_validated:
            self._validate_labels()
            self._vars_validated = True

    def filter_data(self):
        self.filtered_data = self.filter_data_function(
            self.data, **self.filter_data_kwargs
        )

    def create_experiment_list(self):
        if not hasattr(self, "filtered_data"):
            self.filter_data()

        self.experiment_list = []
        exp_dict = {
            "experiment_number": [
                x for x in range(len(self.filtered_data["filtered_x"].dropna()))
            ],
            "x": self.filtered_data["filtered_x"].dropna().values,
            "y": self.filtered_data["filtered_y"].dropna().values,
        }
        self.df_exp = pd.DataFrame.from_dict(exp_dict).set_index("experiment_number")

        for num in self.df_exp.index:
            self.experiment_list.append(
                IXExperiment(
                    self.df_exp,
                    num,
                    self.initial_guess,
                    self.build_function,
                    self.build_kwargs,
                    xlabel=self.xlabel,
                    ylabel=self.ylabel,
                    thetas=self.thetas,
                )
            )

    def run_parmest(self):

        if not hasattr(self, "experiment_list"):
            self.create_experiment_list()

        self.pestimator = parmest.Estimator(
            self.experiment_list, obj_function=self.obj_function, **self.parmest_kwargs
        )
        self.parmest_obj, self.parmest_theta = self.pestimator.theta_est()

        print("\nTHETA ESTIMATES:")
        for k, v in self.parmest_theta.items():
            n = k.split(".")[-1]
            if "e" in repr(v):
                print(f"  {n}: {v:.4e}")
            else:
                print(f"  {n}: {v:.4f}")

    def test_initial_guess(self):

        if not hasattr(self, "filtered_data"):
            self.filter_data()

        self.test_initial_guess_results = defaultdict(list)
        # self.build_kwargs["theta_dict"] = self.initial_guess
        # self.build_model()
        for x in self.filtered_data["filtered_x"].dropna().values:
            self.build_model()
            xvar = self.model.find_component(self.xlabel)
            xvar.fix(x)
            assert degrees_of_freedom(self.model) == 0
            try:
                self.model.fs.unit.initialize()
                results = solver.solve(self.model)
                assert_optimal_termination(results)
            except:
                print(f"Warning: could not solve for x = {x}")
                continue
            self.test_initial_guess_results = update_rd(
                self.model.fs.unit, self.test_initial_guess_results
            )

        self.test_initial_guess_results = pd.DataFrame.from_dict(
            self.test_initial_guess_results
        )
        self.test_initial_guess_results["flag"] = "test_initial_guess"

    def test_theta(self):

        if not hasattr(self, "parmest_theta"):
            self.run_parmest()

        self.test_theta_results = defaultdict(list)
        self.build_kwargs["theta_dict"] = {} 

        for x, y in zip(
            self.filtered_data["filtered_x"].dropna().values,
            self.filtered_data["filtered_y"].dropna().values,
        ):
            self.build_model()
            xvar = self.model.find_component(self.xlabel)
            xvar.fix(x)
            results = solver.solve(self.model)
            for k, v in self.parmest_theta.items():
                self.model.find_component(k).fix(v)
            # yvar = self.model.find_component(self.ylabel)
            # yvar.set_value(y)
            assert degrees_of_freedom(self.model) == 0
            try:
                # self.model.fs.unit.initialize()
                results = solver.solve(self.model)
                assert_optimal_termination(results)
            except:
                try:
                    # self.model.fs.unit.initialize()
                    results = solver.solve(self.model)
                    assert_optimal_termination(results)
                except:
                    print(f"\n\nWarning: could not solve for x = {x}\n")
                continue
            self.test_theta_results = update_rd(
                self.model.fs.unit, self.test_theta_results
            )
            # # save the x and y values that resulted in successful solve
            self.test_theta_results["filtered_x"].append(x)
            self.test_theta_results["filtered_y"].append(y)

        self.test_theta_results = pd.DataFrame.from_dict(self.test_theta_results)
        self.test_theta_results["flag"] = "test_theta"


    def test_theta_full(self, min_x=0.02, max_x=0.98, num_points=11):

        """
        Test theta across the full range of concentration ratios.
        """

        if not hasattr(self, "parmest_theta"):
            self.run_parmest()

        xs = np.linspace(min_x, max_x, num_points)

        self.test_theta_full_results = defaultdict(list)
        self.build_kwargs["theta_dict"] = self.parmest_theta

        for x in xs:
            self.build_model()
            xvar = self.model.find_component(self.xlabel)
            xvar.fix(x)
            results = solver.solve(self.model)
            for k, v in self.parmest_theta.items():
                self.model.find_component(k).fix(v)
            # yvar = self.model.find_component(self.ylabel)
            # yvar.set_value(y)
            assert degrees_of_freedom(self.model) == 0
            try:
                # self.model.fs.unit.initialize()
                results = solver.solve(self.model)
                assert_optimal_termination(results)
            except:
                try:
                    # self.model.fs.unit.initialize()
                    results = solver.solve(self.model)
                    assert_optimal_termination(results)
                except:
                    print(f"\n\nWarning: could not solve for x = {x}\n")
                continue
            self.test_theta_full_results = update_rd(
                self.model.fs.unit, self.test_theta_full_results
            )
            # save the x value that resulted in successful solve
            self.test_theta_full_results["filtered_x"].append(x)

        self.test_theta_full_results = pd.DataFrame.from_dict(self.test_theta_full_results)
        self.test_theta_full_results["flag"] = "test_theta_full"

    def make_results(self):

        if not hasattr(self, "test_initial_guess_results"):
            self.test_initial_guess()
        if not hasattr(self, "test_theta_results"):
            self.test_theta()
        if not hasattr(self, "test_theta_full_results"):
            self.test_theta_full()
        if not hasattr(self, "fit_statistics"):
            self.compute_fit_statistics()

        self.results = pd.concat(
            [self.test_initial_guess_results, self.test_theta_results, self.test_theta_full_results],
            ignore_index=True,
        )

        for col in self.filtered_data.columns:
            self.results[col] = self.filtered_data[col]

        for k, v in self.fit_statistics.items():
            self.results[k] = v

    def compute_fit_statistics(self):
        if not hasattr(self, "filtered_data"):
            self.filter_data()
        if not hasattr(self, "test_theta_results"):
            self.test_theta()

        self.fit_statistics = dict()
        y_exp = self.test_theta_results["filtered_y"].dropna().values
        y_pred = self.test_theta_results[self.ylabel].values
        n = len(y_exp)
        self.residuals = y_exp - y_pred
        self.SSE = np.sum(self.residuals**2)
        self.SS_tot = np.sum((y_exp - np.mean(y_exp)) ** 2)
        self.R_squared = 1 - (self.SSE / self.SS_tot)
        self.RMSE = np.sqrt(self.SSE / n)
        self.MAE = np.mean(np.abs(self.residuals))
        self.MSE = np.mean(self.residuals**2)
        self.MAPE = np.mean(np.abs(self.residuals / y_exp)) * 100

        self.fit_statistics["R^2"] = self.R_squared
        self.fit_statistics["MAE"] = self.MAE
        self.fit_statistics["MSE"] = self.MSE
        self.fit_statistics["MAPE (%)"] = self.MAPE
        self.fit_statistics["RMSE"] = self.RMSE
        self.fit_statistics["SS Residual"] = self.SSE
        self.fit_statistics["SS Total"] = self.SS_tot
        self.fit_statistics["n"] = n

    def _validate_labels(self):
        """
        Validate that the x and y labels exist on the unit model.
        If they are indexed, assume there is only one index and use the first index.
        This is a bit of a hack, but it works for now.
        """

        xvar = self.model.find_component(self.xlabel)
        if xvar is None:
            raise ValueError(
                f"Could not find x-variable {self.xlabel} on {self.model.name}."
            )
        if not isinstance(xvar, Var):
            raise ValueError(
                f"x-variable {self.xlabel} on {self.model.name} is not a Var."
            )
        # However it will find indexed variables with out the index specified
        if xvar.is_indexed():
            # assume there is only one index and it is the one we want
            if not len(xvar.index_set()) == 1:
                # I am not even sure if this is possible
                raise ValueError(
                    f"x-variable {self.xlabel} on {self.model.name} is indexed with more than one index."
                )
            idx = xvar.index_set().first()
            xvar = xvar[idx]
            self.xlabel = f"{self.xlabel}[{idx}]"

        # Repeat for ylabel
        yvar = self.model.find_component(self.ylabel)
        if yvar is None:
            raise ValueError(
                f"Could not find y-variable {self.ylabel} on {self.model.name}."
            )
        if not isinstance(yvar, Var):
            raise ValueError(
                f"y-variable {self.ylabel} on {self.model.name} is not a Var."
            )
        if yvar.is_indexed():
            if not len(yvar.index_set()) == 1:
                raise ValueError(
                    f"y-variable {self.ylabel} on {self.model.name} is indexed with more than one index."
                )
            idx = yvar.index_set().first()
            yvar = yvar[idx]
            self.ylabel = f"{self.ylabel}[{idx}]"


class IXExperiment(Experiment):
    def __init__(
        self,
        data,
        experiment_number,
        initial_guess,
        build_function,
        build_kwargs,
        xlabel="c_norm",
        ylabel="bv",
        thetas=["bv_50", "mass_transfer_coeff", "freundlich_n"],
        unit_blk=None,
    ):
        self.data = data
        self.experiment_number = experiment_number
        self.thetas = thetas
        self.initial_guess = initial_guess
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.unit_blk = unit_blk
        self.build_function = build_function
        self.build_kwargs = build_kwargs

        self.data_exp = self.data.loc[experiment_number]

    def create_model(self):
        self.model = self.build_function(
            initial_guess=self.initial_guess, **self.build_kwargs
        )
        # print(f"\ndof = {degrees_of_freedom(self.model)}")
        if self.unit_blk is None:
            self.unit_blk = self.model.fs.unit
        else:
            self.unit_blk = self.model.find_component(self.unit_blk)

        return self.model

    def finalize_model(self):

        model = self.model

        xval = self.data_exp["x"].astype(float)
        yval = self.data_exp["y"].astype(float)

        xvar = self.model.find_component(self.xlabel)
        if xvar is None:
            raise ValueError(
                f"Could not find x-variable {self.xlabel} on {self.unit_blk.name}."
            )
        if xvar.is_indexed():
            assert len(xvar.index_set()) == 1
            idx = xvar.index_set().first()
            xvar = xvar[idx]

        yvar = self.model.find_component(self.ylabel)
        if yvar is None:
            raise ValueError(
                f"Could not find y-variable {self.ylabel} on {self.unit_blk.name}."
            )
        if yvar.is_indexed():
            assert len(yvar.index_set()) == 1
            idx = yvar.index_set().first()
            yvar = yvar[idx]

        xvar.fix(xval)
        yvar.set_value(yval)

        print(f"\ndof = {degrees_of_freedom(model)}")
        try:
            self.unit_blk.initialize()
            results = solver.solve(model)
            # assert_optimal_termination(results)
        except:
            try:
                self.unit_blk.initialize()
                results = solver.solve(model)
                # assert_optimal_termination(results)
            except:
                pass
        # results = solver.solve(model)
        # assert_optimal_termination(results)

        for theta in self.thetas:
            v = self.unit_blk.find_component(theta)
            if v is None:
                raise ValueError(
                    f"Could not find theta {theta} on {self.unit_blk.name}"
                )
            v.unfix()

        print(f"\ndof = {degrees_of_freedom(model)}")

        return self.model

    def label_model(self):

        xval = self.data_exp["x"].astype(float)
        yval = self.data_exp["y"].astype(float)
        xvar = self.model.find_component(self.xlabel)
        yvar = self.model.find_component(self.ylabel)

        # if yvar.is_indexed():
        #     assert len(yvar.index_set()) == 1
        #     idx = yvar.index_set().first()
        #     yvar = yvar[idx]

        # if xvar.is_indexed():
        #     assert len(xvar.index_set()) == 1
        #     idx = xvar.index_set().first()
        #     xvar = xvar[idx]

        self.model.experiment_outputs = Suffix(direction=Suffix.LOCAL)
        self.model.experiment_outputs.update(
            [
                (
                    xvar,
                    xval,
                ),
                (yvar, yval),
            ]
        )
        self.model.unknown_parameters = Suffix(direction=Suffix.LOCAL)
        self.model.unknown_parameters.update(
            [
                (k, ComponentUID(k))
                for k in [getattr(self.unit_blk, theta) for theta in self.thetas]
            ]
        )
        return self.model

    def get_labeled_model(self):
        m = self.create_model()
        m = self.finalize_model()
        m = self.label_model()

        return m
