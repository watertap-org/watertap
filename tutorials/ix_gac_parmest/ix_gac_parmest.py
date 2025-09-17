import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment
from pyomo.environ import (
    Var,
    ComponentUID,
    ConcreteModel,
    Suffix,
    Reals,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
import idaes.core.util.scaling as iscale

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap.costing import WaterTAPCosting
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import (
    IonExchangeClark,
    GAC,
)

solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

marker_dict = {
    "ect2_sorbix_lc4": "o",
    "calgon_calres_2301": "x",
    "a694e": "d",
    "evoqua_psr2+": "^",
    "calgon_f600": "o",
    "calgon_reactivated": "s",
    "cabot_norit": "D",
    "calgon_virgin": "^",
    "evoqua_ultracarb": "v",
    "cabot_hydrodarco": "<",
    "jacobi": ">",
    "evoqua_aquacarb": "x",
}

color_dict = {"PFOA": "green", "PFBS": "blue"}

with open(f"{__location__}/data/pfas_properties.yaml", "r") as f:
    pfas_properties = yaml.load(f, Loader=yaml.FullLoader)

with open(f"{__location__}/data/resin_properties.yaml", "r") as f:
    resin_properties = yaml.load(f, Loader=yaml.FullLoader)

with open(f"{__location__}/data/gac_properties.yaml", "r") as f:
    gac_properties = yaml.load(f, Loader=yaml.FullLoader)


__all__ = [
    "BreakthroughExperiment",
    "AdsorptionParamEst",
    "build_ix_ocwd_pilot",
    "build_gac_ocwd_pilot",
    "filter_data",
    "plot_curve1",
    "plot_curve2",
    "pfas_properties",
    "gac_properties",
    "resin_properties",
]


class BreakthroughExperiment(Experiment):
    def __init__(
        self,
        data,
        experiment_number,
        initial_guess,
        build_function,
        build_kwargs,
        xlabel="",
        ylabel="",
        thetas=[],
    ):
        self.data = data
        self.experiment_number = experiment_number
        self.thetas = thetas
        self.initial_guess = initial_guess
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.build_function = build_function
        self.build_kwargs = build_kwargs

        self.data_exp = self.data.loc[experiment_number]

    def create_model(self):
        self.model = self.build_function(
            initial_guess=self.initial_guess, **self.build_kwargs
        )
        self.unit_blk = self.model.fs.unit

        return self.model

    def finalize_model(self):

        model = self.model

        xval = self.data_exp["x"].astype(float)
        yval = self.data_exp["y"].astype(float)

        xvar = self.model.find_component(self.xlabel)
        yvar = self.model.find_component(self.ylabel)

        xvar.fix(xval)
        yvar.set_value(yval)

        print(f"\nDOF After Initialization = {degrees_of_freedom(model)}")
        try:
            # sometimes it can take a few tries
            results = solver.solve(model)
            assert_optimal_termination(results)
        except:
            try:
                results = solver.solve(model)
                assert_optimal_termination(results)
            except:
                pass

        for theta in self.thetas:
            v = self.unit_blk.find_component(theta)
            if v is None:
                raise ValueError(
                    f"Could not find theta {theta} on {self.unit_blk.name}"
                )
            v.unfix()

        print(f"\nDOF Into Parmest = {degrees_of_freedom(model)}")

        return self.model

    def label_model(self):

        xval = self.data_exp["x"].astype(float)
        yval = self.data_exp["y"].astype(float)
        xvar = self.model.find_component(self.xlabel)
        yvar = self.model.find_component(self.ylabel)

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


class AdsorptionParamEst:

    def __init__(
        self,
        data=None,
        build_function=None,
        obj_function=None,
        filter_data_function=None,
        initial_guess=dict(),
        thetas=list(),
        xlabel="",
        ylabel="",
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

    def build_model(self):
        self.model = self.build_function(
            initial_guess=self.initial_guess, **self.build_kwargs
        )

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
                BreakthroughExperiment(
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

    def test_theta(self, xs=None, ys=None, theta_dict=None):

        # if not hasattr(self, "parmest_theta"):
        #     self.run_parmest()

        if xs is None:
            xs = self.filtered_data["filtered_x"].dropna().values
        if ys is None:
            ys = self.filtered_data["filtered_y"].dropna().values
        if theta_dict is None:
            theta_dict = self.parmest_theta

        self.test_theta_results = defaultdict(list)
        self.build_kwargs["theta_dict"] = theta_dict

        for x, y in zip(xs, ys):
            self.build_model()
            xvar = self.model.find_component(self.xlabel)
            xvar.fix(x)
            results = solver.solve(self.model)
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

            yvar = self.model.find_component(self.ylabel)
            self.test_theta_results[self.xlabel].append(x)
            self.test_theta_results["filtered_y"].append(y)
            self.test_theta_results[self.ylabel].append(value(yvar))

        self.test_theta_results = pd.DataFrame.from_dict(self.test_theta_results)

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

        self.fit_statistics["SSE"] = self.SSE
        self.fit_statistics["R^2"] = self.R_squared
        self.fit_statistics["MAE"] = self.MAE
        self.fit_statistics["MSE"] = self.MSE
        self.fit_statistics["MAPE (%)"] = self.MAPE
        self.fit_statistics["RMSE"] = self.RMSE
        self.fit_statistics["SS Total"] = self.SS_tot
        self.fit_statistics["n"] = n

        print(f"\n\n")
        header = f"{'Metric':<20} | {'Value':<20} |"
        print(header)
        print("=" * len(header))
        for metric, val in self.fit_statistics.items():
            valf = f"{val:.2f}"
            print(f"{metric:<20} | {valf:<20} |")


def build_ix_ocwd_pilot(
    species="PFOS", resin="calgon_calres_2301", theta_dict=dict(), **kwargs
):
    """
    Build for OCWD IX Pilot System
    """
    # PILOT SYSTEM INFO
    # FROM OCWD REPORT

    bed_depth = 29 * pyunits.inch
    flow_rate = 0.2 * pyunits.gallon / pyunits.minute
    ebct = 2.07 * pyunits.minute

    pfas_data = pfas_properties[species]
    resin_data = resin_properties[resin]

    ion_props = {
        "solute_list": [species],
        "mw_data": {
            "H2O": 0.018,
            species: pfas_data.get("mw", 0.350),
        },
        "molar_volume_data": {
            ("Liq", species): pfas_data.get("molar_volume", 0.0004),
        },
        "diffus_calculation": "HaydukLaudie",
        "charge": {species: -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": species,
        "add_steady_state_approximation": False,
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    # Set pilot feed conditions
    pfas_mol_flow = (
        pfas_data["c0"]
        * pyunits.ng
        / pyunits.L
        / (pfas_data["mw"] * pyunits.kg / pyunits.mol)
        * flow_rate
    )
    h2o_mol_flow = 55.5 * pyunits.mol / pyunits.L * flow_rate
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(h2o_mol_flow)
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", species].fix(
        pfas_mol_flow
    )
    ix.process_flow.properties_in[0].pressure.fix(101325)
    ix.process_flow.properties_in[0].temperature.fix(298)

    # Adjust bounds for pilot system
    ix.bed_depth.setlb(0)
    ix.bed_diameter.setlb(0)
    ix.ebct.setlb(0)
    ix.loading_rate.setlb(0)
    ix.loading_rate.setub(1)
    m.fs.unit.freundlich_n.setlb(1.05)
    m.fs.unit.freundlich_n.setub(100)
    m.fs.unit.bv.setlb(None)

    ix.resin_density.fix(resin_data["density"])
    ix.resin_diam.fix(resin_data["diameter"])
    ix.bed_depth.fix(bed_depth)
    ix.ebct.fix(ebct)
    ix.number_columns.fix(1)
    ix.c_norm.fix(0.5)

    for theta, val in theta_dict.items():
        ixv = m.find_component(theta)
        ixv.fix(val)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1 / value(h2o_mol_flow), index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1 / value(pfas_mol_flow), index=("Liq", species)
    )

    iscale.calculate_scaling_factors(m)

    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


def build_gac_ocwd_pilot(
    species="PFBS",
    media="calgon_f600",
    theta_dict=dict(),
    **kwargs,
):

    def apply_scaling(m):

        pfas_mol_flow = pyunits.convert(
            (pfas_data["c0"] * pyunits.ng / pyunits.L)
            / (pfas_data["mw"] * pyunits.kg / pyunits.mol)
            * flow_rate,
            to_units=pyunits.mol / pyunits.s,
        )
        rho = 1000 * pyunits.kg / pyunits.m**3
        h2o_mol_flow = pyunits.convert(
            (flow_rate * rho) / (0.018 * pyunits.kg / pyunits.mol),
            to_units=pyunits.mol / pyunits.s,
        )
        print(
            "H2O Mol Flow:",
            value(pyunits.convert(h2o_mol_flow, to_units=pyunits.mol / pyunits.s)),
        )
        print(
            "PFAS Mol Flow:",
            value(pyunits.convert(pfas_mol_flow, to_units=pyunits.mol / pyunits.s)),
        )
        iscale.set_scaling_factor(m.fs.unit.bed_area, 1 / value(bed_area))
        iscale.set_scaling_factor(m.fs.unit.bed_diameter, 1 / value(bed_diameter))
        iscale.set_scaling_factor(m.fs.unit.bed_volume, 1 / value(bed_volume))
        iscale.set_scaling_factor(m.fs.unit.bed_length, 1 / value(bed_depth))
        iscale.set_scaling_factor(m.fs.unit.ebct, 1 / value(ebct))
        iscale.set_scaling_factor(m.fs.unit.velocity_sup, 1 / value(loading_rate))
        iscale.set_scaling_factor(m.fs.unit.equil_conc, 1e5)
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].flow_vol_phase["Liq"],
            1 / value(flow_rate),
        )
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"],
            1 / value(h2o_mol_flow),
        )
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", species],
            1 / value(pfas_mol_flow),
        )
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].flow_mass_phase_comp["Liq", "H2O"],
            1 / value(flow_rate * rho),
        )
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].flow_mass_phase_comp[
                "Liq", species
            ],
            1 / value(pfas_mol_flow / pfas_data["mw"]),
        )
        iscale.set_scaling_factor(
            m.fs.unit.process_flow.properties_in[0].conc_mass_phase_comp[
                "Liq", species
            ],
            1 / value(pfas_data["c0"] * 1e-9),
        )
        iscale.constraint_scaling_transform(
            m.fs.unit.eq_equilibrium_concentration[0.0, species], 1e2
        )
        iscale.constraint_scaling_transform(m.fs.unit.eq_mass_adsorbed[species], 1e5)
        iscale.constraint_scaling_transform(m.fs.unit.eq_number_bi, 1e10)
        iscale.constraint_scaling_transform(m.fs.unit.eq_bed_area[0.0], 1e5)
        iscale.constraint_scaling_transform(m.fs.unit.eq_dg[0.0, species], 1e5)
        iscale.constraint_scaling_transform(m.fs.unit.eq_number_bi, 1e10)
        iscale.constraint_scaling_transform(
            m.fs.unit.eq_minimum_operational_time_cps, 1e-6
        )
        iscale.constraint_scaling_transform(m.fs.unit.eq_operational_time, 1e-6)
        iscale.constraint_scaling_transform(m.fs.unit.eq_bed_volumes_treated, 1e-6)
        iscale.constraint_scaling_transform(m.fs.unit.eq_residence_time, 0.1)
        iscale.constraint_scaling_transform(m.fs.unit.eq_min_residence_time_cps, 0.1)

    bed_diameter = pyunits.convert(3.04 * pyunits.inch, to_units=pyunits.m)
    bed_depth = pyunits.convert(52 * pyunits.inch, to_units=pyunits.m)
    bed_volume = pyunits.convert(377.93 * pyunits.inch**3, to_units=pyunits.m**3)
    bed_area = pyunits.convert(7.27 * pyunits.inch**2, to_units=pyunits.m**2)
    flow_rate = pyunits.convert(
        0.1625 * pyunits.gallon / pyunits.minute, to_units=pyunits.m**3 / pyunits.s
    )
    loading_rate = pyunits.convert(
        3.25 * pyunits.gallon / pyunits.minute / pyunits.ft**2,
        to_units=pyunits.m**3 / pyunits.s / pyunits.m**2,
    )
    ebct = pyunits.convert(10.07 * pyunits.minute, to_units=pyunits.s)

    pfas_data = pfas_properties[species]
    media_data = gac_properties[media]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=[species],
        mw_data={"H2O": 0.018, species: pfas_data.get("mw", 0.250)},
        ignore_neutral_charge=True,
        diffus_calculation="HaydukLaudie",
        molar_volume_data={
            ("Liq", species): pfas_data.get("molar_volume", 0.000250),
        },
    )
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="calculated",
        surface_diffusion_coefficient_type="fixed",
        cphsdm_calaculation_method="surrogate",
        add_trapezoidal_effluent_approximation=False,
    )
    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e12, index=("Liq", species)
    )
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.flow_vol_phase
    unit_feed.conc_mass_phase_comp
    iscale.calculate_scaling_factors(m)

    # feed specifications
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed = m.fs.unit.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 2e-8,
            ("conc_mass_phase_comp", ("Liq", species)): 4e-9,
        },
        hold_state=True,
    )

    m.fs.unit.freund_k.fix(3.7)
    m.fs.unit.freund_ninv.fix(0.8316)
    m.fs.unit.ds.fix(2e-13)
    m.fs.unit.particle_dens_app.fix(722)
    m.fs.unit.particle_dia.fix(0.00106)
    m.fs.unit.ebct.fix(300)  # seconds
    m.fs.unit.bed_voidage.fix(0.449)
    m.fs.unit.bed_length.fix(6)  # assumed
    m.fs.unit.shape_correction_factor.fix()

    m.fs.unit.conc_ratio_replace.fix(0.50)
    for v in m.fs.unit.component_objects(Var):
        v.domain = Reals
        v.setlb(None)
    print(f"dof = {degrees_of_freedom(m)}")

    assert degrees_of_freedom(m) == 0

    m.fs.unit.initialize()
    results = solver.solve(m)
    assert_optimal_termination(results)
    m.fs.properties.mw_comp[species] = pfas_data["mw"]
    m.fs.properties.molar_volume_phase_comp["Liq", species] = pfas_data["molar_volume"]

    # recalculate desired inlet flow and effluent target, based on conditions and new mw
    m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].unfix()
    m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", species].unfix()
    m.fs.unit.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_rate,
            ("conc_mass_phase_comp", ("Liq", species)): pfas_data["c0"] * 1e-9,
        },
        hold_state=True,
    )

    m.fs.unit.particle_dens_app.fix(media_data["density"])
    m.fs.unit.particle_dia.fix(media_data["diameter"])
    m.fs.unit.ebct.fix(ebct)
    m.fs.unit.bed_length.fix(bed_depth)
    m.fs.unit.velocity_sup.set_value(bed_depth / ebct)

    apply_scaling(m)

    for v, ig in theta_dict.items():
        vv = m.find_component(v)
        vv.fix(ig)

    # m.fs.unit.ds.setlb(1e-16)
    # m.fs.unit.ds.setub(1e-13)

    m.fs.unit.freund_ninv.setlb(0.01)
    m.fs.unit.freund_ninv.setub(0.99)

    m.fs.unit.freund_k.setlb(1e-8)
    m.fs.unit.freund_k.setub(100)

    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


def filter_data(
    data,
    c0=None,
    bv_min=100,
    cnorm_min_thresh=0.005,
    cnorm_max_thresh=1.2,
):
    """
    Method to filter extracted breakthrough curve
    Loops through provided BVs and Cbs and keeps (bv, cb) pairs that meet all conditions:
    1. Current cb is greater than the last kept cb (last_cb)
    2. Current cb is greater than the minimum threshold
    3. Current cb is less than the maximum threshold
    4. Associated BV is not zero
    5. (if the current cb is not the last) The next cb (i + 1) is some threshold greater than the current cb
    """

    if c0 is None:
        c0 = data.c0.iloc[0]

    last_cb = -1e6
    d = defaultdict(list)

    data.loc[data["c_norm"] < cnorm_min_thresh, "c_norm"] = cnorm_min_thresh
    data["cb"] = data["c_norm"] * data["c0"]

    d["all_bvs"] = data.bv.to_list()
    d["all_cbs"] = data.cb.to_list()
    d["all_cnorms"] = data.c_norm.to_list()

    for i, cb in enumerate(data.cb):

        # if cb == last_cb:
        #     d["excl_cbs"].append(cb)
        #     d["excl_bvs"].append(data.bv.iloc[i])
        #     continue

        if i != len(data.cb) - 1:
            if (
                cb >= last_cb
                and cb < c0
                and cb >= cnorm_min_thresh * c0
                and cb <= cnorm_max_thresh * c0
                and data.bv.iloc[i] > bv_min
            ):
                last_cb = cb
                d["keep_cbs"].append(cb)
                d["keep_bvs"].append(data.bv.iloc[i])
            else:
                d["excl_cbs"].append(cb)
                d["excl_bvs"].append(data.bv.iloc[i])
        else:
            if (
                cb >= last_cb
                and cb < c0
                # cb < c0
                and cb >= cnorm_min_thresh * c0
                and cb <= cnorm_max_thresh * c0
                and data.bv.iloc[i] >= bv_min
            ):
                last_cb = cb
                d["keep_cbs"].append(cb)
                d["keep_bvs"].append(data.bv.iloc[i])
            else:
                d["excl_cbs"].append(cb)
                d["excl_bvs"].append(data.bv.iloc[i])

    assert 0 not in d["keep_cbs"]

    d["keep_cnorms"] = [cb / c0 for cb in d["keep_cbs"] if not np.isnan(cb)]
    d["excl_cnorms"] = [cb / c0 for cb in d["excl_cbs"] if not np.isnan(cb)]

    d["filtered_y"] = d["keep_bvs"]
    d["filtered_x"] = [cb / c0 for cb in d["keep_cbs"]]

    filtered_data = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in d.items()]))
    filtered_data["c0"] = c0

    return filtered_data


def plot_curve1(
    data,
    data2=None,
    x="filtered_x",
    y="filtered_y",
    yleg=None,
    set_dict=dict(),
    fig=None,
    ax=None,
):
    if yleg is None:
        yleg = y
    if data2 is None:
        data2 = data.copy()

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots()
    ax.plot(
        data["all_bvs"],
        data["all_cnorms"],
        marker=None,
        # label="All extracted data",
        color="grey",
        alpha=0.25,
    )
    ax.scatter(
        data["keep_bvs"],
        data["keep_cnorms"],
        marker=".",
        label="Filtered Data",
        color="blue",
    )
    ax.scatter(
        data["excl_bvs"],
        data["excl_cnorms"],
        marker="x",
        label="Excluded Data",
        color="k",
    )
    ax.set_xlabel("Bed Volumes")
    ax.set_ylabel("C/C0")
    ax.set(**set_dict)
    ax.legend()
    ax.grid(zorder=0)
    return fig, ax


def plot_curve2(
    data,
    data2=None,
    x="filtered_x",
    y="filtered_y",
    yleg=None,
    set_dict=dict(),
    fig=None,
    ax=None,
):
    if yleg is None:
        yleg = y
    if data2 is None:
        data2 = data.copy()

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots()
    ax.plot(
        data["all_bvs"],
        data["all_cnorms"],
        marker=None,
        label="All extracted data",
        color="grey",
        alpha=0.25,
    )
    ax.scatter(
        data["keep_bvs"],
        data["keep_cnorms"],
        marker=".",
        label="Filtered data",
        color="blue",
    )
    ax.scatter(
        data["excl_bvs"],
        data["excl_cnorms"],
        marker="x",
        label="Excluded data",
        color="k",
    )
    ax.scatter(
        data2[x],
        data2[y],
        marker=".",
        label=yleg,
        color="green",
    )
    ax.set_xlabel("Bed Volumes")
    ax.set_ylabel("C/C0")
    ax.set(**set_dict)
    # ax.legend()
    ax.grid(zorder=0)
    return fig, ax
