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
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

import numpy as np
import pandas as pd
import pyomo.contrib.parmest.parmest as parmest
from os.path import join, abspath, dirname

__author__ = "Hunter Barber"

# obtain experimental RSSCT data from csv
source_name = "Serrano Water District"
data = pd.read_csv("F400_PFOA.csv")

# initial guess for regressed variables
guess_freund_k = 10
guess_freund_ninv = 0.7

# set up solver
solver = get_solver()
solver_options = {
    "max_iter": 50000,
}
solver.options = solver_options


def main():

    print(
        "########################################## model initial build ##########################################"
    )
    # testing model at parameter values used to initialize
    m = model_init_build()
    m.fs.gac.display()
    assert False
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
    def SSE(model, data):
        expr = (float(data[f"{source_name}_X"]) - model.fs.gac.bed_volumes_treated) ** 2
        return expr

    print(
        "########################################## model regression ##########################################"
    )
    # Create an instance of the parmest estimator
    pest = parmest.Estimator(
        model_parmest_build,
        data,
        theta_names,
        SSE,
        tee=True,
        solver_options=solver_options,
        diagnostic_mode=True,
    )

    pest.objective_at_theta(
        theta_values=theta_values,
        initialize_parmest_model=True,
    )

    # Parameter estimation
    solver_est = SolverFactory("ef_ipopt")
    solver_est.options = solver_options
    obj, theta = pest.theta_est(solver=solver_est)

    print("The SSE at the optimal solution is %0.6f" % obj)
    print("The values for the parameters are as follows:")
    for k, v in theta.items():
        print(k, "=", v)

    print(
        "########################################## model rebuild ##########################################"
    )

    # rebuild model across CP to view regression results
    # theta = [guess_freund_k, guess_freund_ninv, guess_ds]
    model_regressed_build(theta)


def model_init_build():

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
        surface_diffusion_coefficient_type="calculated",
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
    m.fs.gac.particle_porosity.fix(0.69)
    m.fs.gac.tort.fix()
    m.fs.gac.spdfr.fix()
    m.fs.gac.shape_correction_factor.fix()
    m.fs.gac.conc_ratio_start_breakthrough = 0.001

    # scaling
    m = _new_scaling_factors(m)
    iscale.calculate_scaling_factors(m)

    # initialization
    optarg = solver.options
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    m.fs.gac.initialize(optarg=optarg)

    _model_solve(m)

    return m


def model_parmest_build(data):

    print(f'running regression case {int(data["data_iter"])}')
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
    m.fs.gac.particle_dens_app.fix(682)
    m.fs.gac.particle_dia.fix(0.000127)
    # adsorber bed specifications
    m.fs.gac.ebct.fix(13.68)
    m.fs.gac.bed_voidage.fix(0.2037)
    m.fs.gac.velocity_sup.fix(0.002209)
    # design spec
    m.fs.gac.conc_ratio_replace.fix(0.20)
    # parameters
    m.fs.gac.shape_correction_factor.fix(1)
    m.fs.gac.a0.fix(6.31579e0)
    m.fs.gac.a1.fix(56.8421)
    m.fs.gac.b0.fix(0.865453)
    m.fs.gac.b1.fix(0.157618)
    m.fs.gac.b2.fix(0.444973)
    m.fs.gac.b3.fix(0.001650)
    m.fs.gac.b4.fix(0.148084)
    m.fs.gac.conc_ratio_start_breakthrough = 0.001

    # scaling
    m = _new_scaling_factors(m)
    iscale.calculate_scaling_factors(m)

    # initialization
    optarg = solver.options
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    m.fs.gac.initialize(optarg=optarg)

    # re-fix to conc_ratio spec from data
    m.fs.gac.conc_ratio_replace.fix(float(data[f"{source_name}_Y"]))

    return m


def model_regressed_build(theta):

    conc_ratio_input = np.linspace(0.01, 0.95, 5)
    conc_ratio_list = []
    bed_volumes_treated_list = []

    df_main = pd.DataFrame(
        {
            "conc_ratio": [],
            "bed_volumes_treated": [],
        }
    )

    for conc_ratio in conc_ratio_input:
        print("running conc_ratio", conc_ratio)
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
        m.fs.gac.freund_k.fix(theta[0])
        m.fs.gac.freund_ninv.fix(theta[1])
        m.fs.gac.ds.fix(theta[2])
        # gac particle specifications
        m.fs.gac.particle_dens_app.fix(682)
        m.fs.gac.particle_dia.fix(0.000127)
        # adsorber bed specifications
        m.fs.gac.ebct.fix(13.68)
        m.fs.gac.bed_voidage.fix(0.2037)
        m.fs.gac.velocity_sup.fix(0.002209)
        # design spec
        m.fs.gac.conc_ratio_replace.fix(conc_ratio)
        # parameters
        m.fs.gac.shape_correction_factor.fix(1)
        m.fs.gac.a0.fix(6.31579e0)
        m.fs.gac.a1.fix(56.8421)
        m.fs.gac.b0.fix(0.865453)
        m.fs.gac.b1.fix(0.157618)
        m.fs.gac.b2.fix(0.444973)
        m.fs.gac.b3.fix(0.001650)
        m.fs.gac.b4.fix(0.148084)
        m.fs.gac.conc_ratio_start_breakthrough = 0.001

        # scaling
        m = _new_scaling_factors(m)
        iscale.calculate_scaling_factors(m)

        # initialization
        optarg = solver.options
        m.fs.feed.initialize(optarg=optarg)
        propagate_state(m.fs.s01)
        m.fs.gac.initialize(optarg=optarg)

        m = _model_solve(m)

        conc_ratio_list.append(m.fs.gac.conc_ratio_replace.value)
        bed_volumes_treated_list.append(m.fs.gac.bed_volumes_treated.value)

    df = pd.DataFrame(
        {
            "conc_ratio": conc_ratio_list,
            "bed_volumes_treated": bed_volumes_treated_list,
        }
    )
    df.to_csv("pfas_regression_results.csv")

    import matplotlib.pyplot as plt

    plt.plot(bed_volumes_treated_list, conc_ratio_list, label="regression results")
    plt.plot(
        data[f"{source_name}_X"], data[f"{source_name}_Y"], label=f"{source_name} data"
    )
    plt.legend()
    plt.show()


def _new_scaling_factors(m):

    set_scaling_factor(m.fs.gac.mass_adsorbed, 1e10)
    set_scaling_factor(m.fs.gac.mass_adsorbed_saturated, 1e9)
    set_scaling_factor(m.fs.gac.N_Bi, 1e-1)
    set_scaling_factor(m.fs.gac.kf, 1e4)
    set_scaling_factor(m.fs.gac.min_N_St, 1e-2)
    set_scaling_factor(m.fs.gac.equil_conc, 1e6)
    set_scaling_factor(m.fs.gac.gac_usage_rate, 1e9)
    set_scaling_factor(m.fs.gac.dg, 1e-5)
    set_scaling_factor(m.fs.gac.bed_length, 1e2)
    set_scaling_factor(m.fs.gac.bed_diameter, 1e3)
    set_scaling_factor(m.fs.gac.bed_area, 1e5)
    set_scaling_factor(m.fs.gac.bed_volume, 1e6)
    set_scaling_factor(m.fs.gac.bed_mass_gac, 1e4)
    set_scaling_factor(m.fs.gac.residence_time, 1)
    set_scaling_factor(m.fs.gac.velocity_int, 1e2)
    set_scaling_factor(m.fs.gac.velocity_sup, 1e2)
    set_scaling_factor(m.fs.gac.conc_ratio_avg, 1e1)

    return m


def _model_solve(model):

    # check model
    assert_units_consistent(model)  # check that units are consistent
    assert degrees_of_freedom(model) == 0

    # solve simulation
    results = solver.solve(model, tee=False)
    term_cond = results.solver.termination_condition
    print("termination condition:", term_cond)

    # if not optimal try resolve
    if not pyo.check_optimal_termination(results):
        iter = 1
        print("Trouble solving model, trying resolve")
        while not pyo.check_optimal_termination(results):
            results = solver.solve(model, tee=False)
            term_cond = results.solver.termination_condition
            print("termination condition:", term_cond)
            if iter == 5:
                print(f"Could not solve model in {iter} iterations")
                break
            iter = iter + 1

    return model


def _model_debug(model):

    check_jac(model)

    model.obj = pyo.Objective(expr=0)

    # initial point
    solver.options["max_iter"] = 0
    solver.solve(model, tee=False)
    dh = DegeneracyHunter(model, solver=pyo.SolverFactory("cbc"))
    dh.check_residuals(tol=1e-8)
    dh.check_variable_bounds(tol=1e-8)

    # solved model
    solver.options["max_iter"] = 10000
    solver.solve(model, tee=False)
    badly_scaled_var_list = badly_scaled_var_generator(model, large=1e1, small=1e-1)
    for x in badly_scaled_var_list:
        print(f"{x[0].name}\t{x[0].value}\tsf: {get_scaling_factor(x[0])}")
    dh.check_residuals(tol=1e-8)
    dh.check_variable_bounds(tol=1e-8)
    dh.check_rank_equality_constraints(dense=True)
    ds = dh.find_candidate_equations(verbose=True, tee=True)
    ids = dh.find_irreducible_degenerate_sets(verbose=True)

    """
    variables_near_bounds_list = variables_near_bounds_generator(model)
    for x in variables_near_bounds_list:
        print(x, x.value)
    """

    return model


def check_jac(model):
    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(model, min_scale=1e-8)
    # cond_number = iscale.jacobian_cond(model, jac=jac_scaled)  # / 1e10
    # print("--------------------------")
    print("Extreme Jacobian entries:")
    extreme_entries = iscale.extreme_jacobian_entries(
        model, jac=jac_scaled, zero=1e-20, large=10
    )
    for val, var, con in extreme_entries:
        print(val, var.name, con.name)
    print("--------------------------")
    print("Extreme Jacobian columns:")
    extreme_cols = iscale.extreme_jacobian_columns(model, jac=jac_scaled)
    for val, var in extreme_cols:
        print(val, var.name)
    print("------------------------")
    print("Extreme Jacobian rows:")
    extreme_rows = iscale.extreme_jacobian_rows(model, jac=jac_scaled)
    for val, con in extreme_rows:
        print(val, con.name)


if __name__ == "__main__":
    main()
