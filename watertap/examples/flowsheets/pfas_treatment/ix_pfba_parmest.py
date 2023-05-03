from pyomo.environ import (
    ConcreteModel,
    units as pyunits,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import *  # model_statistics
from idaes.core.util.scaling import *
from idaes.core.util.model_diagnostics import DegeneracyHunter

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
    # IsothermType,
    # DiffusionControlType
)
import pyomo.contrib.parmest.parmest as parmest

from watertap.core.util.model_diagnostics.infeasible import *

import pandas as pd


"""
Clark model from Croll
Table 1, Table 2
PFBA

"""

f = "/Users/ksitterl/Documents/Python/watertap/watertap/kurby_watertap/ix/croll_pfba2.csv"

flow_in = 1.73 * pyunits.liter / pyunits.minute
flow_in = pyunits.convert(flow_in, to_units=pyunits.m**3 / pyunits.s)

c0 = 910.9 * pyunits.ng / pyunits.liter
c0 = pyunits.convert(c0, to_units=pyunits.kg / pyunits.m**3)

c_out = c0 * 0.5

bed_depth = 98.1 * pyunits.cm
bed_depth = pyunits.convert(bed_depth, to_units=pyunits.m)

col_diam = 7.73 * pyunits.cm
col_diam = pyunits.convert(col_diam, to_units=pyunits.m)

vel_bed = 36.9 * pyunits.cm / pyunits.min
vel_bed = pyunits.convert(vel_bed, to_units=pyunits.m / pyunits.s)

resin_diam = 675 * pyunits.um
resin_diam = pyunits.convert(resin_diam, to_units=pyunits.m)

resin_dens = 0.72 * pyunits.g / pyunits.milliliter
resin_dens = pyunits.convert(resin_dens, to_units=pyunits.kg / pyunits.liter)

bed_diam = 7.73 * pyunits.cm
bed_diam = pyunits.convert(bed_diam, to_units=pyunits.m)

guess_small_r = 1e-6
guess_A = 100
guess_freundlich_n = 2
guess_bv_50 = 17300

# Test data
cb = 534.3  # c_breakthru
bv = 20329  # bv at breakthru
tb = 3242830.34  # t_breakthru

target_ion = "PFBA_-"
mw = 214.04e-3


def build_it(cb=None, bv=None):
    global prop_in, prop_out, prop_regen, ix, pf
    target_ion = "PFBA_-"
    mw = 214.04e-3

    ion_props = {
        "solute_list": [target_ion],
        "diffusivity_data": {("Liq", target_ion): 0.49e-9},
        "mw_data": {"H2O": 0.018, target_ion: mw},
        "charge": {target_ion: -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "isotherm": "clark",
        "diffusion_control": "combination",
    }
    m.fs.ix = ix = IonExchange0D(**ix_config)
    pf = ix.process_flow
    prop_in = pf.properties_in[0]
    prop_out = pf.properties_out[0]
    prop_regen = ix.regeneration_stream[0]
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 10, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e9, index=("Liq", target_ion)
    )
    pf.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,
            ("conc_mass_phase_comp", ("Liq", target_ion)): c0,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    cb = cb * pyunits.ng / pyunits.liter
    cb = pyunits.convert(cb, to_units=pyunits.kg / pyunits.m**3)

    ix.resin_bulk_dens.fix(resin_dens)
    ix.bed_diam.fix(bed_diam)
    ix.c_breakthru[target_ion].fix(cb())
    ix.bed_depth.fix(bed_depth)
    ix.vel_bed.fix(vel_bed)
    ix.freundlich_exp[target_ion].set_value(2)
    ix.freundlich_exp[target_ion].setlb(2)
    ix.freundlich_exp[target_ion].setub(10)

    calculate_variable_from_constraint(ix.small_r, ix.eq_kinetic_param[target_ion])
    calculate_variable_from_constraint(ix.t_breakthru, ix.eq_bv)
    calculate_variable_from_constraint(ix.A, ix.eq_A[target_ion])
    ix.kinetic_param.setub(0.05)
    ix.kinetic_param.setlb(5e-9)
    ix.small_r.setub(1e-4)
    ix.small_r.setlb(1e-9)
    ix.bv_50.setub(100000)
    ix.bv_50.setlb(5000)
    ix.A.setlb(0.1)
    ix.A.setub(200)

    ix.bv.set_value(bv)
    # ix.bv.setub(bv*1.5)
    ix.bv_50.set_value(17394)

    m.fs.ix.small_r.fix(guess_small_r)
    m.fs.ix.A.fix(guess_A)
    m.fs.ix.freundlich_exp[target_ion].fix(guess_freundlich_n)

    print(f"DOF w all constraints active = {degrees_of_freedom(m)}")

    return m


def scale_it(m):
    calculate_scaling_factors(m)
    set_scaling_factor(ix.small_r, 1e7)
    set_scaling_factor(ix.t_breakthru, 1e-6)
    set_scaling_factor(ix.kinetic_param, 1e2)
    set_scaling_factor(ix.vel_bed, 1e2)
    set_scaling_factor(ix.freundlich_base, 1e2)
    set_scaling_factor(ix.bed_vol, 1e3)
    set_scaling_factor(ix.bv_50, 1e-5)
    set_scaling_factor(ix.bv, 1 / ix.bv())
    set_scaling_factor(ix.bed_diam, 1 / ix.bed_diam())
    set_scaling_factor(ix.c_breakthru[target_ion], 1 / ix.c_breakthru[target_ion]())
    set_scaling_factor(ix.A, 1e-1)
    set_scaling_factor(ix.mass_removed[target_ion], 1e3)
    constraint_scaling_transform(ix.eq_main2[target_ion], 1e6)
    constraint_scaling_transform(ix.eq_main[target_ion], 1e6)
    constraint_scaling_transform(
        ix.process_flow.properties_in[0.0].eq_mass_frac_phase_comp["Liq", "H2O"], 1e6
    )
    constraint_scaling_transform(ix.eq_mass_transfer_term[target_ion], 1e7)


def deactivate_it(m, dont_calc_effluent=True):
    ix = m.fs.ix
    ix.eq_A.deactivate()
    # ix.eq_main2.deactivate()
    # ix.eq_kinetic_param.deactivate()
    if dont_calc_effluent:
        # ix.eq_mass_transfer_target.deactivate()
        ix.eq_mass_transfer_term.deactivate()
        prop_out.deactivate()
        prop_regen.deactivate()
        pf.material_balances.deactivate()
        prop_out.flow_mol_phase_comp["Liq", "H2O"].fix(
            prop_in.flow_mol_phase_comp["Liq", "H2O"]()
        )
        prop_out.flow_mol_phase_comp["Liq", target_ion].fix(
            prop_in.flow_mol_phase_comp["Liq", target_ion]()
        )


def print_it(m):
    ix = m.fs.ix
    print(f"A = {ix.A()}")
    print(f"r = {ix.small_r()}")
    print(f"KF = {ix.freundlich_base[target_ion]()}")
    print(f"n = {ix.freundlich_exp[target_ion]()}")
    print(f"BV = {ix.bv()}")
    print(f"BV50 = {ix.bv_50()}")
    print(f"t_breakthru = {ix.t_breakthru()}")
    print(f"c_breakthru = {ix.c_breakthru[target_ion]()}")
    print(f"k_T = {ix.kinetic_param[target_ion]()}")


def test_single_build(cb=cb, bv=bv):
    m = build_it(cb=cb, bv=bv)

    deactivate_it(m)
    scale_it(m)

    print(f"DOF = {degrees_of_freedom(m)}")

    solver = get_solver()
    solver.options["max_iter"] = 10000
    solver.options["halt_on_ampl_error"] = "yes"
    results = solver.solve(
        m,
        symbolic_solver_labels=True,
        tee=False,
    )
    print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")
    print_it(m)
    print_infeasible_constraints(m)
    # model_debug(m, jacob=True)


def test_all_builds():

    df = pd.read_csv(f)
    bvs = df.bv.to_list()
    cbs = df.cb.to_list()
    tbs = df.t_breakthru.to_list()

    for bv, cb, tb in zip(bvs, cbs, tbs):

        print(f"\n------------------------\nBV= {bv}\nC_b = {cb}\nt_b = {tb}\n")
        m = build_it(cb=cb, bv=bv)
        deactivate_it(m)
        scale_it(m)

        # ix.initialize()
        print(f"DOF = {degrees_of_freedom(m)}")
        solver = get_solver()
        results = solver.solve(m)
        print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")
        print_it(m)


def ix_parm_est():
    theta_names = ["fs.ix.small_r", "fs.ix.A", f"fs.ix.freundlich_exp[{target_ion}]"]
    theta_values = pd.DataFrame(
        data=[[guess_small_r, guess_A, guess_freundlich_n]], columns=theta_names
    )
    # theta_names = ["fs.ix.small_r", "fs.ix.A", f"fs.ix.freundlich_exp[{target_ion}]", "fs.ix.bv_50"]
    # theta_values = pd.DataFrame(
    #     data=[[guess_small_r, guess_A, guess_freundlich_n, guess_bv_50]], columns=theta_names
    # )
    # theta_names = ["fs.ix.small_r", "fs.ix.A", "fs.ix.bv_50"]
    # theta_values = pd.DataFrame(
    #     data=[[guess_small_r, guess_A, guess_bv_50]], columns=theta_names
    # )

    def parmest_regression(data):
        cb = data.cb.to_list()[0]
        bv = data.bv.to_list()[0]
        print(f"\nRUNNING FOR:\n\tCB = {cb}\n\tBV = {bv}\n")
        m = build_it(cb=cb, bv=bv)

        deactivate_it(m)
        scale_it(m)
        # ix.initialize()

        print(f"DOF = {degrees_of_freedom(m)}")
        return m

    df = pd.read_csv(f)
    df = df[["bv", "cb", "t_breakthru"]].copy()

    expr_sf = 1e-6

    def SSE(m, df):
        expr = (float(df.bv) - m.fs.ix.bv) ** 2
        return expr * expr_sf

    pest = parmest.Estimator(
        parmest_regression,
        df,
        theta_names,
        SSE,
        tee=False,
        diagnostic_mode=False,
        # solver_options={"bound_push": 1e-8},
        solver_options={
            "max_iter": 6000,
            "bound_relax_factor": 0.0,
            "tol": 1e-08,
            "constr_viol_tol": 1e-08,
        },
    )

    pest.objective_at_theta(theta_values=theta_values, initialize_parmest_model=False)

    obj, theta = pest.theta_est()

    for k, v in theta.items():
        print(k, "=", v)


test_single_build(cb=cb, bv=bv)
test_all_builds()
ix_parm_est()
