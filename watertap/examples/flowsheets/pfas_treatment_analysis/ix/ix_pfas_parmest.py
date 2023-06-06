from pyomo.environ import (
    ConcreteModel,
    units as pyunits,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import *  # model_statistics
from idaes.core.util.scaling import *

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
)
import pyomo.contrib.parmest.parmest as parmest

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.costing import WaterTAPCosting

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from idaes.core.util.model_diagnostics import DegeneracyHunter

solver = get_solver()
from IPython.display import clear_output
from copy import deepcopy
from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)

import seaborn as sns

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

resin_dens_dict = {
    "a694": 0.72,
    "a600": 0.72,
    "ect_sorbix_a3f": 0.7,  # no data, assumed
    "ect2_sorbix_lc4": 0.7,  # no data, assumed
    "calgon_calres_2301": 0.67,
    "evoqua_psr2": 0.69,
}

resin_diam_dict = {
    "a694": 675e-6,
    "a600": 570e-6,
    "ect_sorbix_a3f": 675e-6,  # no data, assumed
    "ect2_sorbix_lc4": 675e-6,  # no data, assumed
    "calgon_calres_2301": 530e-6,  # "effective size" + std
    "evoqua_psr2": 700e-6,
}

import os

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
f = os.path.join(__location__, "ix_pfas_data.csv")

def main():
    global target_ion, resin_dens, resin_diam, mw, calc_from_constr, bv_50, set_bounds, flow_in, c0, vel_bed, bed_depth, ebct, conc_units, conc_units_str, fix_initial_guess, theta, ref, resin, compound, keep_bvs, keep_cbs, bvs, cbs, bv_pred, cb_pred, bv_infeas, cb_infeas
    df_pfas = pd.read_csv(f)
    refs = df_pfas.ref.unique()
    refs = ["croll_2023"]
    # refs = ["ocwd_2021"]
    # refs = ["franke_2021"]
    guess_kinetic_param = 1e-6
    guess_freundlich_n = 5

    guess_bv_50 = 20000
    fix_initial_guess = {
        "kinetic_param": guess_kinetic_param,
        f"freundlich_n": guess_freundlich_n,
    }
    set_bounds = {
        "kinetic_param": [5e-8, 1e-5],
        "freundlich_n": [1, 10],
    }
    calc_from_constr = {
        "mass_transfer_coeff": "eq_mass_transfer_coeff",
        "c_norm": "eq_c_breakthru",
    }
    # calc_from_constr = {}
    theta_dict = {}
    for ref in refs:
        theta_dict[ref] = {}
        print(ref)
        df_ref = df_pfas[df_pfas.ref == ref].copy()

        compounds = df_ref.compound.unique()
        for compound in compounds:
            theta_dict[ref][compound] = {}
            target_ion = compound + "_-"
            df_comp = df_ref[df_ref.compound == compound].copy()
            resins = df_comp.resin.unique()
            for resin in resins:
                theta_dict[ref][compound][resin] = {}
                df_resin = df_comp[
                    (df_comp.resin == resin) & (df_comp.c_norm != 0.5)
                ].copy()

                (
                    flow_in,
                    c0,
                    vel_bed,
                    bed_depth,
                    ebct,
                    conc_units,
                    conc_units_str,
                ) = get_ix_ref_conditions(df_resin)

                resin_dens = resin_dens_dict[resin]
                resin_diam = resin_diam_dict[resin]

                # print(ref, compound, resin, resin_dens, resin_diam, mw)

                mw = mw_dict[target_ion] * pyunits.kg / pyunits.mol
                bvs = df_resin.bv.to_list()
                cbs = df_resin.cb.to_list()
                c_norms = df_resin.c_norm.to_list()

                keep_bvs, keep_cbs, keep_cnorms = clean_bvs_cbs(bvs, cbs)
                if len(keep_bvs) <= 2:
                    keep_bvs, keep_cbs, keep_cnorms = clean_bvs_cbs(
                        bvs, cbs, round2=True
                    )
                set_bounds["bv_50"] = [0, max(keep_bvs)]
                guess_bv_50 = np.mean(keep_bvs)
                parmest_df = pd.DataFrame.from_dict(
                    {"bv": keep_bvs, "cb": keep_cbs, "c0": [c0() for _ in keep_cbs]}
                )
                cb_init = np.mean(keep_cbs)
                set_bounds["bv_50"] = [min(keep_bvs), max(keep_bvs)]
                obj, theta = ix_parmest(
                    target_ion=target_ion,
                    df=parmest_df,
                    cb_init=cb_init,
                    guess_kinetic_param=guess_kinetic_param,
                    guess_freundlich_n=guess_freundlich_n,
                    guess_bv_50=guess_bv_50,
                    set_bounds=set_bounds,
                    calc_from_constr=calc_from_constr,
                    conc_units=conc_units,
                    expr_sf=1e-9,
                )
                theta_dict[ref][compound][resin]["theta"] = {
                    "kinetic_param": theta[0],
                    "freundlich_n": theta[1],
                    "bv_50": theta[2],
                }
                theta_dict[ref][compound][resin]["obj"] = obj
                check_theta(theta, bv=np.mean(keep_bvs))


def build_it(
    target_ion=None,
    bv_50=None,
    flow_in=None,  # m3/s
    mw=None,
    charge=None,
    c0=None,
    cb=None,
    bv=None,
    tb=None,
    fix_initial_guess=None,
    set_bounds=None,
    calc_from_constr=None,
    dont_calc_effluent=True,
):

    ion_props = {
        "solute_list": [target_ion],
        "diffusivity_data": {("Liq", target_ion): 0.49e-9},
        "mw_data": {"H2O": 0.018, target_ion: mw},
        "charge": {target_ion: charge},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "isotherm": "freundlich",
    }
    m.fs.ix = ix = IonExchange0D(**ix_config)

    pf = ix.process_flow

    c0_flow_mol = pyunits.convert(c0 / mw * flow_in, to_units=pyunits.mol / pyunits.s)

    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1 / c0_flow_mol(), index=("Liq", target_ion)
    )

    c0 = pyunits.convert(c0, to_units=pyunits.kg / pyunits.m**3)
    pf.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,
            ("conc_mass_phase_comp", ("Liq", target_ion)): c0,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.underdrain_h.set_value(0)
    ix.distributor_h.set_value(0)
    ix.bed_expansion_frac_A.set_value(0)
    ix.bw_rate.set_value(0)

    # cb = cb * pyunits.ng / pyunits.liter
    cb = pyunits.convert(cb, to_units=pyunits.kg / pyunits.m**3)
    ix.c_breakthru[target_ion].fix(cb())

    ix.bv_50.fix(bv_50)
    ix.resin_bulk_dens.fix(resin_dens)
    ix.regen_dose.fix()
    ix.bed_porosity.fix()

    ix.bed_depth.fix(bed_depth)
    ix.vel_bed.fix(vel_bed)
    ix.resin_diam.fix(resin_diam)
    ix.number_columns.fix(1)
    ix.t_contact.setlb(0)

    if set_bounds is not None:
        if not isinstance(set_bounds, dict):
            raise Exception("set_bounds must be a dict")
        for k, v in set_bounds.items():
            print("set_bounds", k, v)
            ixv = getattr(ix, k)
            ixv.setlb(v[0])
            ixv.setub(v[1])

    if bv is None:
        ix.t_breakthru.set_value(tb)
        calculate_variable_from_constraint(ix.bv, ix.eq_bv)
    elif tb is None:
        ix.bv.setlb(0)
        ix.bv.setub(None)
        ix.bv.set_value(bv)
        calculate_variable_from_constraint(ix.t_breakthru, ix.eq_bv)

    if fix_initial_guess is not None:
        for k, v in fix_initial_guess.items():
            print("fix_initial_guess", k, v)
            ixv = getattr(ix, k)
            ixv.fix(v)

    if calc_from_constr is not None:
        for k, v in calc_from_constr.items():
            print("calc_from_constr", k, v)
            ixv = getattr(ix, k)
            ixc = getattr(ix, v)
            if all(c.is_indexed() for c in [ixv, ixc]):
                idx = [*ixv.index_set()]
                for i in idx:
                    calculate_variable_from_constraint(ixv[i], ixc[i])
            elif ixv.is_indexed():
                idx = [*ixv.index_set()]
                for i in idx:
                    calculate_variable_from_constraint(ixv[i], ixc)
            elif ixc.is_indexed():
                idx = [*ixc.index_set()]
                for i in idx:
                    calculate_variable_from_constraint(ixv, ixc[i])
            else:
                calculate_variable_from_constraint(ixc, ixv)

    if dont_calc_effluent:
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
        pf.mass_transfer_term[0, "Liq", target_ion].fix()

    print(f"DOF = {degrees_of_freedom(m)}")

    return m


def scale_it(m):
    ix = m.fs.ix
    target_ion = ix.config.target_ion
    prop_out = m.fs.ix.process_flow.properties_out[0]
    prop_regen = ix.regeneration_stream[0]
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", target_ion], 1e11)
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], 1)
    set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", target_ion], 1e15)

    calculate_scaling_factors(m)

    set_scaling_factor(ix.kinetic_param, 1e7)
    set_scaling_factor(ix.t_breakthru, 1e-6)
    set_scaling_factor(ix.mass_transfer_coeff, 1e4)
    set_scaling_factor(ix.vel_bed, 1e2)
    set_scaling_factor(ix.bed_vol_tot, 1e3)
    set_scaling_factor(ix.bv_50, 1e-5)
    set_scaling_factor(ix.bv, 1 / ix.bv())
    set_scaling_factor(ix.c_breakthru[target_ion], 1 / ix.c_breakthru[target_ion]())
    set_scaling_factor(ix.bed_capacity_param, 1e-1)
    set_scaling_factor(ix.c_traps, 1)
    constraint_scaling_transform(ix.eq_clark_2[target_ion], 1e6)
    constraint_scaling_transform(ix.eq_clark_1[target_ion], 1e6)
    constraint_scaling_transform(ix.eq_mass_transfer_regen[target_ion], 1e7)
    constraint_scaling_transform(ix.eq_mass_transfer_target_fr[target_ion], 1e7)


def get_ix_ref_conditions(df):
    conc_units_str = df.conc_units.unique()[0]
    numer = getattr(pyunits, conc_units_str.split("/")[0])
    denom = getattr(pyunits, conc_units_str.split("/")[1])
    conc_units = numer / denom
    c0 = df.c0.unique()[0]
    flow_in = df.flow_in.unique()[0]
    vel_bed = df.vel_bed.unique()[0]
    bed_depth = df.bed_depth.unique()[0]
    ebct = df.ebct.unique()[0]
    ebct = pyunits.convert(ebct * pyunits.min, to_units=pyunits.second)
    if np.isnan(vel_bed):
        vel_bed = bed_depth / ebct()
    return (
        flow_in * (pyunits.m**3 / pyunits.s),
        c0 * conc_units,
        vel_bed * (pyunits.m / pyunits.s),
        bed_depth * pyunits.m,
        ebct,
        conc_units,
        conc_units_str,
    )


def test_initial_guess(keep_bvs, keep_cbs, keep_cnorms):
    bv_pred = []
    cb_pred = []
    bv_infeas = []
    cb_infeas = []
    for bv, cb, c_norm in zip(keep_bvs, keep_cbs, keep_cnorms):
        print(bv, c_norm, cb)
        m = build_it(
            target_ion=target_ion,
            bv_50=bv_50,
            flow_in=flow_in,  # m3/s
            mw=mw,
            charge=-1,
            c0=c0,
            cb=cb * conc_units,
            bv=bv,
            tb=None,
            fix_initial_guess=fix_initial_guess,
            set_bounds=set_bounds,
            calc_from_constr=calc_from_constr,
            dont_calc_effluent=True,
        )
        ix = m.fs.ix
        print(f"DOF = {degrees_of_freedom(m)}")

        scale_it(m)
        ix.initialize()

        # solver = get_debug_solver(solver="ipopt")
        solver = get_solver()
        results = solver.solve(m, symbolic_solver_labels=True, tee=False)

        tc = results.solver.termination_condition
        if tc != "optimal":
            bv_infeas.append(bv)
            cb_infeas.append(cb)
            print(f"\n\n{tc.swapcase()}\n{bv, c_norm, cb}")
            print_infeasible_constraints(m)
            print("\n\n")
            # continue
        else:
            # break
            bv_pred.append(ix.bv())
            cb_pred.append(cb)
        print(f"\nSOLVE = {tc.swapcase()}\n")
    clear_output(wait=True)
    plot_bt_ig()
    plot_parity()


def clean_bvs_cbs(bvs, cbs, c0_thresh=0.01, plot_it=True, round2=False):
    last_cb = 0
    keep_cbs = []
    keep_bvs = []
    keep_cnorms = []
    excl_cbs = []
    excl_bvs = []
    # bvs = sorted(bvs)
    # cbs = sorted(cbs)
    for i, cb in enumerate(cbs):
        # print(i)
        if not round2:
            if i != len(cbs) - 1:
                if (
                    cb > last_cb
                    and cb < c0()
                    and cb >= c0_thresh * c0()
                    and bvs[i] > 0
                    # and cbs[i + 1] > cb * 0.9
                ):
                    last_cb = cb
                    keep_cbs.append(cb)
                    keep_bvs.append(bvs[i])
                else:
                    excl_cbs.append(cb)
                    excl_bvs.append(bvs[i])
            else:
                if cb > last_cb and cb < c0() and cb >= c0_thresh * c0() and bvs[i] > 0:
                    last_cb = cb
                    keep_cbs.append(cb)
                    keep_bvs.append(bvs[i])
                else:
                    excl_cbs.append(cb)
                    excl_bvs.append(bvs[i])
        else:
            if cb < c0() and bvs[i] > 0:
                last_cb = cb
                keep_cbs.append(cb)
                keep_bvs.append(bvs[i])
            else:
                excl_cbs.append(cb)
                excl_bvs.append(bvs[i])

    if plot_it:
        _, ax = plt.subplots()
        ax.scatter(excl_bvs, excl_cbs, marker="x", color="k", label="Excluded data")
        # ax.scatter(keep_bv, keep_cb, color="r", marker=".", label="Fitted data")
        ax.plot(
            keep_bvs, keep_cbs, color="r", marker=".", alpha=0.25, label="Fitted data"
        )
        ax.plot(bvs, [c0() for _ in bvs], linestyle="-.", alpha=0.25, label="Influent")
        title = (
            compound.swapcase()
            + " "
            + ref.replace("_", " ").title()
            + " "
            + resin.replace("_", " ").title()
        )
        ax.set_xlabel("BV")
        ylabe = compound.swapcase() + f" [{conc_units_str}]"
        ax.set_ylabel(ylabe)
        ax.set_title(title)
        ax.set_ylim([-0.01, c0() * 1.02])
        ax.legend()

    keep_cnorms = [cb / c0() for cb in keep_cbs]
    return keep_bvs, keep_cbs, keep_cnorms


def check_theta(theta, num_pts=25, test_cb_min=0.05, test_cb_max=0.99, bv=10000):
    global test_bvs, test_cbs
    kinetic_param = theta[0]
    n = theta[1]
    bv_50 = theta[2]

    test_cbs = [c0() * x for x in np.linspace(test_cb_min, test_cb_max, num_pts)]
    # test_bvs = np.linspace(min(bvs), max(bvs), num_pts).tolist()
    test_bvs = np.linspace(min(keep_bvs), max(keep_bvs), num_pts).tolist()

    # test_bvs.append(bv_50)
    # test_cbs.append(c0() * 0.5)

    # test_bvs = sorted(test_bvs)
    # test_cbs = sorted(test_cbs)

    bv_pred = []
    cb_pred = []
    bv_skip = []
    cb_skip = []

    m = build_it(
        target_ion=target_ion,
        flow_in=flow_in,
        charge=-1,
        c0=c0,
        mw=mw,
        bv_50=bv_50,
        cb=c0() * 0.5 * conc_units,
        bv=bv_50,
        fix_initial_guess={
            "kinetic_param": kinetic_param,
            "freundlich_n": n,
            "bv_50": bv_50,
        },
        calc_from_constr=calc_from_constr,
        set_bounds=set_bounds,
        dont_calc_effluent=True,
    )
    ix = m.fs.ix
    scale_it(m)
    ix.bv.setub(max(bvs))
    ix.initialize()
    prop_in = ix.process_flow.properties_in[0]

    # for bv, cb in zip(bvs, cbs):
    for bv, cb in zip(test_bvs, test_cbs):
        # bv = 10000
        print(f"\n------------------------\nBV= {bv}\nC_b = {cb}")
        ix.c_breakthru[target_ion].fix(cb * conc_units)
        calculate_variable_from_constraint(
            ix.c_norm[target_ion], ix.eq_c_breakthru[target_ion]
        )
        ix.c_norm[target_ion].fix()
        ix.c_breakthru[target_ion].unfix()
        ix.bv.set_value(bv)
        ix.initialize()
        print(f"DOF = {degrees_of_freedom(m)}")
        solver = get_solver()
        results = solver.solve(m)
        tc = results.solver.termination_condition
        print(f"MODEL SOLVE = {tc.swapcase()}")
        if tc != "optimal":
            print_infeasible_constraints(m)
            print_variables_close_to_bounds(m)
            bv_skip.append(bv)
            cb_skip.append(cb)
            continue

        bv_pred.append(ix.bv())
        cb_pred.append(cb)
        # ix.display()
    # clear_output(wait=True)
    fig, ax = plt.subplots()
    ax.scatter(bv_pred, cb_pred, marker="o", color="green", label="Predicted")
    ax.scatter(bv_skip, cb_skip, marker="x", color="black", label="Infeasible")
    ax.plot(keep_bvs, keep_cbs, marker=".", color="red", alpha=0.1, label="Extracted")
    ax.plot(
        [0, max(keep_bvs)],
        [c0() for _ in [0, max(keep_bvs)]],
        linestyle="-.",
        alpha=0.25,
        label="Influent",
    )
    ax.plot(test_bvs, test_cbs, linestyle="-.", alpha=0.25, label="Test Data")
    title = (
        compound.swapcase()
        + " "
        + ref.replace("_", " ").title()
        + " "
        + resin.replace("_", " ").title()
    )
    ax.set_xlabel("BV")
    ylabe = compound.swapcase() + f" [{conc_units_str}]"
    ax.set_ylabel(ylabe)
    ax.set_title(title)
    ax.set_ylim([-0.01, c0() * 1.02])
    ax.legend()
    ax.legend()
    plt.tight_layout()
    file = f"figs/theta_test_{ref}_{compound}_{resin}.png"
    # fig.savefig(file)


def ix_parmest(
    target_ion=None,
    df=None,
    cb_init=None,
    guess_kinetic_param=None,
    guess_freundlich_n=None,
    guess_bv_50=None,
    set_bounds={},
    calc_from_constr={},
    conc_units=None,
    expr_sf=1e-6,
):

    fix_initial_guess = {
        "kinetic_param": guess_kinetic_param,
        f"freundlich_n": guess_freundlich_n,
        "bv_50": guess_bv_50,
    }
    theta_names = ["fs.ix.kinetic_param", f"fs.ix.freundlich_n", "fs.ix.bv_50"]
    theta_values = pd.DataFrame(
        data=[[guess_kinetic_param, guess_freundlich_n, guess_bv_50]],
        columns=theta_names,
    )

    def parmest_regression(data):

        c0 = data.c0.to_list()[0] * conc_units
        cb = data.cb.to_list()[0] * conc_units
        bv = data.bv.to_list()[0]

        print(
            f"\nRUNNING FOR:\n\tREF = {ref}\n\tCOMPOUND = {target_ion}\n\tC0 = {c0}\n\tCB = {cb}\n\tBV = {bv}\n"
        )

        m = build_it(
            target_ion=target_ion,
            flow_in=flow_in,
            charge=-1,
            c0=c0,
            mw=mw,
            bv_50=guess_bv_50,
            cb=cb_init * conc_units,
            bv=bv,
            fix_initial_guess=fix_initial_guess,
            calc_from_constr=calc_from_constr,
            set_bounds=set_bounds,
        )

        ix = m.fs.ix
        ix.resin_bulk_dens.fix(resin_dens)
        scale_it(m)

        ix.initialize()
        # ix.display()

        print(f"DOF = {degrees_of_freedom(m)}")

        return m

    def SSE(m, data):
        expr = (float(data.bv) - m.fs.ix.bv) ** 2
        return expr * expr_sf

    pest = parmest.Estimator(
        parmest_regression,
        df,
        theta_names,
        SSE,
        tee=False,
        diagnostic_mode=False,
    )

    # pest.obj_function
    # clear_output(wait=True)
    pest.objective_at_theta(theta_values=theta_values, initialize_parmest_model=True)
    # clear_output(wait=True)
    obj, theta = pest.theta_est()
    # clear_output(wait=True)
    for k, v in theta.items():
        print(k, "=", v)
    return obj, theta


def plot_bt_ig(save_it=False):
    t = target_ion.split("_")[0].upper()
    fig, ax = plt.subplots()

    ax.plot(
        # [b for b in bvs if b not in bv_infeas],
        # [c for c in cbs if c not in cb_infeas],
        bvs,
        cbs,
        color="b",
        ls=":",
        marker=".",
        alpha=0.2,
        label="Actual",
    )
    ax.scatter(bv_pred, cb_pred, color="r", label="Predicted")
    ax.scatter(bv_infeas, cb_infeas, color="k", marker="x", label="Infeasible")

    # ax.scatter(bv_keep, bv_pred, color="r", label="Predicted")
    # ax.scatter(bv_skip, bv_skip, color="k", marker="x", label="Infeasible")
    ax.set_xlabel("BV")
    ax.set_ylabel(f"{t} Conc. [{conc_units_str}]")
    ax.legend()
    # ax.set_ylabel("Predicted BV")
    # ax.set_xlabel("Actual BV")
    # ax.set_title(f"{t}\nPredicted vs. Actual BVs with Initial Guess")
    title = f'{ref.replace("_", " ").title()} - {t}\nInitial Guess'
    ax.set_title(title)
    plt.tight_layout()
    if save_it:
        fig.savefig(f"figs/{ref}_{t.lower()}_bv_initial_guess.png")


def plot_parity():
    t = target_ion.split("_")[0].upper()
    fig, ax = plt.subplots()
    ax.plot(
        bvs,
        bvs,
        ":",
    )
    ax.scatter(
        [bvs[i] for i, _ in enumerate(bv_pred)], bv_pred, color="r", label="Predicted"
    )
    ax.scatter(
        [bvs[i] for i, _ in enumerate(bv_infeas)],
        bv_infeas,
        color="k",
        marker="x",
        label="Infeasible",
    )
    ax.legend()
    title = f'{ref.replace("_", " ").title()} - {t}\nInitial Guess Parity'
    ax.set_title(title)


if __name__ == "__main__":
    main()
