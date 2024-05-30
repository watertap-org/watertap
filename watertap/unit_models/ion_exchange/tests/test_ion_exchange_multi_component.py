import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools
from pyomo.environ import (
    ConcreteModel,
    Expression,
    value,
    Var,
    Param,
    Constraint,
    Block,
    exp,
    SolverFactory,
    units as pyunits,
    assert_optimal_termination,
    check_optimal_termination,
)
import pyomo.contrib.parmest.parmest as parmest


from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models import IonExchangeMultiComp

# from watertap.unit_models
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.costing import WaterTAPCosting
from copy import deepcopy

import matplotlib.pyplot as plt
from idaes.core.util.constants import Constants


from watertap.core.solvers import get_solver

solver = get_solver()


def build_ix(
    flow_in=0.04381,
    fix_vars=dict(),
    set_vars=dict(),
    target_component="",
    solute_list=list(),
    reactive_ions=list(),
    c0_dict={},
    diff_data_dict={},
    charge_dict={},
    mass_transfer_coeff_dict={},
    freundlich_n_dict={},
    bv_50_dict={},
    number_traps=5,
    regenerant="NaCl",
    deactivate_constr=[],
    dont_calc_effluent=True,
):

    global c0_mol_flow_sf_dict
    assert target_component in solute_list

    mw_dict = {"comp1": 0.21404, "comp2": 0.164031, "inert": 0.04}

    # assert all(d.keys() == solute_list for d in [c0_dict, diff_data_dict, charge_dict])

    flow_in = flow_in * pyunits.m**3 / pyunits.s
    rho = 1000 * pyunits.kg / pyunits.m**3
    mw_water = 0.018 * pyunits.kg / pyunits.mol

    c0_mol_flow_dict = dict()
    c0_mol_flow_sf_dict = dict()
    mw_data = {"H2O": mw_water}

    for s in solute_list:
        mw_data[s] = mw_dict[s] * pyunits.kg / pyunits.mol
        c0_mol_flow_dict[s] = pyunits.convert(
            c0_dict[s] * (pyunits.kg / pyunits.m**3) / mw_data[s] * flow_in,
            to_units=pyunits.mol / pyunits.s,
        )
        c0_mol_flow_sf_dict[s] = 1 / c0_mol_flow_dict[s]()

    water_mol_flow = pyunits.convert(
        (flow_in * rho) / mw_water, to_units=pyunits.mol / pyunits.s
    )
    water_mol_flow_sf = 1 / water_mol_flow()

    ion_props = {
        "solute_list": solute_list,
        "diffusivity_data": diff_data_dict,
        "mw_data": mw_data,
        "charge": charge_dict,
    }

    ix_config = {
        "target_component": target_component,
        "reactive_ions": reactive_ions,
        "regenerant": regenerant,
        "number_traps": number_traps,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # print(*ion_props.items(), sep="\n")
    m.fs.properties = MCASParameterBlock(**ion_props)

    ix_config["property_package"] = m.fs.properties

    m.fs.ix = ix = IonExchangeMultiComp(**ix_config)

    pf = ix.process_flow
    pf.properties_in[0].pressure.fix(101325)
    pf.properties_in[0].temperature.fix(298)

    for s in solute_list:
        pf.properties_in[0].flow_mol_phase_comp["Liq", s].fix(c0_mol_flow_dict[s]())
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", c0_mol_flow_sf_dict[s], index=("Liq", s)
        )
    pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(water_mol_flow)
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", water_mol_flow_sf, index=("Liq", "H2O")
    )

    for k, v in fix_vars.items():
        ixv = getattr(ix, k)
        if ixv.is_indexed():
            ixv[target_component].fix(v)
        else:
            ixv.fix(v)

    for k, v in mass_transfer_coeff_dict.items():
        ix.mass_transfer_coeff[k].fix(v)

    for k, v in freundlich_n_dict.items():
        ix.freundlich_n[k].fix(v)

    for k, v in bv_50_dict.items():
        ix.bv_50[k].fix(v)

    for c in deactivate_constr:
        print(f"deactivate {c}")
        ixc = getattr(ix, c)
        ixc.deactivate()

    for k, v in set_vars.items():
        ixv = getattr(ix, k)
        ixv.set_value(v)

    ix.bv.set_value(
        value(ix.bv_50[target_component]) / value(ix.c_norm[target_component])
    )
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
        if ix.config.regenerant != "single_use":
            # ix.eq_mass_transfer_target_fr.deactivate()
            ix.eq_mass_transfer_regen.deactivate()
        pf.mass_transfer_term[0, "Liq", target_component].fix()
        ix.eq_mass_transfer_reactive_ions.deactivate()

    print(f"DOF: {degrees_of_freedom(m)}")

    return m


def plot_traps(m, x_plot=None, run_curve=None, cmap="viridis"):
    global tb_traps_dict, bv_traps_dict, c_norm_dict

    for c in m.fs.component_objects():
        if isinstance(c, IonExchangeMultiComp):
            ix = c
            break
    comps = ix.config.property_package.solute_set
    target_component = ix.config.target_component
    cmap = plt.cm.get_cmap(cmap)
    rgbas = itertools.cycle([cmap(x) for x in np.linspace(0, 0.9, len(comps))])
    tb_traps_dict = dict()
    bv_traps_dict = dict()
    c_norm_dict = dict()
    fig, ax = plt.subplots()

    if x_plot == "bv":
        ax.scatter(
            ix.bv_50[target_component](),
            0.5,
            marker="X",
            color="red",
            label="BV$_{50}$, Target",
            zorder=100,
        )
        ax.scatter(
            ix.bv(),
            ix.c_norm[target_component](),
            marker="*",
            color="blueviolet",
            label="C$_{b}$, Target",
            zorder=150,
        )
    else:
        ax.scatter(
            ((ix.bv_50[target_component]() * ix.bed_depth()) / ix.loading_rate())
            / 3600
            / 24,
            0.5,
            marker="X",
            color="red",
            label="BV$_{50}$, Target",
            zorder=100,
        )
        ax.scatter(
            ((ix.bv() * ix.bed_depth()) / ix.loading_rate()) / 3600 / 24,
            ix.c_norm[target_component](),
            marker="*",
            color="blueviolet",
            label="C$_{b}$, Target",
            zorder=100,
        )
    for c in ix.reactive_ion_set:
        print(c)
        tb_traps_dict[c] = list()
        bv_traps_dict[c] = list()
        c_norm_dict[c] = list()
        color = next(rgbas)

        tb_traps = []
        bv_traps = []
        c_traps = []
        if x_plot == "bv":
            # print(c, "bv_traps")
            for (i, j), v in ix.bv_traps.items():
                if i != c:
                    continue
                # print(i, j, v)
                bv_traps.append(ix.bv_traps[i, j]())
                ax.vlines(
                    ix.bv_traps[i, j](), 0, ix.c_traps[i, j](), color=color, alpha=0.25
                )
        else:
            # print(c, "tb_traps")
            for (i, j), v in ix.tb_traps.items():
                if i != c:
                    continue
                tb_traps.append(
                    pyunits.convert(ix.tb_traps[i, j], to_units=pyunits.day)()
                )
                ax.vlines(tb_traps[-1], 0, ix.c_traps[i, j](), color=color, alpha=0.25)
        for (i, j), v in ix.c_traps.items():
            # print(i, j, v())
            if i != c:
                continue
            c_traps.append(v())
        if x_plot == "bv":
            x = bv_traps
            xlabe = "BV"
        else:
            x = tb_traps
            xlabe = "Breakthrough time [day]"
        tb_traps_dict[c] = tb_traps
        bv_traps_dict[c] = bv_traps
        c_norm_dict[c] = c_traps

        ax.plot(
            x,
            c_traps,
            marker=".",
            color=color,
            label=f"{c} - Breakthrough",
        )

        ax.set_xlabel(xlabe)
        ax.set_ylabel("C/C$_0$")
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.legend()
        ax.set_ylim([0, 1.02])
        # plt.tight_layout()
    #     return results_dict


def scale_model(
    m, autoscale_fixed=True, scale_from_value=list(), custom_scaling=dict()
):
    """
    Scale the model

    """

    def _get_comp_list(blk, comp=Var, skip_list=[]):
        cs = []
        split_name = blk.name + "."
        skip_list += ["ref", "process_flow", "regeneration"]
        for c in blk.component_objects(comp):
            if any(s in c.name for s in skip_list):
                continue
            cs.append(c.name.split(split_name)[1])
        return cs

    ix = m.fs.ix

    calculate_scaling_factors(m)

    if autoscale_fixed:
        print("autoscale_fixed")
        for v in _get_comp_list(ix):
            # print(v)
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

    if isinstance(scale_from_value, list):
        for v in scale_from_value:
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

    for k, v in custom_scaling.items():
        print(f"{k, v} custom_scaling")
        ixv = getattr(ix, k)
        set_scaling_factor(ixv, v)

    for i, c in ix.eq_clark.items():

        constraint_scaling_transform(c, 1e-2)


## Dummy data for model
def main():
    solute_list = ["comp1", "comp2", "inert"]
    reactive_ions = ["comp1", "comp2"]

    c0_dict = {"comp1": 1e-6, "comp2": 8e-6, "inert": 1e-3}
    diff_data_dict = {
        ("Liq", "comp1"): 10e-10,
        ("Liq", "comp2"): 5e-9,
        ("Liq", "inert"): 3e-9,
    }
    charge_dict = {"comp1": -1, "comp2": -1}

    mass_transfer_coeff_dict = {"comp1": 0.75, "comp2": 0.6}
    freundlich_n_dict = {"comp1": 5.1, "comp2": 1.2}
    bv_50_dict = {"comp1": 75000, "comp2": 90000}
    fix_vars = {
        "ebct": 210,
        "resin_density": 0.72,
        "resin_diam": 700e-6,
        "bed_porosity": 0.44,
        "loading_rate": 0.006,
        "number_columns": 6,
        "c_norm": 0.75,
        "number_columns_redundant": 2,
    }

    target_component = "comp1"
    flow_in = 0.04381
    number_traps = 5

    m = build_ix(
        flow_in=flow_in,
        target_component=target_component,
        solute_list=solute_list,
        reactive_ions=reactive_ions,
        fix_vars=fix_vars,
        c0_dict=c0_dict,
        mass_transfer_coeff_dict=mass_transfer_coeff_dict,
        freundlich_n_dict=freundlich_n_dict,
        bv_50_dict=bv_50_dict,
        diff_data_dict=diff_data_dict,
        number_traps=number_traps,
        charge_dict=charge_dict,
        dont_calc_effluent=False,
    )

    ix = m.fs.ix

    scale_model(m, scale_from_value=list())

    ix.c_trap_min["comp2"].set_value(1e-10)
    print(f"DOF: {degrees_of_freedom(m)}")
    set_scaling_factor(ix.c_traps, 1e3)
    set_scaling_factor(ix.traps, 1e6)
    set_scaling_factor(ix.c_norm_avg["comp2"], 1e6)
    try:
        ix.initialize()
    except:
        print_infeasible_constraints(m)
        print_variables_close_to_bounds(m)

    solver = get_solver()
    results = solver.solve(m)
    print(f"\ntermination = {results.solver.termination_condition.swapcase()}\n")

    plot_traps(m, x_plot="bv", run_curve=True, cmap="turbo")
    plt.show()


if __name__ == "__main__":
    main()
