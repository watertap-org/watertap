import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    Expression,
    TransformationFactory,
    value,
    assert_optimal_termination,
    check_optimal_termination,
    log,
    log10,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.scaling import *
from idaes.core.util.testing import initialization_tester
from idaes.models.unit_models import Mixer, Product, Feed

from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    DSPMDEStateBlock,
)
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
)
from watertap.core.util.infeasible import *
from watertap.costing import WaterTAPCosting


def ix_build(ions, target_ion=None, hazardous_waste=False, regenerant=None):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    ion_prop = get_ion_props([ions])

    m.fs.properties = DSPMDEParameterBlock(**ion_prop)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.regen = Product(property_package=m.fs.properties)

    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.product.properties[0].conc_mass_phase_comp[...]
    m.fs.regen.properties[0].conc_mass_phase_comp[...]

    m.fs.feed.properties[0].flow_vol_phase[...]
    m.fs.product.properties[0].flow_vol_phase[...]
    m.fs.regen.properties[0].flow_vol_phase[...]

    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "hazardous_waste": hazardous_waste,
    }

    if regenerant is not None:
        ix_config["regenerant"] = regenerant

    m.fs.ion_exchange = ix = IonExchange0D(default=ix_config)

    m.fs.feed_to_ix = Arc(source=m.fs.feed.outlet, destination=ix.inlet)
    m.fs.ix_to_product = Arc(source=ix.outlet, destination=m.fs.product.inlet)
    m.fs.ix_to_regen = Arc(source=ix.regen, destination=m.fs.regen.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)
    # ix.display()
    return m


def ix_scaling(m, sf=1e4, est_removal=0.99, est_recov=0.95, sf_mass=1e4):

    ix = m.fs.ion_exchange
    props = m.fs.properties
    prop_in = ix.properties_in[0]
    prop_out = ix.properties_out[0]
    prop_regen = ix.properties_regen[0]
    target_ion = ix.config.target_ion
    ions = ix.config.property_package.ion_set
    sf_flow = 1 / prop_in.flow_mol_phase_comp["Liq", "H2O"].value
    # sf_flow_out = 1 / (
    #     prop_in.flow_mol_phase_comp["Liq", "H2O"].value * (1 - est_removal)
    # )
    set_scaling_factor(prop_in.flow_mol_phase_comp["Liq", "H2O"], sf_flow)
    set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", "H2O"], sf_flow)
    set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", "H2O"], sf)

    for ion in ions:
        props.set_default_scaling("flow_mol_phase_comp", 10, index=("Liq", ion))
        set_scaling_factor(
            prop_in.flow_mol_phase_comp["Liq", ion],
            1 / prop_in.flow_mol_phase_comp["Liq", ion].value,
        )
        if ion == target_ion:
            sf_targ = 1 / (
                prop_in.flow_mol_phase_comp["Liq", ion].value * (1 - est_removal)
            )
            set_scaling_factor(prop_out.flow_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_equiv_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mol_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_out.conc_mass_phase_comp["Liq", ion], sf_targ)
            set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", ion], sf)

        else:
            set_scaling_factor(
                prop_out.flow_mol_phase_comp["Liq", ion],
                1 / (prop_in.flow_mol_phase_comp["Liq", ion].value),
            )
            set_scaling_factor(prop_regen.conc_equiv_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_regen.flow_mol_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_regen.conc_mol_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_regen.conc_mass_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_regen.mass_frac_phase_comp["Liq", ion], sf)
            set_scaling_factor(prop_regen.flow_mass_phase_comp["Liq", ion], sf)
            # set_scaling_factor(ix.mass_out[ion], sf_mass_in)
    if hasattr(m.fs, "costing"):
        regen_flow_var = getattr(m.fs.costing, f"aggregate_flow_{ix.regen_chem}")
        set_scaling_factor(regen_flow_var, 1 / sf)
        set_scaling_factor(ix.costing.capital_cost, 1e-6)
        set_scaling_factor(ix.costing.capital_cost_vessel, 1e-6)
        set_scaling_factor(ix.costing.capital_cost_resin, 1e-6)
        set_scaling_factor(ix.costing.capital_cost_backwash_tank, 1e-6)
        set_scaling_factor(ix.costing.capital_cost_regen_tank_constraint, 1e-5)

    props.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    
    
    calculate_scaling_factors(m)

    set_scaling_factor(ix.mass_in[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_removed[target_ion], 1 / sf_mass)
    set_scaling_factor(ix.mass_out[target_ion], sf_mass)
    # return m


def get_ion_props(ions):
    if not isinstance(ions, list):
        ions = [ions]
    diff_data = {
        "Na_+": 1.33e-9,
        "Ca_2+": 9.2e-10,
        "Cl_-": 2.03e-9,
        "Mg_2+": 0.706e-9,
        "SO4_2-": 1.06e-9,
        "PFAS_-": 0.49e-9,
        "Hardness_2+": 0.706e-9,
    }
    mw_data = {
        "Na_+": 23e-3,
        "Ca_2+": 40e-3,
        "Cl_-": 35e-3,
        "Mg_2+": 24e-3,
        "SO4_2-": 96e-3,
        "PFAS_-": 414.1e-3,
        "Hardness_2+": 100.0869e-3,
    }
    charge_data = {
        "Na_+": 1,
        "Ca_2+": 2,
        "Cl_-": -1,
        "Mg_2+": 2,
        "SO4_2-": -2,
        "PFAS_-": -1,
        "Hardness_2+": 2,
    }
    ix_in = {
        "solute_list": [],
        "diffusivity_data": {},
        "mw_data": {"H2O": 18e-3},
        "charge": {},
    }
    for ion in ions:
        ix_in["solute_list"].append(ion)
        ix_in["diffusivity_data"][("Liq", ion)] = diff_data[ion]
        ix_in["mw_data"][ion] = mw_data[ion]
        ix_in["charge"][ion] = charge_data[ion]
    return ix_in

def set_operating_conditions(m, feed_mass_frac={}, mass_flow_in=43.81264, solver=None):
    if solver is None:
        solver = get_solver()
    ix = m.fs.ion_exchange
    feed = m.fs.feed.properties[0]
    props = m.fs.properties
    mass_flow_in = mass_flow_in * (pyunits.kg / pyunits.s)

    for ion, x in feed_mass_frac.items():
        mol_flow = (x * (pyunits.kg / pyunits.kg) * mass_flow_in / ix.config.property_package.mw_comp[ion])
        feed.flow_mol_phase_comp["Liq", ion].fix(mol_flow)

    h2o_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    h2o_mol_flow = (h2o_mass_frac *  (pyunits.kg / pyunits.kg) * mass_flow_in / ix.config.property_package.mw_comp["H2O"])
    feed.flow_mol_phase_comp["Liq", "H2O"].fix(h2o_mol_flow)

    

target_ion = "Ca_2+"
target_ion_mass_frac = 2.5e-4
inert_ion_mass_frac = None
ions = [target_ion]
feed_mass_frac = {}
for ion in ions:
    if ion == target_ion:
        feed_mass_frac[ion] = target_ion_mass_frac
    else:
        feed_mass_frac[ion] = inert_ion_mass_frac
m = ix_build(target_ion, target_ion=target_ion)
ix_scaling(m)
set_operating_conditions(m, feed_mass_frac=feed_mass_frac)
