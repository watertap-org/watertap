#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    units as pyunits,
    Objective,
    assert_optimal_termination,
    Var,
    Constraint,
    NonNegativeReals,
    Reals,
    Param,
)


from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
)
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType

from watertap.unit_models.nanofiltration_DSPMDE_0D import (
    NanofiltrationDSPMDE0D,
)

from watertap.unit_models.pressure_changer import Pump

from pyomo.environ import ConcreteModel, TransformationFactory

import math

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
import idaes.core.util.scaling as iscale

from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting


def main():
    solver = get_solver()
    m = build()
    initialize(m, solver)
    unfix_opt_vars(m)
    optimize(m, solver)
    print("Optimal cost", m.fs.costing.LCOW.value)
    print("Optimal NF pressure (Bar)", m.fs.NF.pump.outlet.pressure[0].value / 1e5)
    print("Optimal area (m2)", m.fs.NF.nfUnit.area.value)
    print(
        "Optimal nf recovery (%)",
        m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"].value * 100,
    )

    print("Feed hardness (mg/L as CaCO3)", m.fs.feed.total_hardness.value)
    print("Product hardness (mg/L as CaCO3)", m.fs.product.total_hardness.value)
    print("Disposal hardness (mg/L as CaCO3)", m.fs.disposal.total_hardness.value)
    return m


def set_default_feed(m, solver):
    conc_mass_phase_comp = {
        "Ca_2+": 258 / 1000,
        "HCO3_-": 385 / 1000,
        "SO4_2-": 1011 / 1000,
        "Cl_-": 870.0 / 1000,
        "Na_+": 739 / 1000,
        "K_+": 9 / 1000,
        "Mg_2+": 90 / 1000,
    }
    set_NF_feed(
        m.fs,
        solver,
        flow_mass_h2o=1,
        conc_mass_phase_comp=conc_mass_phase_comp,
        feed_mass_frac=None,
        mole_frac=None,
    )


def define_feed_comp():
    default = {
        "solute_list": [
            "Ca_2+",
            "SO4_2-",
            "HCO3_-",
            "Na_+",
            "Cl_-",
            "K_+",
            "Mg_2+",
        ],
        "diffusivity_data": {
            ("Liq", "Ca_2+"): 9.2e-10,
            ("Liq", "SO4_2-"): 1.06e-9,
            ("Liq", "HCO3_-"): 1.19e-9,
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
            ("Liq", "K_+"): 1.957e-9,
            ("Liq", "Mg_2+"): 0.706e-9,
        },
        "mw_data": {
            "H2O": 18e-3,
            "Ca_2+": 40e-3,
            "HCO3_-": 61.0168e-3,
            "SO4_2-": 96e-3,
            "Na_+": 23e-3,
            "Cl_-": 35e-3,
            "K_+": 22.989769e-3,
            "Mg_2+": 24.305e-3,
        },
        "stokes_radius_data": {
            "Ca_2+": 0.309e-9,
            "HCO3_-": 2.06e-10,
            "SO4_2-": 0.230e-9,
            "Cl_-": 0.121e-9,
            "Na_+": 0.184e-9,
            "K_+": 0.125e-9,
            "Mg_2+": 0.347e-9,
        },
        "charge": {
            "Ca_2+": 2,
            "HCO3_-": -1,
            "SO4_2-": -2,
            "Na_+": 1,
            "Cl_-": -1,
            "K_+": 1,
            "Mg_2+": 2,
        },
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
    }
    return default


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = WaterTAPCosting()
    default = define_feed_comp()
    m.fs.properties = MCASParameterBlock(**default)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.feed.properties[0].flow_mass_phase_comp[...]
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    add_hardness_constraint(m.fs.feed)
    add_hardness_constraint(m.fs.product)
    add_hardness_constraint(m.fs.disposal)
    m.fs.NF = FlowsheetBlock(dynamic=False)
    build_nf_block(m, m.fs.NF)
    m.fs.feed_to_nfUnit_feed = Arc(
        source=m.fs.feed.outlet, destination=m.fs.NF.feed.inlet
    )
    m.fs.nfUnit_retentate_to_disposal = Arc(
        source=m.fs.NF.retentate.outlet,
        destination=m.fs.disposal.inlet,
    )
    m.fs.nfUnit_product_to_product = Arc(
        source=m.fs.NF.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)

    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def build_nf_block(m, blk):
    # setting up state junction
    blk.feed = StateJunction(property_package=m.fs.properties)
    m.fs.NF.product = StateJunction(property_package=m.fs.properties)
    blk.retentate = StateJunction(property_package=m.fs.properties)

    blk.pump = Pump(property_package=m.fs.properties)
    blk.pump.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    # nfunit
    blk.nfUnit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)
    blk.nfUnit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.inlet)
    blk.pump_to_nf = Arc(source=blk.pump.outlet, destination=blk.nfUnit.inlet)
    blk.nf_to_product = Arc(source=blk.nfUnit.permeate, destination=blk.product.inlet)
    blk.nf_to_retentate = Arc(
        source=blk.nfUnit.retentate, destination=blk.retentate.inlet
    )


def fix_init_vars(m):
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(293.15)
    # splitter guiess

    # pump operation
    m.fs.NF.pump.outlet.pressure[0].fix(3.0 * 1e5)
    m.fs.NF.pump.efficiency_pump[0].fix(0.75)
    iscale.set_scaling_factor(m.fs.NF.pump.control_volume.work, 1e-4)
    # NF unit operation init values
    m.fs.NF.nfUnit.area.fix(50)
    m.fs.NF.nfUnit.spacer_porosity.fix(0.85)
    m.fs.NF.nfUnit.spacer_mixing_efficiency.fix()
    m.fs.NF.nfUnit.spacer_mixing_length.fix()
    m.fs.NF.nfUnit.channel_height.fix(1e-3)
    m.fs.NF.nfUnit.velocity[0, 0].fix(0.25)
    m.fs.NF.nfUnit.mixed_permeate[0].pressure.fix(101325)
    # NF membrane props for NF270
    m.fs.NF.nfUnit.radius_pore.fix(0.5e-9)
    m.fs.NF.nfUnit.membrane_thickness_effective.fix(8.598945196055952e-07)
    m.fs.NF.nfUnit.membrane_charge_density.fix(-680)
    m.fs.NF.nfUnit.dielectric_constant_pore.fix(41.3)


def unfix_opt_vars(m):
    m.fs.NF.pump.outlet.pressure[0].unfix()
    m.fs.NF.nfUnit.area.unfix()
    m.fs.product.max_hardness.fix(200)


# m.fs.product.max_hardness.fix(500)


def add_objective(m):
    if m.find_component("fs.LCOW_objective") is None:
        m.fs.LCOW_objective = Objective(expr=m.fs.costing.LCOW)


def optimize(m, solver, **kwargs):
    add_objective(m)
    print("Optimizing with {} DOFS".format(degrees_of_freedom(m)))
    result = solver.solve(m, tee=True)
    return result


def initialize(m, solver=None, **kwargs):
    set_default_feed(m, solver)
    fix_init_vars(m)
    init_system(m, solver)
    # solve box problem
    print("initalized, DOFS:", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    solver.solve(m, tee=True)
    print("Solved box problem")


def init_system(m, solver):
    if solver == None:
        solver = get_solver()
    m.fs.feed.initialize(optarg=solver.options)

    propagate_state(m.fs.feed_to_nfUnit_feed)

    init_nf_block(m, m.fs.NF, solver)

    propagate_state(m.fs.nfUnit_retentate_to_disposal)
    propagate_state(m.fs.nfUnit_product_to_product)
    m.fs.NF.product.initialize(optarg=solver.options)
    m.fs.NF.retentate.initialize(optarg=solver.options)


def init_nf_block(m, blk, solver):
    blk.feed.initialize(optarg=solver.options)
    propagate_state(blk.feed_to_pump)
    blk.pump.initialize(optarg=solver.options)
    propagate_state(blk.pump_to_nf)
    blk.nfUnit.initialize(
        optarg=solver.options,
        automate_rescale=False,
        initialize_guess={
            "solvent_recovery": 0.01,
            "solute_recovery": 0.01,
            "cp_modulus": 1,
        },
    )
    propagate_state(blk.nf_to_product)
    propagate_state(blk.nf_to_retentate)
    blk.product.initialize(optarg=solver.options)
    blk.retentate.initialize(optarg=solver.options)


def set_NF_feed(
    blk,
    solver,
    flow_mass_h2o,
    conc_mass_phase_comp=None,
    feed_mass_frac=None,
    mole_frac=None,
):
    mass_flow_in = flow_mass_h2o * pyunits.kg / pyunits.s
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_h2o)

    if feed_mass_frac is not None:
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / blk.feed.properties[0].mw_comp[ion]
            )
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / blk.feed.properties[0].mw_comp["H2O"]
        )
        blk.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)

    if mole_frac is not None:
        for ion, x in mole_frac.items():
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(x)
    if conc_mass_phase_comp is not None:
        for ion, x in conc_mass_phase_comp.items():
            blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].fix(x)
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].unfix()
        result = solver.solve(blk.feed)
        blk.feed.properties[0].conc_mass_phase_comp["Liq", "H2O"].fix()
        for ion, x in conc_mass_phase_comp.items():
            blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].unfix()
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix()
            blk.feed.properties[0].flow_mass_phase_comp["Liq", ion].unfix()
        blk.feed.properties[0].conc_mass_phase_comp["Liq", "H2O"].unfix()
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        blk.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix()
    blk.feed.properties[0].display()
    set_NF_feed_scaling(blk)
    blk.feed.properties[0].assert_electroneutrality(
        defined_state=True,
        adjust_by_ion="Cl_-",
        get_property="flow_mol_phase_comp",
    )

    blk.feed.properties[0].temperature.fix(298.15)
    iscale.calculate_scaling_factors(blk)
    # switching to concetration for ease of adjusting in UI
    for ion, x in conc_mass_phase_comp.items():
        blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].unfix()
        blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix()


def calc_scale(value):
    return -1 * math.floor(math.log(value, 10))


def set_NF_feed_scaling(blk):
    _add = 0
    for index in blk.feed.properties[0].flow_mol_phase_comp:
        scale = calc_scale(blk.feed.properties[0].flow_mol_phase_comp[index].value)
        print(f"{index} flow_mol_phase_comp scaling factor = {10**(scale+_add)}")
        blk.properties.set_default_scaling(
            "flow_mol_phase_comp", 10 ** (scale + 0), index=index
        )


def add_hardness_constraint(stream):
    comp_list = stream.properties[0.0].component_list
    # print("comp_list", comp_list)
    stream.total_hardness = Var(
        initialize=0.5,
        bounds=(0, 10000),
        domain=NonNegativeReals,
        units=pyunits.mg / pyunits.kg,
        doc="Water Hardness as CO3",
    )
    stream.total_hardness.unfix()
    iscale.set_scaling_factor(stream.total_hardness, 1)
    stream.Ca_hardness = Var(
        initialize=0.5,
        bounds=(0, 10000),
        domain=NonNegativeReals,
        units=pyunits.mg / pyunits.kg,
        doc="Water Hardness as CaCO3",
    )
    stream.Ca_hardness.unfix()
    iscale.set_scaling_factor(stream.Ca_hardness, 1)
    stream.Mg_hardness = Var(
        initialize=0,
        bounds=(0, 10000),
        domain=NonNegativeReals,
        units=pyunits.mg / pyunits.kg,
        doc="Water Hardness as MgCO3",
    )
    stream.Mg_hardness.unfix()
    iscale.set_scaling_factor(stream.Mg_hardness, 1)
    stream.CaCO3_mw = Param(
        initialize=100.1,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.mol,
        doc="System Water Recovery",
    )
    iscale.set_scaling_factor(stream.CaCO3_mw, 0.01)
    stream.max_hardness = Var(
        initialize=10000,
        domain=NonNegativeReals,
        units=pyunits.g / pyunits.mol,
        doc="System Water Recovery",
    )
    stream.max_hardness.unfix()
    if "Ca_2+" in comp_list:
        stream.Ca_hardness_constraint = Constraint(
            expr=stream.Ca_hardness
            == (
                stream.properties[0].flow_mol_phase_comp["Liq", "Ca_2+"]
                * pyunits.convert(stream.CaCO3_mw, pyunits.mg / pyunits.mol)
            )
            / stream.properties[0].flow_vol_phase["Liq"]
            / 1000
        )
    else:
        stream.Ca_hardness.fix(0)
    if "Mg_2+" in comp_list:
        stream.Mg_hardness_constraint = Constraint(
            expr=stream.Mg_hardness
            == (
                stream.properties[0].flow_mol_phase_comp["Liq", "Mg_2+"]
                * pyunits.convert(stream.CaCO3_mw, pyunits.mg / pyunits.mol)
                / pyunits.mol
            )
            / stream.properties[0].flow_vol_phase["Liq"]
            / 1000
        )
    else:
        stream.Mg_hardness.fix(0)

    stream.total_hardness_constraint = Constraint(
        expr=stream.total_hardness == stream.Ca_hardness + stream.Mg_hardness
    )

    stream.max_hardness_constraint = Constraint(
        expr=stream.total_hardness <= stream.max_hardness
    )


if __name__ == "__main__":
    main()
