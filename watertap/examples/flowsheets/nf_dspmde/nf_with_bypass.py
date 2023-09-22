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

from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
)

from pyomo.environ import (
    units as pyunits,
    Var,
    assert_optimal_termination,
)

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
)
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType

from pyomo.environ import ConcreteModel, TransformationFactory

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.examples.flowsheets.nf_dspmde import nf

from watertap.costing import WaterTAPCosting


def main():
    solver = get_solver()
    m = build()
    initialize(m, solver)
    unfix_opt_vars(m)
    nf.add_objective(m)
    results = optimize(m, solver)
    assert_optimal_termination(results)
    print("Optimal cost", m.fs.costing.LCOW.value)
    print("Optimal NF pressure (Bar)", m.fs.NF.pump.outlet.pressure[0].value / 1e5)
    print("Optimal area (m2)", m.fs.NF.nfUnit.area.value)
    print(
        "Optimal nf recovery (%)",
        m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"].value * 100,
    )
    print("bypass (%)", m.fs.by_pass_splitter.split_fraction[0, "bypass"].value * 100)

    print("Feed hardness (mg/L as CaCO3)", m.fs.feed.properties[0].total_hardness.value)
    print(
        "Product hardness (mg/L as CaCO3)",
        m.fs.product.properties[0].total_hardness.value,
    )
    print(
        "Disposal hardness (mg/L as CaCO3)",
        m.fs.disposal.properties[0].total_hardness.value,
    )
    return m


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    default = nf.define_feed_comp()
    m.fs.properties = MCASParameterBlock(**default)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].conc_mass_phase_comp[...]
    m.fs.feed.properties[0].flow_mass_phase_comp[...]

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    nf.add_hardness_constraint(m.fs.product)

    m.fs.by_pass_splitter = Separator(
        property_package=m.fs.properties,
        outlet_list=["nf_stage", "bypass"],
    )
    # NF UNIT BLOCK
    m.fs.NF = FlowsheetBlock(dynamic=False)

    nf.build_nf_block(m, m.fs.NF)

    m.fs.total_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["bypass", "nf_stage"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.minimize,
    )
    m.fs.total_product_mixer.mixed_state[0.0].temperature.fix(293.15)
    m.fs.feed_to_splitter = Arc(
        source=m.fs.feed.outlet, destination=m.fs.by_pass_splitter.inlet
    )

    m.fs.splitter_to_nfUnit_feed = Arc(
        source=m.fs.by_pass_splitter.nf_stage, destination=m.fs.NF.feed.inlet
    )

    m.fs.splitter_to_mixer = Arc(
        source=m.fs.by_pass_splitter.bypass, destination=m.fs.total_product_mixer.bypass
    )

    m.fs.nfUnit_product_to_mixer = Arc(
        source=m.fs.NF.product.outlet,
        destination=m.fs.total_product_mixer.nf_stage,
    )

    m.fs.nfUnit_retentate_to_disposal = Arc(
        source=m.fs.NF.retentate.outlet,
        destination=m.fs.disposal.inlet,
    )
    m.fs.mixer_to_product = Arc(
        source=m.fs.total_product_mixer.outlet, destination=m.fs.product.inlet
    )
    m.fs.costing.disposal_cost = Var(
        initialize=0.1,
        bounds=(0, None),
        doc="disposal cost",
        units=pyunits.USD_2020 / pyunits.m**3,
    )
    m.fs.costing.disposal_cost.fix()
    m.fs.costing.add_defined_flow("disposal cost", m.fs.costing.disposal_cost)
    m.fs.costing.cost_flow(
        pyunits.convert(
            m.fs.disposal.properties[0].flow_vol_phase["Liq"],
            pyunits.m**3 / pyunits.s,
        ),
        "disposal cost",
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    iscale.set_scaling_factor(m.fs.costing.aggregate_flow_costs["disposal cost"], 1)
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def fix_init_vars(m):
    # apply defaults for normal NF init
    nf.fix_init_vars(m)
    # fix initial guess for splitter
    m.fs.by_pass_splitter.split_fraction[0, "bypass"].fix(0.9)
    m.fs.by_pass_splitter.split_fraction[0, "bypass"].setlb(0.05)
    m.fs.by_pass_splitter.split_fraction[0, "bypass"].setub(None)


def initialize(m, solver=None):
    if solver is None:
        solver = get_solver()
    # use standard nf default feed
    nf.set_default_feed(m, solver)
    fix_init_vars(m)

    init_system(m, solver)

    # solve box problem
    print("initialized, DOFs:", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    print("Solved box problem")


def init_system(m, solver):
    if solver is None:
        solver = get_solver()
    m.fs.feed.initialize(optarg=solver.options)

    propagate_state(m.fs.feed_to_splitter)

    m.fs.by_pass_splitter.mixed_state.initialize(optarg=solver.options)
    m.fs.by_pass_splitter.initialize(optarg=solver.options)
    propagate_state(m.fs.splitter_to_mixer)
    propagate_state(m.fs.splitter_to_nfUnit_feed)

    nf.init_nf_block(m.fs.NF, solver)

    propagate_state(m.fs.nfUnit_product_to_mixer)

    m.fs.total_product_mixer.mixed_state.initialize(optarg=solver.options)
    m.fs.total_product_mixer.initialize(optarg=solver.options)

    propagate_state(m.fs.nfUnit_retentate_to_disposal)
    propagate_state(m.fs.mixer_to_product)
    m.fs.NF.product.initialize(optarg=solver.options)
    m.fs.NF.retentate.initialize(optarg=solver.options)

    # seq = SequentialDecomposition(tear_solver="cbc")
    # seq.options.iterLim = 10
    #
    # def func_initialize(unit):
    #     unit.initialize(optarg=solver.options)
    #
    # seq.run(m, func_initialize)

    m.fs.costing.initialize()


def unfix_opt_vars(m):
    nf.unfix_opt_vars(m)
    m.fs.by_pass_splitter.split_fraction[0, "bypass"].unfix()


def optimize(m, solver=None):
    if solver is None:
        solver = get_solver()
    result = nf.optimize(m, solver)
    return result


if __name__ == "__main__":
    m = main()
