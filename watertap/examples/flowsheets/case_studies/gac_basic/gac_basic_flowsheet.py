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

import math
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale

from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import (
    Feed,
    Product,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import assert_degrees_of_freedom

__author__ = "Hunter Barber"


def main():

    # testing ui functions
    m = build()
    initialize_model(m)
    res = solve_model(m)
    print("solver termination condition:", res.solver.termination_condition)
    m.fs.display()

    return m, res


def build():

    # blocks
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["solute"],
        mw_data={"H2O": 0.018, "solute": 0.08},
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.gac = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type=FilmTransferCoefficientType.fixed,
        surface_diffusion_coefficient_type=SurfaceDiffusionCoefficientType.fixed,
    )
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.adsorbed_removed = Product(property_package=m.fs.properties)

    # streams
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    m.fs.s02 = Arc(source=m.fs.gac.outlet, destination=m.fs.product.inlet)
    m.fs.s03 = Arc(source=m.fs.gac.adsorbed, destination=m.fs.adsorbed_removed.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # build costing blocks
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyo.units.USD_2021
    m.fs.gac.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"contactor_type": "pressure"},
    )

    # add flowsheet level blocks
    m.fs.costing.cost_process()
    treated_flow = m.fs.product.properties[0].flow_vol
    m.fs.costing.add_annual_water_production(treated_flow)
    m.fs.costing.add_LCOW(treated_flow)
    m.fs.costing.add_specific_energy_consumption(treated_flow)

    return m


def initialize_model(m, solver=None):

    if solver is None:
        solver = get_solver()

    # touch properties and default scaling
    water_sf = 10 ** -math.ceil(
        math.log10(abs(0.043813 * 1000 / m.fs.properties.mw_comp["H2O"].value))
    )
    solute_sf = 10 ** -math.ceil(
        math.log10(abs(0.043813 * 0.1 / m.fs.properties.mw_comp["solute"].value))
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", water_sf, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", solute_sf, index=("Liq", "solute")
    )
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]

    # automated scaling with unit model
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].unfix()
    m.fs.feed.properties[0].flow_mol_phase_comp["Liq", "solute"].unfix()
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.043813,
            ("conc_mass_phase_comp", ("Liq", "solute")): 0.1,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # gac specifications
    # performance parameters
    m.fs.gac.freund_k.fix(10)
    m.fs.gac.freund_ninv.fix(0.9)
    m.fs.gac.kf.fix(5e-5)
    m.fs.gac.ds.fix(2e-13)
    # gac particle specifications
    m.fs.gac.particle_dens_app.fix(750)
    m.fs.gac.particle_dia.fix(0.001)
    # adsorber bed specifications
    m.fs.gac.ebct.fix(600)
    m.fs.gac.bed_voidage.fix(0.4)
    m.fs.gac.bed_length.fix(6)
    # design spec
    m.fs.gac.conc_ratio_replace.fix(0.50)
    # parameters
    m.fs.gac.a0.fix(3.68421)
    m.fs.gac.a1.fix(13.1579)
    m.fs.gac.b0.fix(0.784576)
    m.fs.gac.b1.fix(0.239663)
    m.fs.gac.b2.fix(0.484422)
    m.fs.gac.b3.fix(0.003206)
    m.fs.gac.b4.fix(0.134987)

    # costing specifications
    m.fs.costing.gac_pressure.regen_frac.fix(0.7)
    m.fs.costing.gac_pressure.num_contactors_op.fix(1)
    m.fs.costing.gac_pressure.num_contactors_redundant.fix(1)

    # check model
    assert_units_consistent(m)
    assert_degrees_of_freedom(m, 0)

    # initialize models
    m.fs.feed.initialize()
    propagate_state(m.fs.s01)
    m.fs.gac.initialize()
    propagate_state(m.fs.s02)
    propagate_state(m.fs.s03)
    m.fs.product.initialize()
    m.fs.adsorbed_removed.initialize()

    # initialize costing
    m.fs.gac.costing.initialize()
    m.fs.costing.initialize()

    return m


def solve_model(m, solver=None):

    if solver is None:
        solver = get_solver()

    # check model
    assert_units_consistent(m)
    assert_degrees_of_freedom(m, 0)

    # solve model
    res = solver.solve(m)
    pyo.assert_optimal_termination(res)

    return res


if __name__ == "__main__":
    m, results = main()
