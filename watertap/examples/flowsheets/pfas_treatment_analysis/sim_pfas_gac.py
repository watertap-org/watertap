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
#
###############################################################################

import pytest
import pyomo.environ as pyo
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from idaes.models.unit_models import Feed, Product
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

# set up solver
solver = get_solver()


def main():
    # set up simulation
    m = build()

    # solve simulation
    results = solver.solve(m, tee=True)
    assert pyo.check_optimal_termination(results)

    m.fs.feed.display()


def build():
    # build models
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["PFOS"], mw_data={"H2O": 0.018, "PFOS": 0.50013}
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.gac = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
    )

    # streams
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling property variables
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e7, index=("Liq", "PFOS")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-3, index=("Liq", "H2O")
    )

    # touch properties
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]

    # propagate scaling factors
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 170 / 60 / 60,  # 170 m**3/hr or 1 MGD
            ("conc_mass_phase_comp", ("Liq", "PFOS")): 1585.1e-9,  # 1585.1 ng/L
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # gac specifications
    # adsorption isotherm
    m.fs.gac.freund_k.fix()  # e-6 * (1e6 ** 0.8316))
    m.fs.gac.freund_ninv.fix()
    # gac particle specifications
    m.fs.gac.particle_porosity.fix()
    m.fs.gac.particle_dens_app.fix()
    m.fs.gac.particle_dia.fix()
    # adsorber bed specifications
    m.fs.gac.ebct.fix()  # seconds
    m.fs.gac.bed_voidage.fix()
    m.fs.gac.bed_length.fix()  # assumed
    # design spec
    m.fs.gac.conc_ratio_replace.fix()
    # parameters
    m.fs.gac.kf.fix()
    m.fs.gac.ds.fix()
    m.fs.gac.a0.fix()
    m.fs.gac.a1.fix()
    m.fs.gac.b0.fix()
    m.fs.gac.b1.fix()
    m.fs.gac.b2.fix()
    m.fs.gac.b3.fix()
    m.fs.gac.b4.fix()

    # initialization
    optarg = solver.options
    m.fs.feed.initialize(optarg=optarg)
    m.fs.gac.initialize(optarg=optarg)

    # check model
    assert_units_consistent(m)  # check that units are consistent
    print("Degrees of freedom:", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0

    return m


if __name__ == "__main__":
    main()
