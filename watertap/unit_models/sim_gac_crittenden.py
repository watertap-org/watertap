###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# 'https://github.com/watertap-org/watertap/'
#
###############################################################################

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    TransformationFactory,
    Objective,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from idaes.core.util.model_diagnostics import DegeneracyHunter

from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from watertap.unit_models.gac import GAC


def main():

    # set solver
    solver = SolverFactory("ipopt")
    solver.options = {"tol": 1e-8, "nlp_scaling_method": "user-scaling"}

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # TODO: Determine whether new property pack is needed to handle lower than 1e-8 ---> 1e-12
    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["TCE"], "mw_data": {"H2O": 18e-3, "TCE": 131.4e-3}}
    )

    # unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.gac = GAC(
        default={
            "property_package": m.fs.properties,
            "film_transfer_rate_type": "calculated",
            "surface_diffusion_coefficient_type": "fixed",
        }
    )
    m.fs.prod = Product(default={"property_package": m.fs.properties})
    m.fs.remo = Product(default={"property_package": m.fs.properties})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    m.fs.s02 = Arc(source=m.fs.gac.outlet, destination=m.fs.prod.inlet)
    m.fs.s03 = Arc(source=m.fs.gac.spent_gac, destination=m.fs.remo.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "TCE")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-6, index=("Liq", "H2O")
    )
    # touch properties
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].dens_mass_phase["Liq"]
    m.fs.feed.properties[0].mw_comp
    # propagate scaling factors
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.89 * 60,  # assumed basis in m**3/s
            ("conc_mass_phase_comp", ("Liq", "TCE")): 500e-6,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # gac specifications
    # trial problem from Crittenden
    # --------------------------------------------------------------------
    # specify two of three of these
    m.fs.gac.target_removal_frac["TCE"].fix(0.5)
    # m.fs.gac.replace_removal_frac["TCE"].fix(0.8)
    m.fs.gac.replace_gac_saturation_frac.fix(0.99)
    # --------------------------------------------------------------------
    m.fs.gac.freund_k.fix(1062e-6 * (1e6**0.48))
    m.fs.gac.freund_ninv.fix(0.48)
    m.fs.gac.ebct.fix(10 * 60)
    m.fs.gac.eps_bed.fix(0.44)
    m.fs.gac.particle_dens_app.fix(803.4)
    m.fs.gac.particle_dp.fix(0.001026)
    # m.fs.gac.kf.fix(3.73e-5)
    m.fs.gac.ds.fix(1.24e-14)
    m.fs.gac.velocity_sup.fix(5 / 3600 / 0.44)
    m.fs.gac.molal_volume.fix(9.81e-5)
    # m.fs.gac.shape_correction_factor.fix(1.5)
    # TODO: Determine whether to embed tabulated data for coefficients
    m.fs.gac.a0.fix(0.8)
    m.fs.gac.a1.fix(0)
    m.fs.gac.b0.fix(0.023)
    m.fs.gac.b1.fix(0.793673)
    m.fs.gac.b2.fix(0.039324)
    m.fs.gac.b3.fix(0.009326)
    m.fs.gac.b4.fix(0.08275)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-6, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "TCE")
    )
    calculate_scaling_factors(m)

    # initialization
    optarg = solver.options
    m.fs.feed.initialize(optarg=optarg)
    m.fs.gac.initialize(optarg=optarg)
    m.fs.prod.initialize(optarg=optarg)
    m.fs.remo.initialize(optarg=optarg)

    # solving
    assert_units_consistent(m)  # check that units are consistent
    print("Degrees of freedom:", degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    """
    # initialization testing
    # solver.options['max_iter'] = 0
    m.obj = Objective(expr=0)  # dummy for DegeneracyHunter
    results = solver.solve(m, tee=False)
    dh = DegeneracyHunter(m, solver=SolverFactory("cbc"))
    # dh.check_residuals(tol=0.01)
    # dh.check_variable_bounds(tol=1e-3)
    m.fs.gac.display()
    """
    # run simulation
    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal
    # displays
    st = create_stream_table_dataframe(
        {"In": m.fs.s01, "Out": m.fs.s02, "Removed": m.fs.s03}
    )
    print(stream_table_dataframe_to_string(st))
    m.fs.gac.display()
    m.fs.gac.molecular_diffusion_coefficient.display()
    m.fs.gac.re.display()
    m.fs.gac.sc.display()
    m.fs.gac.kf.display()
    # """


if __name__ == "__main__":
    main()
