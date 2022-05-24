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

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.gac import GAC


def main():

    # set solver
    solver = SolverFactory("ipopt")
    solver.options = {"tol": 1e-8, "nlp_scaling_method": "user-scaling"}

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["DCE"], "mw_data": {"H2O": 18, "DCE": 98.96}}
    )  # TODO: Check MW units

    # unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.gac = GAC(default={"property_package": m.fs.properties})
    m.fs.prod = Product(default={"property_package": m.fs.properties})
    m.fs.remo = Product(default={"property_package": m.fs.properties})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    m.fs.s02 = Arc(source=m.fs.gac.outlet, destination=m.fs.prod.inlet)
    m.fs.s03 = Arc(source=m.fs.gac.spent_gac, destination=m.fs.remo.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "DCE")
    )
    # touch properties
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mw_comp
    # propagate scaling factors
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1e-3,  # assumed basis in m**3/s
            ("conc_mass_phase_comp", ("Liq", "DCE")): 23.2,  # 2.32e-5,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # TODO: Switch specifications to Feed model and arc to GAC, test GAC embedded touch of needed properties
    m.fs.gac.solute_absorb[0, "Liq", "DCE"].fix(1e-4)

    # TODO: Add variable scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "DCE"))
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

    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    results = solver.solve(m, tee=False)
    # assert results.solver.termination_condition == TerminationCondition.optimal

    # displays
    # m.fs.feed.properties[0].conc_mass_phase_comp.display()
    st = create_stream_table_dataframe(
        {"In": m.fs.s01, "Out": m.fs.s02, "Removed": m.fs.s03}
    )
    print(stream_table_dataframe_to_string(st))
    m.fs.gac.treatwater.mass_transfer_term.display()
    m.fs.gac.solute_absorb.display()


if __name__ == "__main__":
    main()
