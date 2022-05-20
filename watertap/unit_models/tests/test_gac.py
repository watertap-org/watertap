from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    TransformationFactory,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.generic_models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors
import idaes.core.util.scaling as iscale

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.gac import GAC


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["DCE"], "mw_data": {"H2O": 18e-3, "DCE": 98.96e-3}}
    )

    # unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    # m.fs.gac = GAC(default={"property_package": m.fs.properties})

    # connections
    # m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.gac.inlet)
    # TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e7, index=("Liq", "DCE")
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
            ("conc_mass_phase_comp", ("Liq", "DCE")): 2.32e-5,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    print("Feed degrees of freedom:", degrees_of_freedom(m.fs.feed))

    # TODO: Switch specifications to Feed model and arc to GAC, test GAC embedded touch of needed properties
    """
    m.fs.gac.treatwater.properties_in[0].conc_mass_phase_comp
    m.fs.gac.treatwater.properties_in[0].flow_vol_phase
    m.fs.gac.treatwater.properties_in[0].mw_comp
    
    m.fs.gac.treatwater.properties_in[0].temperature.fix(273.15 + 25)
    m.fs.gac.treatwater.properties_in[0].pressure.fix(101325)
    m.fs.gac.treatwater.properties_in[0].flow_vol_phase.fix(1)  # assumed basis in m**3/s
    m.fs.gac.treatwater.properties_in[0].conc_mass_phase_comp["Liq", "DCE"].fix(2.32e-5)


    reduction = 0.50*m.fs.gac.treatwater.properties_in[0].flow_mol_phase_comp['Liq', 'DCE'].value
    m.fs.gac.mass_solute_absorb[0, "Liq", "DCE"].fix(1e-8)
    """
    # TODO: Add variable scaling

    # solving
    assert_units_consistent(m)  # check that units are consistent
    # print("GAC model degrees of freedom:", degrees_of_freedom(m.fs.gac))
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    solver = SolverFactory("ipopt")
    solver.options = {"tol": 1e-8, "nlp_scaling_method": "user-scaling"}

    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal

    # m.fs.gac.display()

    # displays
    m.fs.feed.properties[0].conc_mass_phase_comp.display()


if __name__ == "__main__":
    main()
