from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.gac import GAC

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = DSPMDEParameterBlock(
    default={"solute_list": ["DCE"], "mw_data": {"H2O": 18e-3, "DCE": 98.96e-3}}
)
m.fs.gac = GAC(default={"property_package": m.fs.properties})


# TODO: Switch specifications to Feed model and arc to GAC, test GAC embedded touch of needed properties
m.fs.gac.treatwater.properties_in[0].conc_mass_phase_comp
m.fs.gac.treatwater.properties_in[0].flow_vol_phase
m.fs.gac.treatwater.properties_in[0].mw_comp

m.fs.gac.treatwater.properties_in[0].temperature.fix(273.15 + 25)
m.fs.gac.treatwater.properties_in[0].pressure.fix(101325)
m.fs.gac.treatwater.properties_in[0].flow_vol_phase.fix(1)  # assumed basis in m**3/s
m.fs.gac.treatwater.properties_in[0].conc_mass_phase_comp["Liq", "DCE"].fix(2.32e-5)

m.fs.gac.mass_solute_absorb[0, "Liq", "DCE"].fix(1e-8)

# solving
assert_units_consistent(m)  # check that units are consistent
print("GAC model degrees of freedom:", degrees_of_freedom(m.fs.gac))
assert (
    degrees_of_freedom(m) == 0
)  # check that the degrees of freedom are what we expect

solver = SolverFactory("ipopt")
solver.options = {"tol": 1e-8, "nlp_scaling_method": "user-scaling"}

results = solver.solve(m, tee=False)
assert results.solver.termination_condition == TerminationCondition.optimal

m.fs.gac.display()
