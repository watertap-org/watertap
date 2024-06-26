# Imports from Pyomo, including "value" for getting the 
# value of Pyomo objects
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    value,
    units as pyunits,
)
# Imports from IDAES
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock
# Import function to get default solver
from watertap.core.solvers import get_solver
# Import function to check degrees of freedom
from idaes.core.util.model_statistics import degrees_of_freedom
# Import utility function for calculating scaling factors
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util import DiagnosticsToolbox

# Imports from WaterTAP
# Import NaCl property model
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
# Import RO model
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
        ConcentrationPolarizationType, MassTransferCoefficient)

# Create a Pyomo concrete model, flowsheet, and NaCl property parameter block.
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = NaClParameterBlock()

# Add an RO unit to the flowsheet.
m.fs.unit = ReverseOsmosis0D(
    property_package=m.fs.properties,
    concentration_polarization_type=ConcentrationPolarizationType.none,
    mass_transfer_coefficient=MassTransferCoefficient.none,
    has_pressure_change=False,
    )

m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(0.035)  # mass flow rate of NaCl (kg/s)
m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.965)   # mass flow rate of water (kg/s)
m.fs.unit.inlet.pressure[0].fix(50e5)                              # feed pressure (Pa)
m.fs.unit.inlet.temperature[0].fix(298.15)                         # feed temperature (K)
m.fs.unit.area.fix(50)                                             # membrane area (m^2)
m.fs.unit.A_comp.fix(4.2e-12)                                      # membrane water permeability (m/Pa/s)
m.fs.unit.B_comp.fix(3.5e-8)                                       # membrane salt permeability (m/s)
m.fs.unit.permeate.pressure[0].fix(101325)                         # permeate pressure (Pa)

# Set scaling factors for component mass flowrates.
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

# Set scaling factor for membrane area.
set_scaling_factor(m.fs.unit.area, 1e-2)

# Calculate scaling factors for all other variables.
calculate_scaling_factors(m)

# dt = DiagnosticsToolbox(m)
# dt.report_structural_issues()
# dt.display_underconstrained_set()
m.fs.unit.feed_side.display()
m.fs.unit.report()
m.fs.unit.initialize()

# Check that degrees of freedom = 0 before attempting simulation.
# This means that the performance of the flowsheet is completely
# determined by the system variables that were fixed above.
assert degrees_of_freedom(m) == 0
# Setup solver
solver = get_solver()
# Run simulation
simulation_results = solver.solve(m)
# Display report, reports include a small subset of the most important variables
m.fs.unit.report()
# Display all results, this shows all variables and constraints
# m.fs.unit.display()