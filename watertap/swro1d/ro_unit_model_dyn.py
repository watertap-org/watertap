# Imports from Pyomo, including "value" for getting the 
# value of Pyomo objects
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    value,
    units as pyunits,
)
import pyomo.environ as pyo
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

import idaes.logger as idaeslog
from idaes.core.solvers import petsc
import idaes.core.util.scaling as iscale

# Imports from WaterTAP
# Import NaCl property model
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
# Import RO model
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
        ConcentrationPolarizationType, MassTransferCoefficient)

# Create a Pyomo concrete model, flowsheet, and NaCl property parameter block.
m = ConcreteModel()
m.fs = FlowsheetBlock(
    dynamic=True,
    time_set=[0, 3],
    time_units=pyunits.s
)
# m.fs = FlowsheetBlock(dynamic=False)
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

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
# dt.display_underconstrained_set()
m.fs.unit.feed_side.properties_interface[0.0,0.0].temperature.fix(298.15) # K
m.fs.unit.initialize()

# Check that degrees of freedom = 0 before attempting simulation.
# This means that the performance of the flowsheet is completely
# determined by the system variables that were fixed above.
# assert degrees_of_freedom(m) == 0
# Setup solver
solver = get_solver()
# Run simulation
simulation_results = solver.solve(m)
# Display report, reports include a small subset of the most important variables
m.fs.unit.display()
# Display all results, this shows all variables and constraints
# m.fs.unit.display()

m.fs.unit.fix_initial_conditions()

time_nfe = len(m.fs.time) - 1

pyo.TransformationFactory("dae.finite_difference").apply_to(
    m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
)
scaling_log = idaeslog.getLogger("idaes.core.util.scaling")
scaling_log.setLevel(idaeslog.ERROR)
iscale.calculate_scaling_factors(m)

idaeslog.solver_log.tee = True
results = petsc.petsc_dae_by_time_element(
    m,
    time=m.fs.time,
    keepfiles=True,
    symbolic_solver_labels=True,
    ts_options={
        "--ts_type": "beuler",
        # "-ts_arkimex_type": "1bee",
        "--ts_dt": 0.1,
        "--ts_rtol": 1e-3,
        # "--ts_adapt_clip":"0.001,3600",
        # "--ksp_monitor":"",
        "--ts_adapt_dt_min": 1e-3,
        "--ts_adapt_dt_max": 3600,
        "--snes_type": "newtontr",
        # "--ts_max_reject": 200,
        "--ts_monitor": "",
        "-ts_adapt_monitor": "",
        # "--snes_monitor":"",
        "-snes_converged_reason": "",
        # "-ksp_monitor_true_residual": "",
        # "-ksp_converged_reason": "",
        # "-snes_test_jacobian": "",
        "snes_grid_sequence": "",
        "-pc_type": "lu",
        # "-mat_view": "",
        "--ts_save_trajectory": 1,
        "--ts_trajectory_type": "visualization",
        "--ts_max_snes_failures": 25,
        # "--show_cl":"",
        "-snes_max_it": 50,
        "-snes_rtol": 0,
        "-snes_stol": 0,
        "-snes_atol": 1e-6,
    },
    skip_initial=True,
    initial_solver="ipopt",
    initial_solver_options={
        "constr_viol_tol": 1e-8,
        "nlp_scaling_method": "user-scaling",
        "linear_solver": "mumps",
        "OF_ma57_automatic_scaling": "yes",
        "max_iter": 300,
        "tol": 1e-8,
        "halt_on_ampl_error": "no",
    },
)
for result in results.results:
    pyo.assert_optimal_termination(result)

print(pyo.value(m.fs.wtank.control_volume.material_accumulation[:, :, :]))