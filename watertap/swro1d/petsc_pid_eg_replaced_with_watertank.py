import numpy as np
import matplotlib.pyplot as plt

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, MaterialBalanceType, EnergyBalanceType
from idaes.models_extra.power_generation.unit_models.watertank import WaterTank
import idaes_examples.mod.dae.petsc.pid_steam_tank as pid
from idaes.models.unit_models import Heater, Valve
from idaes.models.properties import iapws95
import pyomo.dae as pyodae
import idaes.core.solvers.petsc as petsc  # petsc utilities module
from idaes.core.util.initialization import propagate_state
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
)
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver, petsc
from idaes.core.util.math import smooth_max, smooth_min
from idaes.core.util import DiagnosticsToolbox

def _valve_pressure_flow_cb(b):
    """
    Callback for valve pressure-flow relation F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    umeta = b.config.property_package.get_metadata().get_derived_units

    b.Cv = pyo.Var(
        initialize=0.1,
        doc="Valve flow coefficent",
        units=umeta("amount") / umeta("time") / umeta("pressure"),
    )
    b.Cv.fix()

    b.flow_var = pyo.Reference(b.control_volume.properties_in[:].flow_mol)
    b.pressure_flow_equation_scale = lambda x: x**2

    @b.Constraint(b.flowsheet().time)
    def pressure_flow_equation(b2, t):
        Po = b2.control_volume.properties_out[t].pressure
        Pi = b2.control_volume.properties_in[t].pressure
        F = b2.control_volume.properties_in[t].flow_mol
        Cv = b2.Cv
        fun = b2.valve_function[t]
        return F**2 == Cv**2 * (Pi**2 - Po**2) * fun**2


def create_model(
    time_set=None,
    time_units=pyo.units.s,
    nfe=5,
    tee=False,
    calc_integ=True,
):
    """Create a test model and solver

    Args:
        time_set (list): The begining and end point of the time domain
        time_units (Pyomo Unit object): Units of time domain
        nfe (int): Number of finite elements argument for the DAE
            transformation.
        calc_integ (bool): If True, calculate in the initial condition for
            the integral term, else use a fixed variable (fs.ctrl.err_i0),
            False is the better option if you have a value from a previous
            time period

    Returns
        (tuple): (ConcreteModel, Solver)
    """
    fs_cfg = {"dynamic": True, "time_set": time_set, "time_units": time_units}
    model_name = "Steam Tank, Dynamic"

    if time_set is None:
        time_set = [0, 3]

    m = pyo.ConcreteModel(name=model_name)
    m.fs = FlowsheetBlock(**fs_cfg)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.LG
    )
    # Create the valve and tank models
    m.fs.valve_1 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.tank = WaterTank(
        property_package=m.fs.prop_water,
        material_balance_type=MaterialBalanceType.componentTotal,
        tank_type='rectangular_tank',
        has_heat_transfer='False',
        has_pressure_change='True'
    )
    m.fs.valve_2 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    # Add a controller
    m.fs.ctrl = PIDController(
        process_var=m.fs.tank.control_volume.properties_out[:].pressure,
        manipulated_var=m.fs.valve_1.valve_opening,
        calculate_initial_integral=calc_integ,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
        controller_type=ControllerType.PI,
    )
    # The control volume block doesn't assume the two phases are in equilibrium
    # by default, so I'll make that assumption here, I don't actually expect
    # liquid to form but who knows. The phase_fraction in the control volume is
    # volumetric phase fraction hence the densities.
    @m.fs.tank.Constraint(m.fs.time)
    def vol_frac_vap(b, t):
        return (
            b.control_volume.properties_out[t].phase_frac["Vap"]
            * b.control_volume.properties_out[t].dens_mol
            / b.control_volume.properties_out[t].dens_mol_phase["Vap"]
        ) == (b.control_volume.phase_fraction[t, "Vap"])

    # Connect the models
    m.fs.v1_to_tank = Arc(source=m.fs.valve_1.outlet, destination=m.fs.tank.inlet)
    m.fs.tank_to_v2 = Arc(source=m.fs.tank.outlet, destination=m.fs.valve_2.inlet)

    # Add the stream constraints and do the DAE transformation
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m.fs, nfe=nfe, wrt=m.fs.time, scheme="BACKWARD"
    )

    # Fix the derivative variables to zero at time 0 (steady state assumption)
    m.fs.fix_initial_conditions()

    # Fix the input variables
    m.fs.valve_1.inlet.enth_mol.fix(50000)
    m.fs.valve_1.inlet.pressure.fix(5e5)
    m.fs.valve_2.outlet.pressure.fix(101325)
    m.fs.valve_1.Cv.fix(0.001)
    m.fs.valve_2.Cv.fix(0.001)
    m.fs.valve_1.valve_opening.fix(1)
    m.fs.valve_2.valve_opening.fix(1)
    # m.fs.tank.heat_duty.fix(0)
    # m.fs.tank.tank_cross_sect_area.fix(150)
    m.fs.tank.tank_width.fix(1)
    m.fs.tank.tank_length.fix(1)
    # m.fs.tank.tank_level.fix(2.5)
    m.fs.tank.control_volume.volume.fix(2.0)
    # print(dir(m.fs.tank.pressure_change_eqn))
    # m.fs.tank.pressure_change_eqn.pprint()
    m.fs.tank.volume_eqn.pprint()

    # Fix controller settings
    m.fs.ctrl.gain_p.fix(1e-6)
    m.fs.ctrl.gain_i.fix(1e-5)
    # m.fs.ctrl.gain_d.fix(1e-6)
    # m.fs.ctrl.derivative_of_error[m.fs.time.first()].fix(0)
    m.fs.ctrl.setpoint.fix(3e5)
    m.fs.ctrl.mv_ref.fix(0)
    m.fs.ctrl.mv_lb = 0.0
    m.fs.ctrl.mv_ub = 1.0

    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize()
    propagate_state(m.fs.v1_to_tank)
    m.fs.tank.initialize()
    propagate_state(m.fs.tank_to_v2)
    # Can't specify both flow and outlet pressure so free the outlet pressure
    # for initialization and refix it after.  Inlet flow gets fixed in init, but
    # is unfixed for the final problem
    m.fs.valve_2.outlet.pressure.unfix()
    m.fs.valve_2.initialize()
    m.fs.valve_2.outlet.pressure.fix(101325)
    m.fs.valve_1.valve_opening.unfix()
    m.fs.valve_1.valve_opening[m.fs.time.first()].fix()
    # Return the model and solver
    return m

tf=12
m = create_model(
    time_set=[0, tf],
    nfe=1,
    calc_integ=True,
)

m.fs.tank.report()

# dt = DiagnosticsToolbox(m)
# dt.report_structural_issues()
# dt.display_components_with_inconsistent_units()
# dt.display_potential_evaluation_errors()

result = petsc.petsc_dae_by_time_element(
    m,
    time=m.fs.time,
    keepfiles=True,
    symbolic_solver_labels=True,
    ts_options={
        "--ts_type": "beuler",
        "--ts_dt": 0.01,
        "--ts_rtol": 1e-3,
        # "--ts_adapt_clip":"0.001,3600",
        # "--ksp_monitor":"",
        "--ts_adapt_dt_min": 1e-3,
        "--ts_adapt_dt_max": 3600,
        "--snes_type": "newtontr",
        # "--ts_max_reject": 200,
        # "--snes_monitor":"",
        "--ts_monitor": "",
        "--ts_save_trajectory": 1,
        "--ts_trajectory_type": "visualization",
        "--ts_max_snes_failures": 25,
        # "--show_cl":"",
        "-snes_max_it": 50,
        "-snes_rtol": 0,
        "-snes_stol": 0,
        "-snes_atol": 1e-6,
    },
    skip_initial=False,
    initial_solver="ipopt",
    initial_solver_options={
        "constr_viol_tol": 1e-8,
        "nlp_scaling_method": "user-scaling",
        "linear_solver": "ma27",
        "OF_ma57_automatic_scaling": "yes",
        "max_iter": 300,
        "tol": 1e-8,
        "halt_on_ampl_error": "yes",
    },
)
tj = result.trajectory  # trajectroy data
res = result.results  # solver status list

print(m.fs.tank.outlet.pressure)
print(dir(m.fs.tank.control_volume.properties_out[:]))

fig, axs = plt.subplots(5, 1, sharex=True)
axs[0].plot(tj.time, tj.get_vec(m.fs.tank.outlet.pressure[tf]), label='Tank outlet pressure')
axs[1].plot(tj.time, tj.get_vec(m.fs.tank.deltaP[tf]), label='Tank pressure drop')
axs[2].plot(tj.time, tj.get_vec(m.fs.tank.outlet.flow_mol[tf]), label='Tank outlet molar flow')
axs[3].plot(tj.time, tj.get_vec(m.fs.valve_1.valve_opening[tf]), label='Valve fraction open')
axs[4].plot(tj.time, tj.get_vec(m.fs.tank.tank_level[tf]), label='Tank level')
# axs[4].plot(tj.time, tj.get_vec(m.fs.tank.temperature[tf]), label='Tank inlet temperature')
axs[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
axs[0].set_xlim(0, tf)
axs[0].set_ylabel('Pressure (Pa)')
axs[1].set_ylabel('Pressure drop (Pa)')
axs[2].set_ylabel('Molar flow (mol/s)')
axs[3].set_ylabel('Valve fraction open')
axs[4].set_ylabel('Tank level (m)')
axs[4].set_ylim(tj.get_vec(m.fs.tank.tank_level[tf])[0]-0.01, tj.get_vec(m.fs.tank.tank_level[tf])[0]+0.01)
# axs[4].set_ylabel('Tank Temperature (K)')
# axs[0].grid(True)
axs[-1].set_xlabel('Time (s)')
fig.tight_layout()
plt.show()