import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from idaes.core.solvers import petsc
import idaes.core.util.scaling as iscale

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from watertap.core.solvers import get_solver
from pyomo.network import Arc

from pyomo.environ import (
    ConcreteModel,
    units as pyunits,
    TransformationFactory,
    assert_optimal_termination,
    value,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    StateBlock,
)

from watertap.core import (
    MembraneChannel0DBlock,
    FrictionFactor,
    ModuleType,
)

from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
)

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from idaes.models.unit_models.mixer import MomentumMixingType

from watertap.unit_models.pressure_changer import Pump

import watertap.property_models.NaCl_prop_pack as props

def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    set_operating_conditions(m)
    initialize_system(m)

    assert_optimal_termination(solve(m, solver=solver))

    # optimize and display
    solve_dynamic(m)

    display_system(m)
    display_design(m)

    return m

def build():
    # flowsheet set up
    m = ConcreteModel()
    num_time_points = 2
    time_set = np.linspace(0, 1400, num_time_points+1)
    m.fs = FlowsheetBlock(
        dynamic=False,
        # time_set=list(time_set),
        # time_units=pyunits.s,
    )
    m.fs.properties = props.NaClParameterBlock()

    # Control volume flow blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets=2,
        momentum_mixing_type=MomentumMixingType.none,
    )

    # --- Feed and recycle pumps ---
    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P2 = Pump(property_package=m.fs.properties)

    # --- Reverse Osmosis Block ---
    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        module_type="spiral_wound",
        has_full_reporting=True,
    )

    # connect unit models
    m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
    m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
    m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

    m.fs.RO_permeate_to_product = Arc(
        source=m.fs.RO.permeate, destination=m.fs.product.inlet
    )
    m.fs.RO_retentate_to_P2 = Arc(
        source=m.fs.RO.retentate, destination=m.fs.P2.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    if m.fs.config.dynamic:
        time_nfe = len(m.fs.time) - 1
        TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
        )

    # set unit model values
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-1)

    return m

def set_operating_conditions(
    m,
    flow_vol=0.9,
    salt_mass_conc=4.0,
    solver=None,
):

    if solver is None:
        solver = get_solver()
    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(170 * pyunits.psi)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix((273.15 + 21.8) * pyunits.degK)  # feed temperature [K]

    # scaling
    # set default property values

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1000 * flow_vol, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-3 / flow_vol / salt_mass_conc, index=("Liq", "NaCl")
    )

    # m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(1.8/60)
    # m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix(1.5e-3)

    # set scaling factors
    iscale.set_scaling_factor(
        m.fs.P1.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
    )
    iscale.set_scaling_factor(
        m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"], 1
    )
    iscale.set_scaling_factor(m.fs.P1.work_fluid[0], 1)
    iscale.set_scaling_factor(m.fs.RO.mass_transfer_phase_comp[0, "Liq", "NaCl"], 1e3)
    iscale.set_scaling_factor(
        m.fs.RO.feed_side.mass_transfer_term[0, "Liq", "NaCl"], 1e3
    )

    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_vol,  # feed volumetric flow rate [m3/s]
            ("conc_mass_phase_comp", ("Liq", "NaCl")): salt_mass_conc,
        },  # feed NaCl mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P2.efficiency_pump.fix(0.80)  # pump efficiency [-]
    # num_time_points = 2
    # time_set = np.linspace(0, 1400, num_time_points+1)
    m.fs.P1.control_volume.properties_out[0].pressure.fix(275 * pyunits.psi)
    # for i in range(len(time_set+1)):
    #     m.fs.P1.control_volume.properties_out[time_set[i]].pressure.fix(
    #         (250 + 1400/num_time_points * 1/7 * i) * 6895
    #     )  # feed pressure (Pa)

    m.fs.RO.area.fix(7.2)  # membrane area (m^2)
    m.fs.RO.A_comp.fix(1.2286540698976616e-11)  # membrane water permeability (m/Pa/s)
    m.fs.RO.B_comp.fix(5.613990582160079e-08)  # membrane salt permeability (m/s)
    m.fs.RO.permeate.pressure[:].fix(101325)  # permeate pressure (Pa)

    m.fs.RO.feed_side.channel_height.fix(0.0008636)
    m.fs.RO.feed_side.spacer_porosity.fix(0.70818)  # 72.9%
    m.fs.RO.length.fix(1.016 - 2 * 0.0267)  # m

    # m.fs.RO.feed_side.material_accumulation[:, :, :].value = 0
    # m.fs.RO.feed_side.material_accumulation[0, :, :].fix(0)

    assert not hasattr(m.fs.RO.feed_side, "energy_accumulation")
    # initialize RO
    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"] = value(
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
    )
    m.fs.RO.feed_side.properties_in[0].temperature = value(
        m.fs.feed.properties[0].temperature
    )
    m.fs.RO.feed_side.properties_in[0].pressure = value(
        m.fs.P1.control_volume.properties_out[0].pressure
    )
    # m.fs.RO.initialize(optarg=solver.options)
    m.fs.P2.efficiency_pump.fix(0.80)


    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )

def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results

def solve_dynamic(m):
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
        skip_initial=False,
        initial_solver="ipopt",
        initial_solver_options={
            "constr_viol_tol": 1e-8,
            "nlp_scaling_method": "user-scaling",
            "linear_solver": "ma27",
            "OF_ma57_automatic_scaling": "yes",
            "max_iter": 300,
            "tol": 1e-8,
            "halt_on_ampl_error": "no",
        },
    )
    for result in results.results:
        assert_optimal_termination(result)

def initialize_mixer(m, guess=True):
    if guess:
        recovery_guess = 0.9
        m.fs.mixer.recycle_state[0.0].flow_mass_phase_comp["Liq", "H2O"].fix(
            (1 - recovery_guess) * 0.994
        )
        m.fs.mixer.recycle_state[0.0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            0.006 / (1 - recovery_guess)
        )
    m.fs.mixer.initialize()
    m.fs.mixer.recycle_state[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.mixer.recycle_state[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

def do_forward_initialization_pass(m, pass_num=1):
    # initialize unit by unit
    m.fs.feed.properties[0].pressure_osm_phase
    m.fs.feed.initialize()
    m.fs.feed.properties[0].dens_mass_phase.display()
    m.fs.feed.properties[0].pressure_osm_phase.display()
    over_press_factor = 1.15
    pressure_guess = value(m.fs.feed.properties[0].pressure_osm_phase["Liq"]) * over_press_factor

    propagate_state(m.fs.feed_to_P1)

    m.fs.P1.control_volume.properties_out[0].pressure.fix(pressure_guess)
    m.fs.P1.initialize()
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.display()
    propagate_state(m.fs.P1_to_mixer)

    m.fs.P2.initialize()
    propagate_state(m.fs.P2_to_mixer)

    propagate_state(m.fs.RO_to_P2)

    if pass_num > 0:
        m.fs.mixer.initialize()
    else:
        initialize_mixer(m)
    propagate_state(m.fs.mixer_to_RO)

    m.fs.RO.inlet.pressure.fix(2e6)
    m.fs.RO.initialize()
    m.fs.RO.inlet.pressure.unfix()
    propagate_state(m.fs.RO_to_P2)
    propagate_state(m.fs.RO_to_prod)

    m.fs.product.initialize()

def initialize_system(m):
    # initialize unit by unit
    for idx, init_pass in enumerate(range(1)):
        print(f"\n\nINITIALIZATION PASS {idx+1}\n\n")
        do_forward_initialization_pass(m, pass_num=idx)
        print_results(m)

if __name__ == "__main__":
    m = main()